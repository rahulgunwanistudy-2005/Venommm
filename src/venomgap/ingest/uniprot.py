"""UniProt/ToxProt ingest: mature toxin sequences per species per toxin family.

Sequences are the only route to cross-reactivity that does not require new laboratory work. UniProt
annotates every reviewed toxin entry with a `protein_families` string ("Three-finger toxin family,
Short-chain subfamily, ..."), which is what maps an entry onto the VenomGap family vocabulary.

Everything fetched is cached to `data/raw/uniprot/` as JSON. With `VENOMGAP_OFFLINE=1` the module
never touches the network and raises `MissingCacheError` if a species is not cached, so CI and the
offline demo are guaranteed to run on exactly the sequences that were committed.
"""

from __future__ import annotations

import json
import logging
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any
from urllib.parse import urlencode

import requests

from venomgap.config import OFFLINE, SEQUENCE_CACHE
from venomgap.errors import MissingCacheError, VenomGapError

logger = logging.getLogger(__name__)

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
REQUEST_TIMEOUT_S = 60
MAX_RETRIES = 3
RETRY_BACKOFF_S = 3.0
PAGE_SIZE = 500
MAX_PAGES = 12

# UniProt keyword KW-0800 is "Toxin", the annotation the ToxProt programme applies to venom
# components. Without it an organism query returns the whole proteome, and for a species with a
# sequenced genome that is thousands of housekeeping proteins for a handful of toxins.
TOXIN_KEYWORD = "KW-0800"

# Under-sequenced species (Naja sagittifera, Bungarus caeruleus) have too few Toxin-keyword entries
# to estimate identity from. For those the query is broadened to the whole organism and the family
# classifier does the filtering instead. The broadening is recorded in the cache file so that it is
# visible in SOURCES.md rather than being an invisible change of method.
MIN_SEQUENCES_BEFORE_BROADENING = 8

FIELDS = ("accession", "id", "protein_name", "protein_families", "organism_name", "organism_id",
          "length", "sequence", "reviewed", "keywordid")

# UniProt's `protein_families` annotation mapped onto the VenomGap vocabulary. Order matters:
# the first pattern that matches wins, so more specific families are listed before broader ones.
# "Venom Kunitz-type family" must be tested before any generic protease-inhibitor phrasing, and
# "Disintegrin" before "metalloproteinase" because PII/PIII SVMPs carry disintegrin-like domains.
FAMILY_PATTERNS: tuple[tuple[str, str], ...] = (
    ("three-finger toxin", "3FTx"),
    ("snake three-finger", "3FTx"),
    ("phospholipase a2", "PLA2"),
    ("venom metalloproteinase", "SVMP"),
    ("reprolysin", "SVMP"),
    ("peptidase m12b", "SVMP"),
    ("disintegrin", "DIS"),
    ("peptidase s1", "SVSP"),
    ("venom serine protease", "SVSP"),
    ("l-amino-acid oxidase", "LAAO"),
    ("l-amino acid oxidase", "LAAO"),
    ("flavin amine oxidase", "LAAO"),
    ("crisp", "CRISP"),
    ("cysteine-rich secretory protein", "CRISP"),
    ("snaclec", "CTL"),
    ("c-type lectin", "CTL"),
    ("venom kunitz-type", "KSPI"),
    ("bpti/kunitz", "KSPI"),
    ("kunitz", "KSPI"),
    ("natriuretic peptide", "NP"),
    ("pdgf/vegf", "VEGF"),
    ("vegf", "VEGF"),
    ("5'-nucleotidase", "NUC"),
    ("nucleotidase", "NUC"),
    ("phosphodiesterase", "PDE"),
    ("nucleotide pyrophosphatase", "PDE"),
    ("type-b carboxylesterase", "AChE"),
    ("cholinesterase", "AChE"),
    ("complement c3", "CVF"),
)


@dataclass(frozen=True, slots=True)
class ToxinSequence:
    """One mature toxin sequence assigned to a VenomGap family."""

    accession: str
    species: str
    organism_id: int
    family: str
    protein_name: str
    sequence: str
    reviewed: bool
    uniprot_families: str

    @property
    def length(self) -> int:
        return len(self.sequence)


def classify_family(uniprot_families: str, protein_name: str = "") -> str | None:
    """Map a UniProt `protein_families` annotation onto the VenomGap vocabulary.

    Returns None when the entry belongs to no family in the vocabulary. Unclassified entries are
    counted and reported rather than silently dropped -- see `fetch_species_sequences`.
    """
    haystack = f"{uniprot_families} {protein_name}".lower()
    for pattern, family in FAMILY_PATTERNS:
        if pattern in haystack:
            return family
    return None


def _cache_path(species: str, reviewed_only: bool) -> Path:
    slug = species.replace(" ", "_").replace("/", "_")
    suffix = "reviewed" if reviewed_only else "all"
    return SEQUENCE_CACHE / f"{slug}__{suffix}.json"


def _request_tsv(query: str) -> list[dict[str, str]]:
    """Every page of a UniProt TSV query, following the cursor in the Link header."""
    params = {
        "query": query,
        "format": "tsv",
        "fields": ",".join(FIELDS),
        "size": str(PAGE_SIZE),
    }
    url: str | None = f"{UNIPROT_SEARCH}?{urlencode(params)}"
    rows: list[dict[str, str]] = []
    pages = 0
    while url and pages < MAX_PAGES:
        response = _get_with_retry(url)
        lines = response.text.splitlines()
        if lines:
            header = lines[0].split("\t")
            rows.extend(
                dict(zip(header, line.split("\t"), strict=False)) for line in lines[1:]
            )
        url = response.links.get("next", {}).get("url")
        pages += 1
    if url:
        logger.warning("stopped after %d pages for query %r; results may be truncated",
                       MAX_PAGES, query)
    return rows


def _get_with_retry(url: str) -> requests.Response:
    last_error: Exception | None = None
    for attempt in range(MAX_RETRIES):
        try:
            response = requests.get(url, timeout=REQUEST_TIMEOUT_S)
            response.raise_for_status()
        except requests.RequestException as exc:
            last_error = exc
            logger.warning(
                "UniProt request failed (attempt %d/%d): %s", attempt + 1, MAX_RETRIES, exc
            )
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_BACKOFF_S * (attempt + 1))
        else:
            return response
    raise VenomGapError(f"UniProt request failed after {MAX_RETRIES} attempts: {last_error}")


def fetch_species_sequences(
    species: str,
    reviewed_only: bool = False,
    refresh: bool = False,
) -> list[ToxinSequence]:
    """All vocabulary-classifiable toxin sequences for one species, cached on disk.

    `reviewed_only=False` includes TrEMBL entries, which matters for under-studied species: krait
    and *Naja sagittifera* have only a handful of reviewed entries each, and a cross-reactivity
    estimate from three sequences is not one worth reporting.
    """
    path = _cache_path(species, reviewed_only)
    if path.exists() and not refresh:
        payload = json.loads(path.read_text())
        return [ToxinSequence(**record) for record in payload["sequences"]]

    if OFFLINE:
        raise MissingCacheError(
            f"offline mode is set and {path} is not cached; "
            f"run `python -m venomgap.cli fetch-sequences` with network access first"
        )

    clause = f'organism_name:"{species}" AND keyword:{TOXIN_KEYWORD}'
    query = f"{clause} AND reviewed:true" if reviewed_only else clause
    rows = _request_tsv(query)
    broadened = False
    if _count_classifiable(rows) < MIN_SEQUENCES_BEFORE_BROADENING:
        broad_clause = f'organism_name:"{species}"'
        query = f"{broad_clause} AND reviewed:true" if reviewed_only else broad_clause
        logger.info(
            "%s has too few Toxin-keyword entries; broadening to the whole organism query", species
        )
        rows = _request_tsv(query)
        broadened = True

    sequences: list[ToxinSequence] = []
    unclassified = 0
    for row in rows:
        seq = (row.get("Sequence") or "").strip()
        if not seq:
            continue
        family = classify_family(
            row.get("Protein families", ""), row.get("Protein names", "")
        )
        if family is None:
            unclassified += 1
            continue
        try:
            organism_id = int(row.get("Organism (ID)", "0") or 0)
        except ValueError:
            organism_id = 0
        sequences.append(
            ToxinSequence(
                accession=row.get("Entry", ""),
                species=species,
                organism_id=organism_id,
                family=family,
                protein_name=row.get("Protein names", ""),
                sequence=seq,
                reviewed=(row.get("Reviewed", "") == "reviewed"),
                uniprot_families=row.get("Protein families", ""),
            )
        )

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {
                "species": species,
                "query": query,
                "broadened_query": broadened,
                "rows_returned": len(rows),
                "sequences_kept": len(sequences),
                "unclassified_dropped": unclassified,
                "sequences": [asdict(s) for s in sequences],
            },
            indent=1,
        )
    )
    logger.info(
        "%s: %d entries, %d classified into %d families, %d unclassified",
        species,
        len(rows),
        len(sequences),
        len({s.family for s in sequences}),
        unclassified,
    )
    return sequences


def _count_classifiable(rows: list[dict[str, str]]) -> int:
    return sum(
        1
        for row in rows
        if (row.get("Sequence") or "").strip()
        and classify_family(row.get("Protein families", ""), row.get("Protein names", ""))
    )


def group_by_family(sequences: list[ToxinSequence]) -> dict[str, list[ToxinSequence]]:
    out: dict[str, list[ToxinSequence]] = {}
    for s in sequences:
        out.setdefault(s.family, []).append(s)
    return out


def cache_summary() -> dict[str, dict[str, Any]]:
    """What is on disk, for the reproducibility section of SOURCES.md."""
    summary: dict[str, dict[str, Any]] = {}
    if not SEQUENCE_CACHE.exists():
        return summary
    for path in sorted(SEQUENCE_CACHE.glob("*.json")):
        payload = json.loads(path.read_text())
        families: dict[str, int] = {}
        for record in payload["sequences"]:
            families[record["family"]] = families.get(record["family"], 0) + 1
        summary[payload["species"]] = {
            "query": payload["query"],
            "sequences_kept": payload["sequences_kept"],
            "broadened_query": payload.get("broadened_query", False),
            "unclassified_dropped": payload.get("unclassified_dropped", 0),
            "families": dict(sorted(families.items())),
        }
    return summary
