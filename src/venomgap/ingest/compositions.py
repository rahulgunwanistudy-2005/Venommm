"""Curated composition corpus: CSV of published proteomes -> validated `VenomPopulation` objects.

The CSV stores abundances as **published percentages**, one column per source family label, exactly
as the paper printed them. This module does three things and states each in the row's flags:

1. **Maps source family labels onto the VenomGap vocabulary.** Papers split three-finger toxins into
   neurotoxic and cytotoxic, report snaclecs separately from C-type lectins, and name minor families
   (vespryn, cystatin, NGF, hyaluronidase) that have no vocabulary slot and fold into `other`.
2. **Accounts for the unreported residual.** Where the listed families sum to less than 100% of
   whole venom the shortfall goes to `other`, not redistributed. Redistributing it would invent
   abundance for families the study did not see.
3. **Normalises to the simplex** and flags the row when the correction was material.

Rows in, rows out, and every exclusion are counted. A study is never silently dropped.
"""

from __future__ import annotations

import logging
from datetime import date
from pathlib import Path
from typing import TypedDict, cast

import pandas as pd

from venomgap.config import COMPOSITIONS_CSV
from venomgap.errors import DataValidationError, ProvenanceError
from venomgap.types import ProteomicMethod, Provenance, ResultFlag, VenomPopulation

logger = logging.getLogger(__name__)

# Source family labels as the papers print them, mapped onto the vocabulary. Several source labels
# legitimately collapse onto one vocabulary family and their percentages are summed.
SOURCE_FAMILY_MAP: dict[str, str] = {
    "n_3ftx": "3FTx",
    "c_3ftx": "3FTx",
    "3ftx": "3FTx",
    "pla2": "PLA2",
    "svmp": "SVMP",
    "svsp": "SVSP",
    "laao": "LAAO",
    "crisp": "CRISP",
    "snaclec": "CTL",
    "ctl": "CTL",
    "kspi": "KSPI",
    "kunitz": "KSPI",
    "dis": "DIS",
    "disintegrin": "DIS",
    "np": "NP",
    "vegf": "VEGF",
    "nuc": "NUC",
    "nucleotidase": "NUC",
    "pde": "PDE",
    "ache": "AChE",
    "cvf": "CVF",
    # No vocabulary slot: folded into `other` with the residual.
    "vespryn": "other",
    "cystatin": "other",
    "ngf": "other",
    "hyaluronidase": "other",
    "plb": "other",
    "serpin": "other",
    "calreticulin": "other",
    "pli": "other",
    "other": "other",
}

VALID_FLAGS: frozenset[str] = frozenset(
    (
        "sequence_imputed",
        "composition_imputed",
        "no_nearby_proteome",
        "single_study_basis",
        "partial_table_renormalised",
        "genus_consensus_sequences",
    )
)

REQUIRED_COLUMNS: tuple[str, ...] = (
    "pop_id",
    "species",
    "genus",
    "locality",
    "state",
    "region",
    "country",
    "lat",
    "lon",
    "doi",
    "table",
    "accessed",
    "method",
    "holdout",
    "is_immunogen_source",
)

# A residual larger than this fraction of whole venom means the study characterised materially less
# than all of it, and the row is flagged so the atlas and the paper can say so.
MATERIAL_RESIDUAL = 0.10

# Beyond this the row is not usable as a composition estimate at all.
MAX_RESIDUAL = 0.45


def load_compositions(path: Path = COMPOSITIONS_CSV) -> list[VenomPopulation]:
    """Read, validate and convert the curated corpus. Raises on any contract violation."""
    if not path.exists():
        raise DataValidationError(f"curated composition file not found: {path}")

    frame = pd.read_csv(path)
    missing = [c for c in REQUIRED_COLUMNS if c not in frame.columns]
    if missing:
        raise DataValidationError(f"{path} is missing required columns: {missing}")

    family_columns = [c for c in frame.columns if c.startswith("pct_")]
    if not family_columns:
        raise DataValidationError(f"{path} has no pct_* abundance columns")
    unknown = [
        c for c in family_columns if c.removeprefix("pct_").lower() not in SOURCE_FAMILY_MAP
    ]
    if unknown:
        raise DataValidationError(
            f"{path} has abundance columns with no vocabulary mapping: {unknown}"
        )

    populations: list[VenomPopulation] = []
    for index, row in frame.iterrows():
        populations.append(_row_to_population(row, family_columns, row_number=int(index) + 2))

    if len(populations) != len(frame):
        raise DataValidationError(
            f"{len(frame)} curated rows produced {len(populations)} populations; "
            f"no row may be dropped silently"
        )
    duplicates = _duplicates([p.pop_id for p in populations])
    if duplicates:
        raise DataValidationError(f"duplicate pop_id values in {path}: {duplicates}")

    logger.info(
        "loaded %d populations across %d species (%d holdout, %d immunogen-source)",
        len(populations),
        len({p.species for p in populations}),
        sum(p.holdout for p in populations),
        sum(p.is_immunogen_source for p in populations),
    )
    return populations


def _row_to_population(
    row: pd.Series[object], family_columns: list[str], row_number: int
) -> VenomPopulation:
    """One CSV row -> one validated VenomPopulation, with the residual accounted for."""
    pop_id = str(row["pop_id"]).strip()

    for field in ("doi", "table", "accessed"):
        value = row[field]
        if pd.isna(value) or not str(value).strip():
            raise ProvenanceError(
                f"row {row_number} ({pop_id}) is missing {field}; every curated row must carry "
                f"a DOI or PMC id, a table reference and an access date"
            )

    raw: dict[str, float] = {}
    for column in family_columns:
        value = row[column]
        if pd.isna(value):
            continue
        percentage = float(value)
        if percentage < 0.0:
            raise DataValidationError(
                f"row {row_number} ({pop_id}) has a negative abundance in {column}"
            )
        family = SOURCE_FAMILY_MAP[column.removeprefix("pct_").lower()]
        raw[family] = raw.get(family, 0.0) + percentage

    listed_total = sum(raw.values())
    if listed_total <= 0.0:
        raise DataValidationError(f"row {row_number} ({pop_id}) has no reported abundance")

    flags: list[ResultFlag] = []
    residual = max(0.0, 100.0 - listed_total)
    if residual > MAX_RESIDUAL * 100.0:
        raise DataValidationError(
            f"row {row_number} ({pop_id}) characterises only {listed_total:.1f}% of whole venom; "
            f"the unreported residual exceeds the {MAX_RESIDUAL:.0%} limit and the row cannot be "
            f"used as a composition estimate"
        )
    if residual > MATERIAL_RESIDUAL * 100.0:
        flags.append("partial_table_renormalised")
    if residual > 0.0:
        raw["other"] = raw.get("other", 0.0) + residual

    total = sum(raw.values())
    composition = {f: raw[f] / total for f in raw if raw[f] > 0.0}
    # Force an exact simplex: floating-point division leaves a residue of order 1e-16 per family,
    # and VenomPopulation validates to 1e-6, so absorb any remainder into the largest component.
    composition = _snap_to_simplex(composition)

    if str(row["method"]).strip() == "other":
        flags.append("composition_imputed")
    extra_flags = _optional_text(row.get("flags"))
    for flag in (f.strip() for f in extra_flags.split(";") if f.strip()):
        if flag not in VALID_FLAGS:
            raise DataValidationError(f"row {row_number} ({pop_id}) has unknown flag {flag!r}")
        flags.append(cast("ResultFlag", flag))

    provenance = Provenance(
        doi=str(row["doi"]).strip(),
        table=str(row["table"]).strip(),
        accessed=date.fromisoformat(str(row["accessed"]).strip()),
        method=_validated_method(row["method"], pop_id),
        imputed=bool(flags and "composition_imputed" in flags),
        note=_optional_text(row.get("note")),
    )

    return VenomPopulation(
        pop_id=pop_id,
        species=str(row["species"]).strip(),
        genus=str(row["genus"]).strip(),
        locality=str(row["locality"]).strip(),
        state=str(row["state"]).strip(),
        region=str(row["region"]).strip(),
        country=str(row["country"]).strip(),
        lat=float(row["lat"]),
        lon=float(row["lon"]),
        composition=composition,
        provenance=provenance,
        holdout=_as_bool(row["holdout"]),
        is_immunogen_source=_as_bool(row["is_immunogen_source"]),
        flags=tuple(dict.fromkeys(flags)),
    )


def _snap_to_simplex(composition: dict[str, float]) -> dict[str, float]:
    total = sum(composition.values())
    largest = max(composition, key=lambda k: composition[k])
    out = {k: v / total for k, v in composition.items()}
    out[largest] += 1.0 - sum(out.values())
    return out


VALID_METHODS: frozenset[str] = frozenset(
    ("LC-MS/MS", "RP-HPLC+MS", "transcriptome-informed", "other")
)


def _validated_method(value: object, pop_id: str) -> ProteomicMethod:
    text = str(value).strip()
    if text not in VALID_METHODS:
        raise DataValidationError(f"row {pop_id} has unknown proteomic method {text!r}")
    return cast("ProteomicMethod", text)


def _optional_text(value: object) -> str:
    """Empty CSV cells arrive as NaN, which str() would turn into the literal "nan"."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    text = str(value).strip()
    return "" if text.lower() == "nan" else text


def _as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"true", "1", "yes", "y"}:
        return True
    if text in {"false", "0", "no", "n", "", "nan"}:
        return False
    raise DataValidationError(f"cannot interpret {value!r} as a boolean")


def _duplicates(values: list[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        if value in seen and value not in out:
            out.append(value)
        seen.add(value)
    return out


class CorpusSummary(TypedDict):
    """Counts reported by the M2 gate. Every field is computed, never asserted."""

    populations: int
    species: int
    states: int
    studies: int
    holdout_populations: list[str]
    immunogen_sources: list[str]
    by_species: dict[str, int]
    by_state: dict[str, int]
    flagged_populations: list[str]
    dominant_family: dict[str, str]


def corpus_summary(populations: list[VenomPopulation]) -> CorpusSummary:
    """Counts used by the M2 gate, the README and SOURCES.md. Every number here is computed."""
    by_species: dict[str, int] = {}
    by_state: dict[str, int] = {}
    for p in populations:
        by_species[p.species] = by_species.get(p.species, 0) + 1
        by_state[p.state] = by_state.get(p.state, 0) + 1
    flagged = [p.pop_id for p in populations if p.flags]
    return CorpusSummary(
        populations=len(populations),
        species=len(by_species),
        states=len(by_state),
        studies=len({p.provenance.doi for p in populations}),
        holdout_populations=sorted(p.pop_id for p in populations if p.holdout),
        immunogen_sources=sorted(p.pop_id for p in populations if p.is_immunogen_source),
        by_species=dict(sorted(by_species.items())),
        by_state=dict(sorted(by_state.items())),
        flagged_populations=sorted(flagged),
        dominant_family={
            p.pop_id: max(p.composition, key=lambda k: p.composition[k]) for p in populations
        },
    )
