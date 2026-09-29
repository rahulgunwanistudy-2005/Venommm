"""Cross-reactivity from sequence identity.

Antibodies raised against toxin family `f` in population `q` recognise family `f` in population `p`
only to the extent that the epitopes are conserved. The identity between the two populations'
family-`f` sequences is the observable proxy:

    x_f(q -> p) = clip( (I_f(q, p) - theta_f) / (1 - theta_f), 0, 1 )

`I_f(q, p)` is the mean pairwise identity between the mature sequences of family `f` in the two
populations, from global Needleman-Wunsch alignment under BLOSUM62. `theta_f` is the family-specific
recognition floor, fitted on the calibration split only (see `model/calibrate.py`).

Identity is computed at species level, because proteomes and sequence records are indexed by
different things: a proteome is a locality, a UniProt record is a species. A population therefore
inherits its species' sequences, and every value so derived is flagged `sequence_imputed` so that
the assumption is visible in every downstream result rather than buried here.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations
from pathlib import Path

import numpy as np
from Bio import Align
from Bio.Align import substitution_matrices
from numpy.typing import NDArray

from venomgap.config import DERIVED_DIR, FAMILIES, FAMILY_INDEX, THETA_BOUNDS
from venomgap.errors import SequenceUnavailableError
from venomgap.ingest.uniprot import ToxinSequence, group_by_family

logger = logging.getLogger(__name__)

IDENTITY_CACHE_JSON: Path = DERIVED_DIR / "family_identity.json"

# Above this many sequences per family per species, sub-sample deterministically: mean pairwise
# identity converges quickly and the alignment count grows quadratically.
MAX_SEQUENCES_PER_FAMILY = 12

# A family with fewer than this many sequences on one side is not a usable identity estimate.
MIN_SEQUENCES_FOR_IDENTITY = 2


@lru_cache(maxsize=1)
def _aligner() -> Align.PairwiseAligner:
    """Global aligner, BLOSUM62, gap open -11 / extend -1, as fixed in the pre-registration."""
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11.0
    aligner.extend_gap_score = -1.0
    return aligner


def pairwise_identity(seq_a: str, seq_b: str) -> float:
    """Fraction of aligned columns that are identical, over the alignment length.

    Dividing by alignment length rather than by the shorter sequence means an insertion counts
    against identity. That is the conservative choice: it cannot inflate cross-reactivity.
    """
    a = _clean(seq_a)
    b = _clean(seq_b)
    if not a or not b:
        return 0.0
    if a == b:
        return 1.0
    alignment = _aligner().align(a, b)[0]
    target, query = alignment[0], alignment[1]
    matches = sum(1 for x, y in zip(target, query, strict=True) if x == y and x != "-")
    return matches / len(target) if target else 0.0


def _clean(seq: str) -> str:
    """Strip whitespace and any residue outside the BLOSUM62 alphabet."""
    allowed = set("ACDEFGHIKLMNPQRSTVWYBZX*")
    return "".join(ch for ch in seq.strip().upper() if ch in allowed)


def _subsample(sequences: list[ToxinSequence]) -> list[ToxinSequence]:
    """Deterministic sub-sample: reviewed entries first, then longest, then accession order."""
    ordered = sorted(sequences, key=lambda s: (not s.reviewed, -s.length, s.accession))
    return ordered[:MAX_SEQUENCES_PER_FAMILY]


def mean_between_group_identity(
    group_a: list[ToxinSequence], group_b: list[ToxinSequence]
) -> float | None:
    """Mean identity over all cross pairs, or None if either side is too thin to estimate."""
    a = _subsample(group_a)
    b = _subsample(group_b)
    if len(a) < MIN_SEQUENCES_FOR_IDENTITY or len(b) < MIN_SEQUENCES_FOR_IDENTITY:
        return None
    values = [pairwise_identity(x.sequence, y.sequence) for x in a for y in b]
    return float(np.mean(values)) if values else None


def mean_within_group_identity(group: list[ToxinSequence]) -> float | None:
    """Mean identity among a single group's own sequences: the within-species reference level."""
    members = _subsample(group)
    if len(members) < MIN_SEQUENCES_FOR_IDENTITY:
        return None
    values = [pairwise_identity(x.sequence, y.sequence) for x, y in combinations(members, 2)]
    return float(np.mean(values)) if values else None


@dataclass(frozen=True, slots=True)
class IdentityTable:
    """Mean family-level sequence identity between every pair of species with sequence data.

    `values[(species_a, species_b, family)]` is symmetric. Missing keys mean the pair could not be
    estimated for that family, which propagates to a `sequence_imputed` flag downstream rather than
    to a silent zero.
    """

    values: dict[tuple[str, str, str], float]
    species: tuple[str, ...]
    within_species: dict[tuple[str, str], float]

    def get(self, species_a: str, species_b: str, family: str) -> float | None:
        if species_a == species_b:
            return 1.0
        key = (species_a, species_b, family) if species_a <= species_b else (
            species_b, species_a, family
        )
        return self.values.get(key)

    def to_json(self, path: Path = IDENTITY_CACHE_JSON) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            json.dumps(
                {
                    "species": list(self.species),
                    "within_species": {f"{k[0]}|{k[1]}": v for k, v in self.within_species.items()},
                    "between_species": {
                        f"{k[0]}|{k[1]}|{k[2]}": v for k, v in sorted(self.values.items())
                    },
                },
                indent=1,
            )
        )

    @classmethod
    def from_json(cls, path: Path = IDENTITY_CACHE_JSON) -> IdentityTable:
        if not path.exists():
            raise SequenceUnavailableError(
                f"{path} does not exist; run `python -m venomgap.cli build-crossreact` first"
            )
        payload = json.loads(path.read_text())
        values: dict[tuple[str, str, str], float] = {}
        for key, value in payload["between_species"].items():
            a, b, family = key.split("|")
            values[(a, b, family)] = float(value)
        within: dict[tuple[str, str], float] = {}
        for key, value in payload["within_species"].items():
            species, family = key.split("|")
            within[(species, family)] = float(value)
        return cls(
            values=values, species=tuple(payload["species"]), within_species=within
        )


def build_identity_table(sequences_by_species: dict[str, list[ToxinSequence]]) -> IdentityTable:
    """Compute mean family identity for every species pair. This is the expensive step; it caches."""
    grouped = {
        species: group_by_family(seqs) for species, seqs in sequences_by_species.items()
    }
    species_list = tuple(sorted(grouped))

    within: dict[tuple[str, str], float] = {}
    for species in species_list:
        for family, group in grouped[species].items():
            value = mean_within_group_identity(group)
            if value is not None:
                within[(species, family)] = value

    values: dict[tuple[str, str, str], float] = {}
    for species_a, species_b in combinations(species_list, 2):
        shared = set(grouped[species_a]) & set(grouped[species_b])
        for family in shared:
            value = mean_between_group_identity(
                grouped[species_a][family], grouped[species_b][family]
            )
            if value is not None:
                values[(species_a, species_b, family)] = value
    logger.info(
        "identity table: %d species, %d within-species values, %d between-species values",
        len(species_list),
        len(within),
        len(values),
    )
    return IdentityTable(values=values, species=species_list, within_species=within)


def crossreactivity(identity: float, theta: float) -> float:
    """clip((I - theta) / (1 - theta), 0, 1).

    Linear above the floor and hard zero below it. The floor is what encodes the observation that
    conformational epitopes of 15-20 residues lose recognition sharply rather than gradually.
    """
    if theta >= 1.0:
        return 0.0
    return float(np.clip((identity - theta) / (1.0 - theta), 0.0, 1.0))


def crossreact_matrix(
    immunogen_species: list[str],
    target_species: str,
    table: IdentityTable,
    theta: NDArray[np.float64],
    genus_of: dict[str, str] | None = None,
) -> tuple[NDArray[np.float64], bool]:
    """x_f(q -> p) for every immunogen q and family f, plus whether anything was imputed.

    Fallback order when a species pair has no identity estimate for a family:
      1. the same genus -> that genus's mean identity across the families it does have,
      2. otherwise the global mean between-genus identity for that family,
      3. otherwise the global mean identity overall.
    Every fallback sets the imputed flag, which becomes a `sequence_imputed` entry on the result.
    """
    n = len(immunogen_species)
    out = np.zeros((n, len(FAMILIES)), dtype=np.float64)
    imputed = False

    genus_of = genus_of or {}
    global_by_family = _global_family_means(table)
    global_mean = float(np.mean(list(global_by_family.values()))) if global_by_family else 0.7

    for q_index, species in enumerate(immunogen_species):
        same_genus = (
            genus_of.get(species) is not None
            and genus_of.get(species) == genus_of.get(target_species)
        )
        for family in FAMILIES:
            f_index = FAMILY_INDEX[family]
            identity = table.get(species, target_species, family)
            if identity is None:
                imputed = True
                identity = _fallback_identity(
                    species, target_species, family, table, global_by_family, global_mean, same_genus
                )
            out[q_index, f_index] = crossreactivity(identity, float(theta[f_index]))
    return out, imputed


def _fallback_identity(
    species_a: str,
    species_b: str,
    family: str,
    table: IdentityTable,
    global_by_family: dict[str, float],
    global_mean: float,
    same_genus: bool,
) -> float:
    """Best available identity when the exact (pair, family) estimate is missing."""
    if species_a == species_b:
        return 1.0
    pair_values = [
        value
        for (a, b, _), value in table.values.items()
        if {a, b} == {species_a, species_b}
    ]
    if pair_values:
        return float(np.mean(pair_values))
    if same_genus:
        # Congeners: use the mean identity observed between congeners for this family, if any.
        congeneric = [
            value
            for (a, b, fam), value in table.values.items()
            if fam == family and a.split()[0] == b.split()[0]
        ]
        if congeneric:
            return float(np.mean(congeneric))
    return global_by_family.get(family, global_mean)


def _global_family_means(table: IdentityTable) -> dict[str, float]:
    buckets: dict[str, list[float]] = {}
    for (_, _, family), value in table.values.items():
        buckets.setdefault(family, []).append(value)
    return {family: float(np.mean(values)) for family, values in buckets.items()}


def theta_vector_from_scalar(theta: float) -> NDArray[np.float64]:
    """A flat theta across families, used as the calibration starting point."""
    low, high = THETA_BOUNDS
    return np.full(len(FAMILIES), float(np.clip(theta, low, high)), dtype=np.float64)
