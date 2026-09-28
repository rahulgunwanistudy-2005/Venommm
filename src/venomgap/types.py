"""Data contracts. Nothing untyped crosses a module boundary.

Every model in this file validates on construction. A composition that does not sum to 1, a
mixture whose weights do not sum to 1, or a provenance record missing a DOI is a construction-time
failure, not a downstream surprise.
"""

from __future__ import annotations

from datetime import date
from typing import Annotated, Literal

import numpy as np
from numpy.typing import NDArray
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from venomgap.config import FAMILIES, FAMILY_INDEX, SIMPLEX_TOL

ProteomicMethod = Literal[
    "LC-MS/MS",
    "RP-HPLC+MS",
    "transcriptome-informed",
    "other",
]

ResultFlag = Literal[
    "sequence_imputed",
    "composition_imputed",
    "no_nearby_proteome",
    "single_study_basis",
    "partial_table_renormalised",
    "genus_consensus_sequences",
]

Fraction = Annotated[float, Field(ge=0.0, le=1.0)]


# --------------------------------------------------------------------------------------
# Provenance
# --------------------------------------------------------------------------------------


class Provenance(BaseModel):
    """Where a single curated number came from. One of these per composition row."""

    model_config = ConfigDict(frozen=True)

    doi: str = Field(min_length=3)
    table: str = Field(min_length=1)
    accessed: date
    method: ProteomicMethod
    imputed: bool = False
    note: str = ""

    @field_validator("doi")
    @classmethod
    def _doi_is_resolvable(cls, v: str) -> str:
        v = v.strip()
        if not (v.startswith("10.") or v.startswith("PMC") or v.startswith("PMID:")):
            raise ValueError(
                f"provenance identifier {v!r} must be a DOI (10.x), a PMC id, or PMID:x"
            )
        return v


# --------------------------------------------------------------------------------------
# Venom populations
# --------------------------------------------------------------------------------------


class VenomPopulation(BaseModel):
    """A geographically localised venom proteome.

    `composition` is the relative abundance of each toxin family, keys a subset of FAMILIES,
    summing to 1 within SIMPLEX_TOL. Absent keys are zero.
    """

    model_config = ConfigDict(frozen=True)

    pop_id: str = Field(min_length=3, pattern=r"^[A-Za-z0-9_]+__[A-Za-z0-9_]+$")
    species: str = Field(min_length=3)
    genus: str = Field(min_length=3)
    locality: str = Field(min_length=2)
    state: str = Field(min_length=2)
    region: str = Field(min_length=2)
    lat: float = Field(ge=-90.0, le=90.0)
    lon: float = Field(ge=-180.0, le=180.0)
    composition: dict[str, float]
    provenance: Provenance
    holdout: bool = False
    is_immunogen_source: bool = False
    flags: tuple[ResultFlag, ...] = ()

    @field_validator("composition")
    @classmethod
    def _valid_simplex(cls, v: dict[str, float]) -> dict[str, float]:
        unknown = set(v) - set(FAMILIES)
        if unknown:
            raise ValueError(f"composition has families outside the vocabulary: {sorted(unknown)}")
        if any(x < 0.0 for x in v.values()):
            raise ValueError("composition has a negative abundance")
        total = sum(v.values())
        if abs(total - 1.0) > SIMPLEX_TOL:
            raise ValueError(f"composition sums to {total:.9f}, expected 1 +/- {SIMPLEX_TOL}")
        return v

    @model_validator(mode="after")
    def _genus_matches_species(self) -> VenomPopulation:
        if not self.species.startswith(self.genus):
            raise ValueError(f"species {self.species!r} does not start with genus {self.genus!r}")
        return self

    def vector(self) -> NDArray[np.float64]:
        """Composition as a dense vector aligned to FAMILIES."""
        out = np.zeros(len(FAMILIES), dtype=np.float64)
        for family, value in self.composition.items():
            out[FAMILY_INDEX[family]] = value
        return out


# --------------------------------------------------------------------------------------
# Immunogen mixtures
# --------------------------------------------------------------------------------------


class Immunogen(BaseModel):
    """An immunising mixture: which populations, and in what proportions."""

    model_config = ConfigDict(frozen=True)

    populations: tuple[str, ...] = Field(min_length=1)
    weights: tuple[float, ...] = Field(min_length=1)
    label: str = ""

    @model_validator(mode="after")
    def _weights_on_simplex(self) -> Immunogen:
        if len(self.populations) != len(self.weights):
            raise ValueError(
                f"{len(self.populations)} populations but {len(self.weights)} weights"
            )
        if len(set(self.populations)) != len(self.populations):
            raise ValueError("immunogen set contains a duplicate population")
        if any(w < 0.0 for w in self.weights):
            raise ValueError("immunogen has a negative mixture weight")
        total = sum(self.weights)
        if abs(total - 1.0) > SIMPLEX_TOL:
            raise ValueError(f"mixture weights sum to {total:.9f}, expected 1")
        return self

    @property
    def size(self) -> int:
        return len(self.populations)

    @classmethod
    def uniform(cls, populations: tuple[str, ...], label: str = "") -> Immunogen:
        n = len(populations)
        if n == 0:
            raise ValueError("cannot build a uniform immunogen from an empty population set")
        return cls(populations=populations, weights=tuple([1.0 / n] * n), label=label)


# --------------------------------------------------------------------------------------
# Fitted parameters
# --------------------------------------------------------------------------------------


class FittedParameters(BaseModel):
    """Output of `model.calibrate`. The only place fitted numbers are allowed to live."""

    model_config = ConfigDict(frozen=True)

    theta: dict[str, float]
    kappa: dict[str, float]
    kappa_scale: float
    ell_km: float | None = None
    calibration_rows: int
    calibration_rmse: float
    calibration_studies: tuple[str, ...] = ()
    excluded_study: str | None = None
    converged: bool = True

    @field_validator("theta", "kappa")
    @classmethod
    def _covers_vocabulary(cls, v: dict[str, float]) -> dict[str, float]:
        missing = set(FAMILIES) - set(v)
        if missing:
            raise ValueError(f"fitted parameter dict is missing families: {sorted(missing)}")
        return v

    def theta_vector(self) -> NDArray[np.float64]:
        return np.array([self.theta[f] for f in FAMILIES], dtype=np.float64)

    def kappa_vector(self) -> NDArray[np.float64]:
        return np.array([self.kappa[f] for f in FAMILIES], dtype=np.float64)


# --------------------------------------------------------------------------------------
# Results
# --------------------------------------------------------------------------------------


class CoverageResult(BaseModel):
    """Per-population coverage under a given immunogen and dose."""

    model_config = ConfigDict(frozen=True)

    pop_id: str
    coverage: Fraction
    deficit: Fraction
    per_family_neutralised: dict[str, float]
    limiting_families: tuple[str, ...]
    uncertainty: float = Field(ge=0.0)
    flags: tuple[ResultFlag, ...] = ()

    @model_validator(mode="after")
    def _deficit_complements_coverage(self) -> CoverageResult:
        if abs(self.coverage + self.deficit - 1.0) > 1e-9:
            raise ValueError(
                f"coverage {self.coverage} and deficit {self.deficit} do not sum to 1"
            )
        return self


class DistrictResult(BaseModel):
    """Per-district coverage, the unit the atlas renders."""

    model_config = ConfigDict(frozen=True)

    district_id: str
    district: str
    state: str
    lat: float
    lon: float
    coverage: Fraction
    deficit: Fraction
    uncertainty: float = Field(ge=0.0)
    status: Literal["estimated", "unknown"]
    burden_weight: float = Field(ge=0.0)
    nearest_proteome_km: float = Field(ge=0.0)
    dominant_species: tuple[str, ...]
    limiting_families: tuple[str, ...]
    flags: tuple[ResultFlag, ...] = ()

    @model_validator(mode="after")
    def _unknown_districts_carry_a_flag(self) -> DistrictResult:
        if self.status == "unknown" and "no_nearby_proteome" not in self.flags:
            raise ValueError(
                "a district with status='unknown' must carry the no_nearby_proteome flag"
            )
        return self


class SitingSolution(BaseModel):
    """A solved collection-site portfolio."""

    model_config = ConfigDict(frozen=True)

    k: int = Field(ge=0)
    sites: tuple[str, ...]
    weights: tuple[float, ...]
    national_coverage: Fraction
    coverage_gain_vs_baseline: float
    districts_moved_out_of_high_deficit: int = Field(ge=0)
    solver: Literal["greedy", "greedy+local", "exact", "baseline"]
    optimality_gap: float | None = None

    @model_validator(mode="after")
    def _sites_match_weights(self) -> SitingSolution:
        if len(self.sites) != len(self.weights):
            raise ValueError("sites and weights differ in length")
        if self.sites and abs(sum(self.weights) - 1.0) > SIMPLEX_TOL:
            raise ValueError("siting solution weights do not sum to 1")
        return self


class CalibrationRecord(BaseModel):
    """One published antivenomics / immunorecognition observation used for fitting."""

    model_config = ConfigDict(frozen=True)

    record_id: str
    study: str
    doi: str
    table: str
    antivenom: str
    pop_id: str
    family: str
    observed_recognition: Fraction
    observable: Literal["immunorecognition_fraction", "neutralisation_fraction"]
    holdout: bool = False
    note: str = ""

    @field_validator("family")
    @classmethod
    def _known_family(cls, v: str) -> str:
        if v not in FAMILIES:
            raise ValueError(f"unknown toxin family {v!r}")
        return v


class RetrodictionCriterion(BaseModel):
    """One pre-registered sub-criterion and its evaluated outcome."""

    model_config = ConfigDict(frozen=True)

    criterion_id: str
    statement: str
    computed: float | None
    threshold: float | None
    comparison: Literal[">=", "<=", "==", "bool", "ill_posed"]
    passed: bool
    detail: str = ""


class RetrodictionOutcome(BaseModel):
    """The full blinded-retrodiction result, written to results/retrodiction.json."""

    model_config = ConfigDict(frozen=True)

    preregistration_commit: str
    data_version: str
    model_variant: str
    criteria: tuple[RetrodictionCriterion, ...]
    targets_passed: int
    targets_total: int
    notes: str = ""
