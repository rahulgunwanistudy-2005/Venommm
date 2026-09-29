"""The national objective: burden-weighted mean coverage across Indian districts.

    maximise   sum_d burden_d * C(mixture_d | S, w, D)

A district's venom is a mixture over the species plausibly present there, weighted by their relative
bite contribution. Coverage of that mixture is evaluated per species and then combined, because
neutralisation is per-bite: a patient is bitten by one snake, not by a district's average snake.
So the district's coverage is the bite-share-weighted mean of per-species coverages, not the
coverage of an averaged composition. Averaging the compositions first would let a well-covered
species mask a badly-covered one inside the same district.

**The objective is non-monotone in |S| and it is not submodular**, because adding a venom to the
mixture moves antibody supply away from families already covered. No (1 - 1/e) guarantee is claimed
anywhere in this package, and `optimize/exact.py` measures the real gap instead.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from venomgap.config import (
    BITE_SHARE_CSV,
    DEFAULT_B_MG,
    DEFAULT_D_MG,
    FAMILIES,
    RANGE_SPECIES_CSV_FALLBACK,
)
from venomgap.errors import DataValidationError, OptimisationError
from venomgap.model.assemble import ImmunogenSpec, ModelAssembly
from venomgap.model.coverage import SIGMA_VECTOR, coverage_many
from venomgap.types import FittedParameters

logger = logging.getLogger(__name__)

N_FAMILIES = len(FAMILIES)


@dataclass(frozen=True, slots=True)
class DistrictSpecies:
    """One (district, species) cell: an interpolated composition and its bite share."""

    district_id: str
    species: str
    composition: NDArray[np.float64]
    bite_share: float
    uncertainty: float
    nearest_km: float
    status: str


@dataclass
class NationalObjective:
    """Everything needed to score an immunogen against the whole country, fast.

    Compositions and cross-reactivity depend on the target, not on the mixture weights, so they are
    built once. Scoring a candidate mixture is then a single einsum.
    """

    assembly: ModelAssembly
    fitted: FittedParameters
    district_ids: tuple[str, ...]
    burden: NDArray[np.float64]
    cell_district_index: NDArray[np.int64]
    cell_species: tuple[str, ...]
    cell_compositions: NDArray[np.float64]
    cell_weights: NDArray[np.float64]
    cell_uncertainty: NDArray[np.float64]
    cell_nearest_km: NDArray[np.float64]
    cell_status: tuple[str, ...]
    B: float = DEFAULT_B_MG
    D: float = DEFAULT_D_MG

    def __post_init__(self) -> None:
        if self.burden.shape != (len(self.district_ids),):
            raise OptimisationError("burden vector does not match the district list")
        total = float(self.burden.sum())
        if total <= 0.0:
            raise OptimisationError("total burden is zero")
        self.burden = self.burden / total
        self._crossreact_cache: dict[tuple[str, ...], NDArray[np.float64]] = {}

    @property
    def n_cells(self) -> int:
        return len(self.cell_species)

    def crossreact_for(self, immunogen: ImmunogenSpec) -> NDArray[np.float64]:
        """(n_cells, n_immunogen, n_families) cross-reactivity, cached by immunogen species tuple.

        Cross-reactivity depends only on which *species* are in the mixture, never on their weights,
        so the cache key is the species tuple and the weight optimiser gets it for free.
        """
        key = immunogen.species
        cached = self._crossreact_cache.get(key)
        if cached is None:
            stack, _ = self.assembly.crossreact_stack(
                list(self.cell_species), immunogen, self.fitted.theta_vector()
            )
            self._crossreact_cache[key] = stack
            cached = stack
        return cached

    def cell_coverage(self, immunogen: ImmunogenSpec) -> NDArray[np.float64]:
        """Coverage of every (district, species) cell under this immunogen."""
        result: NDArray[np.float64] = coverage_many(
            self.cell_compositions,
            immunogen.compositions,
            self.crossreact_for(immunogen),
            immunogen.weights,
            self.fitted.kappa_vector(),
            SIGMA_VECTOR,
            B=self.B,
            D=self.D,
        )
        return result

    def district_coverage(self, immunogen: ImmunogenSpec) -> NDArray[np.float64]:
        """Bite-share-weighted coverage per district."""
        cells = self.cell_coverage(immunogen)
        numerator = np.bincount(
            self.cell_district_index,
            weights=cells * self.cell_weights,
            minlength=len(self.district_ids),
        )
        denominator = np.bincount(
            self.cell_district_index,
            weights=self.cell_weights,
            minlength=len(self.district_ids),
        )
        out: NDArray[np.float64] = np.zeros_like(numerator)
        present = denominator > 0.0
        out[present] = numerator[present] / denominator[present]
        return out

    def national_coverage(self, immunogen: ImmunogenSpec) -> float:
        """The objective. Burden-weighted mean district coverage across India."""
        return float(np.dot(self.burden, self.district_coverage(immunogen)))

    def districts_with_no_estimate(self) -> NDArray[np.bool_]:
        """Districts where every contributing cell is `unknown`: grey on the map, never green."""
        known = np.bincount(
            self.cell_district_index,
            weights=np.array([s == "estimated" for s in self.cell_status], dtype=float),
            minlength=len(self.district_ids),
        )
        unknown: NDArray[np.bool_] = known <= 0.0
        return unknown


def build_national_objective(
    assembly: ModelAssembly,
    fitted: FittedParameters,
    B: float = DEFAULT_B_MG,
    D: float = DEFAULT_D_MG,
    uniform_within_state: bool = False,
) -> NationalObjective:
    """Assemble the district x species grid once, from the curated tables and the spatial model."""
    from venomgap.ingest.burden import district_burden, load_districts, load_state_burden
    from venomgap.ingest.gbif import fetch_species_occurrences, range_mask

    districts = district_burden(
        load_districts(), load_state_burden(), uniform_within_state=uniform_within_state
    )
    bite_shares = _load_bite_shares()
    species_list = sorted({p.species for p in assembly.populations if p.country == "India"})

    lats = districts["lat"].to_numpy(dtype=float)
    lons = districts["lon"].to_numpy(dtype=float)
    district_ids = tuple(districts["district_id"].astype(str))
    index_of = {d: i for i, d in enumerate(district_ids)}

    published_ranges = _load_published_ranges()

    cells: list[DistrictSpecies] = []
    for species in species_list:
        mask = _presence_mask(
            species, districts, lats, lons, published_ranges, fetch_species_occurrences, range_mask
        )
        if not mask.any():
            logger.warning("species %s is present in no district; it contributes nothing", species)
            continue
        share = bite_shares.get(species)
        if share is None:
            raise DataValidationError(
                f"no bite share for {species}; every species in the corpus needs a row in "
                f"{BITE_SHARE_CSV}, even if the row records a uniform assumption"
            )
        for position in np.flatnonzero(mask):
            estimate = assembly.spatial.estimate(
                float(lats[position]), float(lons[position]), species
            )
            cells.append(
                DistrictSpecies(
                    district_id=district_ids[position],
                    species=species,
                    composition=estimate.composition,
                    bite_share=share,
                    uncertainty=estimate.uncertainty,
                    nearest_km=estimate.nearest_km,
                    status=estimate.status,
                )
            )

    if not cells:
        raise OptimisationError("no (district, species) cell was built; check ranges and corpus")

    logger.info(
        "national objective: %d districts, %d species, %d cells",
        len(district_ids), len(species_list), len(cells),
    )
    return NationalObjective(
        assembly=assembly,
        fitted=fitted,
        district_ids=district_ids,
        burden=districts["burden_weight"].to_numpy(dtype=float),
        cell_district_index=np.array([index_of[c.district_id] for c in cells], dtype=np.int64),
        cell_species=tuple(c.species for c in cells),
        cell_compositions=np.vstack([c.composition for c in cells]),
        cell_weights=np.array([c.bite_share for c in cells], dtype=float),
        cell_uncertainty=np.array([c.uncertainty for c in cells], dtype=float),
        cell_nearest_km=np.array([c.nearest_km for c in cells], dtype=float),
        cell_status=tuple(c.status for c in cells),
        B=B,
        D=D,
    )


def _presence_mask(  # type: ignore[no-untyped-def]
    species, districts, lats, lons, published_ranges, fetch, mask_fn
) -> NDArray[np.bool_]:
    """Range mask for a species: published description where GBIF is too sparse, else occurrences.

    GBIF occurrences are used only to exclude districts, never to weight a species. Where a
    published range description exists it takes precedence, because a presence-only record set with
    a handful of points would otherwise shrink an island endemic to nothing.
    """
    from venomgap.ingest.gbif import MIN_RECORDS_FOR_RANGE

    states = published_ranges.get(species)
    try:
        occurrences = fetch(species)
    except Exception as exc:
        logger.warning("no occurrence cache for %s (%s); using the published range", species, exc)
        occurrences = []

    if states is not None:
        by_state: NDArray[np.bool_] = np.asarray(
            districts["state"].isin(states).to_numpy()
            | districts["geo_state"].isin(states).to_numpy(),
            dtype=bool,
        )
        if len(occurrences) < MIN_RECORDS_FOR_RANGE:
            return by_state
        # Both sources available: union, so a published range can add districts GBIF missed but
        # cannot remove ones it recorded.
        combined: NDArray[np.bool_] = by_state | mask_fn(lats, lons, occurrences)
        return combined
    from_occurrences: NDArray[np.bool_] = mask_fn(lats, lons, occurrences)
    return from_occurrences


def _load_bite_shares() -> dict[str, float]:
    if not BITE_SHARE_CSV.exists():
        raise DataValidationError(f"bite attribution file not found: {BITE_SHARE_CSV}")
    frame = pd.read_csv(BITE_SHARE_CSV)
    if {"species", "relative_bite_share"} - set(frame.columns):
        raise DataValidationError(f"{BITE_SHARE_CSV} needs species and relative_bite_share columns")
    return dict(zip(frame["species"], frame["relative_bite_share"].astype(float), strict=True))


def _load_published_ranges() -> dict[str, list[str]]:
    if not RANGE_SPECIES_CSV_FALLBACK.exists():
        return {}
    frame = pd.read_csv(RANGE_SPECIES_CSV_FALLBACK)
    return {
        str(row["species"]): [s.strip() for s in str(row["states"]).split(";") if s.strip()]
        for _, row in frame.iterrows()
    }
