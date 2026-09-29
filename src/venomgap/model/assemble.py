"""Wiring: corpus + sequences + parameters -> `CoverageInputs` the model can evaluate.

Everything upstream produces pieces (a composition, an identity table, a fitted theta); this module
is the single place they are assembled, so that calibration, retrodiction, the optimiser and the API
all evaluate exactly the same model. A second assembly path would be a second model.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

from venomgap.config import (
    BIG_FOUR_SPECIES,
    FAMILIES,
    IRULA_LAT,
    IRULA_LON,
    KAPPA_PRIOR,
    THETA_PRIOR,
)
from venomgap.errors import DataValidationError
from venomgap.model.coverage import SIGMA_VECTOR, CoverageInputs
from venomgap.model.crossreact import IdentityTable, crossreact_matrix
from venomgap.model.spatial import SpatialCompositionModel
from venomgap.types import ResultFlag, VenomPopulation

logger = logging.getLogger(__name__)

N_FAMILIES = len(FAMILIES)
KAPPA_PRIOR_VECTOR: NDArray[np.float64] = np.array(
    [KAPPA_PRIOR[f] for f in FAMILIES], dtype=np.float64
)


@dataclass(frozen=True, slots=True)
class ImmunogenSpec:
    """An immunising mixture described by what goes into it, not by a population id.

    A component is either a curated population (`pop_id` set) or a species interpolated at a
    coordinate (`species` + `lat`/`lon`), which is how the Big Four immunogen is represented: no
    complete proteome exists at the collection locality, so it is estimated there like any district.
    """

    labels: tuple[str, ...]
    species: tuple[str, ...]
    compositions: NDArray[np.float64]
    weights: NDArray[np.float64]
    uncertainties: NDArray[np.float64]
    flags: tuple[ResultFlag, ...] = ()

    def __post_init__(self) -> None:
        n = len(self.labels)
        if not n:
            raise DataValidationError("immunogen spec has no components")
        if len(self.species) != n or self.compositions.shape != (n, N_FAMILIES):
            raise DataValidationError("immunogen spec components are inconsistent")
        if abs(float(self.weights.sum()) - 1.0) > 1e-6:
            raise DataValidationError(
                f"immunogen weights sum to {float(self.weights.sum()):.9f}, expected 1"
            )

    @property
    def size(self) -> int:
        return len(self.labels)

    def with_weights(self, weights: NDArray[np.float64]) -> ImmunogenSpec:
        return ImmunogenSpec(
            labels=self.labels,
            species=self.species,
            compositions=self.compositions,
            weights=np.asarray(weights, dtype=np.float64),
            uncertainties=self.uncertainties,
            flags=self.flags,
        )


@dataclass
class ModelAssembly:
    """Everything needed to evaluate coverage for any target against any immunogen."""

    populations: list[VenomPopulation]
    identity: IdentityTable
    spatial: SpatialCompositionModel
    genus_of: dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.by_id = {p.pop_id: p for p in self.populations}
        if not self.genus_of:
            self.genus_of = {p.species: p.genus for p in self.populations}
            # Subspecies and their nominate species share a genus, which is what lets
            # `Echis carinatus sochureki` count as a congener of `Echis carinatus`.
            for species in list(self.genus_of):
                nominate = " ".join(species.split()[:2])
                self.genus_of.setdefault(nominate, self.genus_of[species])

    # ---------------------------------------------------------------- immunogens

    def immunogen_from_pop_ids(
        self, pop_ids: tuple[str, ...], weights: NDArray[np.float64] | None = None
    ) -> ImmunogenSpec:
        """Immunogen built from curated populations, which is what the optimiser selects."""
        missing = [p for p in pop_ids if p not in self.by_id]
        if missing:
            raise DataValidationError(f"unknown immunogen populations: {missing}")
        n = len(pop_ids)
        chosen = [self.by_id[p] for p in pop_ids]
        flags: list[ResultFlag] = []
        for p in chosen:
            flags.extend(p.flags)
        return ImmunogenSpec(
            labels=pop_ids,
            species=tuple(p.species for p in chosen),
            compositions=np.vstack([p.vector() for p in chosen]),
            weights=np.full(n, 1.0 / n) if weights is None else np.asarray(weights, float),
            uncertainties=np.zeros(n),
            flags=tuple(dict.fromkeys(flags)),
        )

    def big_four_immunogen(
        self,
        weights: NDArray[np.float64] | None = None,
        substitute: dict[str, str] | None = None,
    ) -> ImmunogenSpec:
        """The Indian polyvalent immunogen, interpolated at the Irula collection locality.

        `substitute` swaps one species for another, which is how the R2b test replaces
        *Naja naja* with *Naja kaouthia* while holding everything else fixed.
        """
        substitute = substitute or {}
        species = [substitute.get(s, s) for s in BIG_FOUR_SPECIES]
        compositions = []
        uncertainties = []
        labels = []
        flags: list[ResultFlag] = ["composition_imputed"]
        for name in species:
            estimate = self.spatial.estimate(IRULA_LAT, IRULA_LON, name)
            compositions.append(estimate.composition)
            uncertainties.append(estimate.uncertainty)
            labels.append(f"{name.replace(' ', '_')}__IrulaTamilNadu")
            if not estimate.is_known:
                flags.append("no_nearby_proteome")
        n = len(species)
        return ImmunogenSpec(
            labels=tuple(labels),
            species=tuple(species),
            compositions=np.vstack(compositions),
            weights=np.full(n, 1.0 / n) if weights is None else np.asarray(weights, float),
            uncertainties=np.array(uncertainties),
            flags=tuple(dict.fromkeys(flags)),
        )

    # ---------------------------------------------------------------- coverage inputs

    def kappa(self, kappa_scale: float) -> NDArray[np.float64]:
        """kappa_f = kappa_scale * (M_ab / M_f). One fitted number; the shape is a priori."""
        return kappa_scale * KAPPA_PRIOR_VECTOR

    def inputs_for(
        self,
        target_composition: NDArray[np.float64],
        target_species: str,
        immunogen: ImmunogenSpec,
        theta: NDArray[np.float64],
        kappa_scale: float,
    ) -> tuple[CoverageInputs, bool]:
        """`CoverageInputs` for one target, plus whether any cross-reactivity was imputed."""
        crossreact, imputed = crossreact_matrix(
            list(immunogen.species), target_species, self.identity, theta, self.genus_of
        )
        return (
            CoverageInputs(
                target_composition=target_composition,
                immunogen_compositions=immunogen.compositions,
                crossreact=crossreact,
                kappa=self.kappa(kappa_scale),
                sigma=SIGMA_VECTOR,
            ),
            imputed,
        )

    def inputs_for_population(
        self,
        pop_id: str,
        immunogen: ImmunogenSpec,
        theta: NDArray[np.float64],
        kappa_scale: float,
    ) -> tuple[CoverageInputs, bool]:
        population = self.by_id[pop_id]
        return self.inputs_for(
            population.vector(), population.species, immunogen, theta, kappa_scale
        )

    def crossreact_stack(
        self,
        target_species: list[str],
        immunogen: ImmunogenSpec,
        theta: NDArray[np.float64],
    ) -> tuple[NDArray[np.float64], NDArray[np.bool_]]:
        """(n_targets, n_immunogen, n_families) cross-reactivity for the vectorised hot path.

        Cross-reactivity depends only on the *species* pair, and a national grid has thousands of
        cells but a handful of species, so it is computed once per distinct species and broadcast.
        """
        distinct: dict[str, tuple[NDArray[np.float64], bool]] = {}
        for species in target_species:
            if species not in distinct:
                distinct[species] = crossreact_matrix(
                    list(immunogen.species), species, self.identity, theta, self.genus_of
                )
        stack = np.empty((len(target_species), immunogen.size, N_FAMILIES), dtype=np.float64)
        imputed = np.zeros(len(target_species), dtype=bool)
        for index, species in enumerate(target_species):
            matrix, was_imputed = distinct[species]
            stack[index] = matrix
            imputed[index] = was_imputed
        return stack, imputed


def theta_from_parameters(
    theta_global: float, deviations: dict[str, float]
) -> NDArray[np.float64]:
    """Assemble the theta vector from a global floor plus per-family deviations.

    Families with no calibration support keep the a priori floor exactly, so the fit cannot move a
    recognition threshold it has no evidence about.
    """
    theta = np.full(N_FAMILIES, float(theta_global), dtype=np.float64)
    for family, delta in deviations.items():
        theta[FAMILIES.index(family)] = float(theta_global + delta)
    return theta


def default_theta() -> NDArray[np.float64]:
    return np.full(N_FAMILIES, THETA_PRIOR, dtype=np.float64)
