"""Spatial interpolation of venom composition, and the honest labelling of where it fails.

Proteomes exist for a few dozen localities; India has hundreds of districts. Composition must be
interpolated, and the interpolation must say where it has no business speaking.

For a target point and species, the estimate is a Gaussian-kernel weighted mean of that species'
sampled proteomes:

    a_hat(x) = sum_j k(d(x, x_j)) a_j / sum_j k(d(x, x_j)),    k(r) = exp(-r^2 / 2 ell^2)

`ell` is fitted by **leave-one-population-out cross-validation on the composition vectors
themselves**. No coverage, antivenom or outcome information enters that fit, so the fitted length
scale is a statement about venom biogeography alone: how far, in kilometres, a venom proteome
generalises. That number is reportable on its own.

Uncertainty is the kernel-weighted spread of the contributing proteomes, inflated by distance to the
nearest sample. Past a cut-off the answer is `unknown`, never `covered`. Grey on the map, not green.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from venomgap.config import (
    EARTH_RADIUS_KM,
    ELL_GRID_KM,
    ELL_PLATEAU_TOLERANCE,
    EPS,
    FAMILIES,
    UNKNOWN_CUTOFF_ELL_MULTIPLE,
    UNKNOWN_CUTOFF_MAX_KM,
)
from venomgap.errors import SpatialModelError
from venomgap.types import VenomPopulation

logger = logging.getLogger(__name__)

N_FAMILIES = len(FAMILIES)

# Genus-level fallback needs at least this many same-genus populations to be worth doing.
MIN_GENUS_POPULATIONS = 2


def haversine_km(
    lat1: NDArray[np.float64] | float,
    lon1: NDArray[np.float64] | float,
    lat2: NDArray[np.float64] | float,
    lon2: NDArray[np.float64] | float,
) -> NDArray[np.float64]:
    """Great-circle distance in km. Vectorised over any broadcastable shapes."""
    phi1, phi2 = np.radians(lat1), np.radians(lat2)
    dphi = phi2 - phi1
    dlambda = np.radians(lon2) - np.radians(lon1)
    a = np.sin(dphi / 2.0) ** 2 + np.cos(phi1) * np.cos(phi2) * np.sin(dlambda / 2.0) ** 2
    return np.asarray(2.0 * EARTH_RADIUS_KM * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0))),
    dtype=np.float64)


def gaussian_kernel(distance_km: NDArray[np.float64], ell_km: float) -> NDArray[np.float64]:
    if ell_km <= 0.0:
        raise SpatialModelError(f"length scale must be positive, got {ell_km}")
    return np.asarray(np.exp(-(distance_km**2) / (2.0 * ell_km**2)), dtype=np.float64)


@dataclass(frozen=True, slots=True)
class CompositionEstimate:
    """An interpolated composition with everything needed to judge whether to trust it."""

    composition: NDArray[np.float64]
    uncertainty: float
    nearest_km: float
    contributors: tuple[str, ...]
    status: str  # "estimated" | "unknown"

    @property
    def is_known(self) -> bool:
        return self.status == "estimated"


class SpatialCompositionModel:
    """Kernel interpolation over the curated corpus, with a fitted length scale.

    Constructed from the corpus and an `ell_km`; use `fit_length_scale` to obtain the latter by
    leave-one-population-out CV rather than choosing it.
    """

    def __init__(self, populations: list[VenomPopulation], ell_km: float) -> None:
        if not populations:
            raise SpatialModelError("cannot build a spatial model from an empty corpus")
        self.populations = populations
        self.ell_km = float(ell_km)
        self.cutoff_km = min(
            UNKNOWN_CUTOFF_ELL_MULTIPLE * self.ell_km, UNKNOWN_CUTOFF_MAX_KM
        )
        self._by_species: dict[str, list[VenomPopulation]] = {}
        self._by_genus: dict[str, list[VenomPopulation]] = {}
        for p in populations:
            self._by_species.setdefault(p.species, []).append(p)
            self._by_genus.setdefault(p.genus, []).append(p)
        # Subspecies inherit their nominate species' populations as neighbours, so that
        # "Echis carinatus sochureki" can borrow from "Echis carinatus" and vice versa.
        for species, group in list(self._by_species.items()):
            nominate = " ".join(species.split()[:2])
            if nominate != species:
                self._by_species.setdefault(nominate, [])
                for p in group:
                    if p not in self._by_species[nominate]:
                        self._by_species[nominate].append(p)

    def neighbours(self, species: str) -> list[VenomPopulation]:
        """Populations usable as neighbours for a species, widening to the genus if needed."""
        direct = list(self._by_species.get(species, []))
        nominate = " ".join(species.split()[:2])
        for p in self._by_species.get(nominate, []):
            if p not in direct:
                direct.append(p)
        if len(direct) >= MIN_GENUS_POPULATIONS:
            return direct
        genus = species.split()[0]
        return list(self._by_genus.get(genus, direct))

    def estimate(
        self,
        lat: float,
        lon: float,
        species: str,
        exclude: str | None = None,
    ) -> CompositionEstimate:
        """Interpolated composition at a point for a species.

        `exclude` drops one `pop_id`, which is how leave-one-population-out CV is done without
        rebuilding the model.
        """
        pool = [p for p in self.neighbours(species) if p.pop_id != exclude]
        if not pool:
            return CompositionEstimate(
                composition=np.full(N_FAMILIES, 1.0 / N_FAMILIES),
                uncertainty=1.0,
                nearest_km=float("inf"),
                contributors=(),
                status="unknown",
            )

        lats = np.array([p.lat for p in pool])
        lons = np.array([p.lon for p in pool])
        distances = haversine_km(lat, lon, lats, lons)
        nearest_km = float(distances.min())

        weights = gaussian_kernel(distances, self.ell_km)
        total = float(weights.sum())
        if total <= EPS:
            # Every neighbour is so far away that the kernel underflows. Fall back to the single
            # nearest population so the answer is defined, and mark it unknown.
            index = int(np.argmin(distances))
            return CompositionEstimate(
                composition=pool[index].vector(),
                uncertainty=1.0,
                nearest_km=nearest_km,
                contributors=(pool[index].pop_id,),
                status="unknown",
            )

        weights = weights / total
        matrix = np.vstack([p.vector() for p in pool])
        composition = weights @ matrix

        # Kernel-weighted standard deviation across contributors, summed over families: how much
        # the proteomes that actually informed this point disagree with each other.
        spread = float(np.sqrt(np.sum(weights @ (matrix - composition) ** 2)))
        # Distance inflation: at the cut-off the penalty is 1.0, so uncertainty saturates where the
        # status flips to unknown and the two signals agree by construction.
        inflation = min(1.0, nearest_km / self.cutoff_km) if self.cutoff_km > 0 else 1.0
        uncertainty = float(np.clip(spread + inflation * (1.0 - spread), 0.0, 1.0))

        status = "estimated" if nearest_km <= self.cutoff_km else "unknown"
        contributors = tuple(
            pool[i].pop_id for i in np.argsort(-weights)[:5] if weights[i] > 0.01
        )
        composition = _snap(composition)
        return CompositionEstimate(
            composition=composition,
            uncertainty=uncertainty,
            nearest_km=nearest_km,
            contributors=contributors,
            status=status,
        )


def _snap(vector: NDArray[np.float64]) -> NDArray[np.float64]:
    """Force a non-negative vector onto the simplex."""
    clipped = np.clip(vector, 0.0, None)
    total = clipped.sum()
    if total <= EPS:
        return np.full(N_FAMILIES, 1.0 / N_FAMILIES)
    return np.asarray(clipped / total, dtype=np.float64)


@dataclass(frozen=True, slots=True)
class LengthScaleFit:
    """Result of the leave-one-population-out CV over the kernel length scale."""

    ell_km: float
    rmse: float
    grid: tuple[float, ...]
    rmse_by_ell: tuple[float, ...]
    folds: int
    per_population_error: dict[str, float]
    plateau_km: tuple[float, float]
    at_grid_edge: bool

    @property
    def plateau_is_wide(self) -> bool:
        """True when the CV curve barely distinguishes length scales.

        A wide plateau is a statement about the data, not a failure of the fit: it means that at
        the present sampling density geographic distance carries little predictive power for venom
        composition, and the best available estimator is close to a species-level mean.
        """
        low, high = self.plateau_km
        return high >= 4.0 * low


def fit_length_scale(
    populations: list[VenomPopulation],
    grid: tuple[float, ...] = ELL_GRID_KM,
) -> LengthScaleFit:
    """Fit `ell` by leave-one-population-out CV on composition vectors only.

    For each candidate length scale, every population is predicted from the others of its species
    and the root-mean-square error over families is accumulated. The minimiser is the distance over
    which a venom proteome generalises. No antivenom, coverage or clinical quantity is involved, so
    this fit cannot leak holdout information -- a test asserts that holdout populations may take
    part here for exactly that reason.
    """
    if len(populations) < 3:
        raise SpatialModelError("length-scale CV needs at least three populations")

    rmse_by_ell: list[float] = []
    best: tuple[float, float, dict[str, float]] | None = None

    for ell in grid:
        model = SpatialCompositionModel(populations, ell)
        errors: dict[str, float] = {}
        for p in populations:
            # A population with no other proteome of its species (or genus) cannot be predicted at
            # all; it contributes no fold rather than a fabricated error of zero.
            pool = [q for q in model.neighbours(p.species) if q.pop_id != p.pop_id]
            if not pool:
                continue
            estimate = model.estimate(p.lat, p.lon, p.species, exclude=p.pop_id)
            errors[p.pop_id] = float(
                np.sqrt(np.mean((estimate.composition - p.vector()) ** 2))
            )
        if not errors:
            raise SpatialModelError("no population could be cross-validated")
        rmse = float(np.sqrt(np.mean(np.square(list(errors.values())))))
        rmse_by_ell.append(rmse)
        if best is None or rmse < best[1]:
            best = (ell, rmse, errors)

    if best is None:
        raise SpatialModelError("length-scale CV produced no candidate")
    ell, rmse, errors = best
    inside = [
        candidate
        for candidate, error in zip(grid, rmse_by_ell, strict=True)
        if error <= rmse * (1.0 + ELL_PLATEAU_TOLERANCE)
    ]
    plateau = (min(inside), max(inside))
    at_edge = ell in (grid[0], grid[-1])
    logger.info(
        "fitted spatial length scale: ell = %.0f km (LOO-CV RMSE %.4f over %d folds); "
        "plateau within %.0f%% of the minimum spans %.0f-%.0f km%s",
        ell,
        rmse,
        len(errors),
        ELL_PLATEAU_TOLERANCE * 100.0,
        plateau[0],
        plateau[1],
        " [AT GRID EDGE]" if at_edge else "",
    )
    return LengthScaleFit(
        ell_km=ell,
        rmse=rmse,
        grid=tuple(grid),
        rmse_by_ell=tuple(rmse_by_ell),
        folds=len(errors),
        per_population_error=errors,
        plateau_km=plateau,
        at_grid_edge=at_edge,
    )


def composition_dissimilarity_vs_distance(
    populations: list[VenomPopulation],
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.bool_]]:
    """Pairwise (distance_km, dissimilarity, same_species) for the spatial-generalisation figure.

    Dissimilarity is Euclidean distance between composition vectors -- a descriptive quantity for
    figure 7, deliberately not the coverage metric.
    """
    distances: list[float] = []
    dissimilarities: list[float] = []
    same: list[bool] = []
    for i, p in enumerate(populations):
        for q in populations[i + 1 :]:
            distances.append(float(haversine_km(p.lat, p.lon, q.lat, q.lon)))
            dissimilarities.append(float(np.linalg.norm(p.vector() - q.vector())))
            same.append(p.species == q.species or p.genus == q.genus)
    return (
        np.array(distances, dtype=np.float64),
        np.array(dissimilarities, dtype=np.float64),
        np.array(same, dtype=bool),
    )
