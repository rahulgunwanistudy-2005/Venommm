"""Greedy site selection and the mixture-size curve that R4 tests.

Greedy is the baseline, not the answer: the objective is non-monotone and not submodular, so greedy
carries no approximation guarantee here. `local_search.py` improves on it and `exact.py` measures
how far both are from optimal on a reduced instance.

`best_subset_per_size` is the R4 machinery. For each mixture size it finds the best subset it can at
that size and records the achievable coverage. The resulting curve is what turns over.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from itertools import combinations

import numpy as np

from venomgap.config import BIG_FOUR_SPECIES, MAX_K
from venomgap.errors import OptimisationError
from venomgap.model.assemble import ModelAssembly
from venomgap.optimize.objective import NationalObjective
from venomgap.types import FittedParameters

logger = logging.getLogger(__name__)

# Above this many candidates an exhaustive sweep at a given size is replaced by greedy growth.
EXHAUSTIVE_CANDIDATE_LIMIT = 14


@dataclass(frozen=True, slots=True)
class GreedyResult:
    sites: tuple[str, ...]
    coverage: float
    trace: tuple[tuple[str, float], ...]


def candidate_populations(assembly: ModelAssembly, indian_only: bool = False) -> tuple[str, ...]:
    """Populations the optimiser may choose as collection sites.

    Holdout populations are eligible: the pre-registration holds out four published *findings*, and
    excluding a real venom population from the candidate pool would make the siting answer worse for
    no scientific reason. What the holdout forbids is fitting a parameter on those studies' results,
    which `calibrate.py` enforces separately.
    """
    return tuple(
        sorted(
            p.pop_id
            for p in assembly.populations
            if (p.country == "India" or not indian_only)
        )
    )


def evaluate_subset(
    objective: NationalObjective,
    assembly: ModelAssembly,
    pop_ids: tuple[str, ...],
    weights: np.ndarray | None = None,
) -> float:
    if not pop_ids:
        return 0.0
    immunogen = assembly.immunogen_from_pop_ids(pop_ids, weights)
    return objective.national_coverage(immunogen)


def greedy_sites(
    objective: NationalObjective,
    assembly: ModelAssembly,
    k: int,
    candidates: tuple[str, ...] | None = None,
    seed_sites: tuple[str, ...] = (),
) -> GreedyResult:
    """Add the single best site at each step, k times.

    Because the objective is non-monotone, a greedy step can *lower* coverage. When the best
    available addition does that, greedy stops early rather than accepting a worse portfolio --
    which is itself the dilution effect showing up inside the solver.
    """
    pool = list(candidates or candidate_populations(assembly))
    chosen = list(seed_sites)
    trace: list[tuple[str, float]] = []
    current = evaluate_subset(objective, assembly, tuple(chosen)) if chosen else 0.0

    for _ in range(k - len(chosen)):
        best: tuple[float, str] | None = None
        for candidate in pool:
            if candidate in chosen:
                continue
            value = evaluate_subset(objective, assembly, (*chosen, candidate))
            if best is None or value > best[0]:
                best = (value, candidate)
        if best is None:
            break
        if chosen and best[0] <= current + 1e-12:
            logger.info(
                "greedy stopped at |S| = %d: the best remaining addition (%s) would not improve "
                "coverage (%.5f -> %.5f). This is the dilution penalty, not a solver failure.",
                len(chosen), best[1], current, best[0],
            )
            break
        current = best[0]
        chosen.append(best[1])
        trace.append((best[1], best[0]))

    return GreedyResult(sites=tuple(chosen), coverage=current, trace=tuple(trace))


def best_subset_at_size(
    objective: NationalObjective,
    assembly: ModelAssembly,
    size: int,
    candidates: tuple[str, ...] | None = None,
) -> tuple[tuple[str, ...], float]:
    """Best subset of exactly `size` populations, exhaustively when that is affordable.

    R4 needs the *best achievable* coverage at each size, not the coverage of whatever greedy
    happened to build, because a curve that turns over only because the solver got worse would prove
    nothing about dilution.
    """
    pool = list(candidates or candidate_populations(assembly))
    if size <= 0 or size > len(pool):
        raise OptimisationError(f"cannot choose {size} of {len(pool)} candidates")

    if len(pool) <= EXHAUSTIVE_CANDIDATE_LIMIT:
        best: tuple[float, tuple[str, ...]] | None = None
        for subset in combinations(pool, size):
            value = evaluate_subset(objective, assembly, subset)
            if best is None or value > best[0]:
                best = (value, subset)
        assert best is not None
        return best[1], best[0]

    # Too many candidates to enumerate: grow greedily to `size` and then improve by single swaps,
    # which is a far better estimate of the achievable optimum than plain greedy.
    from venomgap.optimize.local_search import improve_by_swaps

    grown = _greedy_exactly(objective, assembly, size, tuple(pool))
    improved = improve_by_swaps(objective, assembly, grown, tuple(pool))
    return improved.sites, improved.coverage


def _greedy_exactly(
    objective: NationalObjective,
    assembly: ModelAssembly,
    size: int,
    pool: tuple[str, ...],
) -> tuple[str, ...]:
    """Greedy growth that does not stop early, so a subset of exactly `size` is always returned."""
    chosen: list[str] = []
    for _ in range(size):
        best: tuple[float, str] | None = None
        for candidate in pool:
            if candidate in chosen:
                continue
            value = evaluate_subset(objective, assembly, (*chosen, candidate))
            if best is None or value > best[0]:
                best = (value, candidate)
        if best is None:
            break
        chosen.append(best[1])
    return tuple(chosen)


def best_subset_per_size(
    objective: NationalObjective,
    assembly: ModelAssembly,
    fitted: FittedParameters,  # noqa: ARG001 - uniform signature with the cosine curve builder
    max_size: int = MAX_K,
    candidates: tuple[str, ...] | None = None,
):  # type: ignore[no-untyped-def]
    """The R4 curve: best achievable national coverage at each mixture size 1..max_size.

    The candidate pool is restricted to one population per Big Four species plus the most
    informative alternatives, so that "mixture size" means what it means in an immunisation
    protocol. Without that restriction, size 6 could be six *Naja* populations, which is a different
    question from the one R4 asks.
    """
    from venomgap.validate.retrodiction import MixtureSizeCurve

    pool = list(candidates or _immunogen_candidate_pool(assembly))
    sizes: list[int] = []
    coverages: list[float] = []
    best_sets: list[tuple[str, ...]] = []
    for size in range(1, min(max_size, len(pool)) + 1):
        subset, value = best_subset_at_size(objective, assembly, size, tuple(pool))
        sizes.append(size)
        coverages.append(value)
        best_sets.append(subset)
        logger.info("|S| = %d: best coverage %.5f with %s", size, value, subset)
    return MixtureSizeCurve(
        sizes=tuple(sizes), coverages=tuple(coverages), best_sets=tuple(best_sets)
    )


def _immunogen_candidate_pool(assembly: ModelAssembly) -> tuple[str, ...]:
    """A pool with one representative per medically important species, plus Big Four alternatives.

    Chosen deterministically: for each species, the Indian population geographically closest to the
    national burden centroid, so the pool is not hand-picked.
    """
    from venomgap.ingest.burden import district_burden, load_districts, load_state_burden
    from venomgap.model.spatial import haversine_km

    districts = district_burden(load_districts(), load_state_burden())
    weights = districts["burden_weight"].to_numpy(dtype=float)
    centroid_lat = float(np.dot(weights, districts["lat"].to_numpy(dtype=float)))
    centroid_lon = float(np.dot(weights, districts["lon"].to_numpy(dtype=float)))

    by_species: dict[str, list] = {}
    for population in assembly.populations:
        if population.country != "India":
            continue
        by_species.setdefault(population.species, []).append(population)

    pool: list[str] = []
    for species in sorted(by_species):
        group = by_species[species]
        distances = [
            float(haversine_km(centroid_lat, centroid_lon, p.lat, p.lon)) for p in group
        ]
        pool.append(group[int(np.argmin(distances))].pop_id)

    # Guarantee the Big Four species are represented even if one has no Indian population.
    represented = {assembly.by_id[p].species for p in pool}
    for species in BIG_FOUR_SPECIES:
        if species in represented:
            continue
        fallback = [p for p in assembly.populations if p.species == species]
        if fallback:
            pool.append(fallback[0].pop_id)
    return tuple(pool)
