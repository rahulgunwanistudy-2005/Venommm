"""Exact solution by enumeration on a reduced instance, to measure the heuristic's gap.

The full problem is not enumerable, and because the objective is non-submodular there is no
approximation guarantee to fall back on. So the gap is *measured*: restrict to a small candidate set
and a small k, enumerate every subset, and compare.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from itertools import combinations

import numpy as np

from venomgap.config import EXACT_INSTANCE_MAX_K, EXACT_INSTANCE_POPULATIONS
from venomgap.errors import OptimisationError
from venomgap.model.assemble import ModelAssembly
from venomgap.optimize.greedy import candidate_populations, evaluate_subset
from venomgap.optimize.objective import NationalObjective

logger = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class ExactResult:
    k: int
    sites: tuple[str, ...]
    coverage: float
    subsets_evaluated: int
    candidates: tuple[str, ...]


def reduced_instance(
    assembly: ModelAssembly, size: int = EXACT_INSTANCE_POPULATIONS
) -> tuple[str, ...]:
    """A deterministic reduced candidate set: the first `size` Indian populations by id.

    Deterministic and boring on purpose. Picking the reduced instance by any quality criterion would
    bias the measured gap in the heuristic's favour.
    """
    pool = candidate_populations(assembly, indian_only=True)
    return pool[:size]


def exact_best(
    objective: NationalObjective,
    assembly: ModelAssembly,
    k: int,
    candidates: tuple[str, ...] | None = None,
) -> ExactResult:
    """Enumerate every subset of size at most k. Feasible only on the reduced instance."""
    pool = candidates or reduced_instance(assembly)
    if k > EXACT_INSTANCE_MAX_K:
        raise OptimisationError(
            f"exact enumeration is limited to k <= {EXACT_INSTANCE_MAX_K}; asked for {k}"
        )
    best: tuple[float, tuple[str, ...]] | None = None
    evaluated = 0
    # Every size up to k, because a smaller set can beat a larger one under a non-monotone
    # objective.
    for size in range(1, k + 1):
        for subset in combinations(pool, size):
            value = evaluate_subset(objective, assembly, subset)
            evaluated += 1
            if best is None or value > best[0]:
                best = (value, subset)
    if best is None:
        raise OptimisationError("exact enumeration found no feasible subset")
    logger.info(
        "exact: k <= %d over %d candidates, %d subsets evaluated, best %.5f with %s",
        k, len(pool), evaluated, best[0], best[1],
    )
    return ExactResult(
        k=k,
        sites=best[1],
        coverage=best[0],
        subsets_evaluated=evaluated,
        candidates=tuple(pool),
    )


def optimality_gap(heuristic_coverage: float, exact_coverage: float) -> float:
    """Relative shortfall of the heuristic against the exact optimum, on the same instance."""
    if exact_coverage <= 0.0:
        return float("nan")
    return float(max(0.0, (exact_coverage - heuristic_coverage) / exact_coverage))


def measure_gap(
    objective: NationalObjective,
    assembly: ModelAssembly,
    k: int = EXACT_INSTANCE_MAX_K,
) -> dict[str, object]:
    """Run both solvers on the same reduced instance and report the gap."""
    from venomgap.optimize.local_search import greedy_plus_local

    pool = reduced_instance(assembly)
    exact = exact_best(objective, assembly, k, pool)
    heuristic = greedy_plus_local(objective, assembly, k, pool)
    gap = optimality_gap(heuristic.coverage, exact.coverage)
    return {
        "k": k,
        "candidates": list(pool),
        "exact_sites": list(exact.sites),
        "exact_coverage": exact.coverage,
        "subsets_evaluated": exact.subsets_evaluated,
        "heuristic_sites": list(heuristic.sites),
        "heuristic_coverage": heuristic.coverage,
        "optimality_gap": gap,
        "heuristic_matched_exact": bool(np.isclose(heuristic.coverage, exact.coverage, atol=1e-9)),
        "note": (
            "The objective is non-monotone and not submodular, so no (1 - 1/e) guarantee applies. "
            "This gap is measured on a reduced instance, not derived."
        ),
    }
