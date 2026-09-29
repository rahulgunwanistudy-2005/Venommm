"""Greedy + local search: the working solver.

Add/drop/swap moves on the chosen site set, accepting any move that improves the objective. This is
the honest solver for a non-submodular objective: no guarantee, but a measurable gap against exact
enumeration on a reduced instance.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

from venomgap.config import LOCAL_SEARCH_MAX_ITERS
from venomgap.model.assemble import ModelAssembly
from venomgap.optimize.greedy import candidate_populations, evaluate_subset, greedy_sites
from venomgap.optimize.objective import NationalObjective

logger = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class LocalSearchResult:
    sites: tuple[str, ...]
    coverage: float
    iterations: int
    moves: tuple[str, ...]


def improve_by_swaps(
    objective: NationalObjective,
    assembly: ModelAssembly,
    sites: tuple[str, ...],
    candidates: tuple[str, ...] | None = None,
    max_iters: int = LOCAL_SEARCH_MAX_ITERS,
) -> LocalSearchResult:
    """Swap moves only, holding |S| fixed. Used where the mixture size is part of the question."""
    pool = list(candidates or candidate_populations(assembly))
    current = list(sites)
    best = evaluate_subset(objective, assembly, tuple(current))
    moves: list[str] = []

    for iteration in range(max_iters):
        improved = False
        for position, incumbent in enumerate(current):
            for candidate in pool:
                if candidate in current:
                    continue
                trial = list(current)
                trial[position] = candidate
                value = evaluate_subset(objective, assembly, tuple(trial))
                if value > best + 1e-12:
                    best, current, improved = value, trial, True
                    moves.append(f"swap {incumbent} -> {candidate} ({value:.5f})")
                    break
            if improved:
                break
        if not improved:
            return LocalSearchResult(
                sites=tuple(current), coverage=best, iterations=iteration, moves=tuple(moves)
            )
    logger.warning("swap search hit the %d-iteration cap without converging", max_iters)
    return LocalSearchResult(
        sites=tuple(current), coverage=best, iterations=max_iters, moves=tuple(moves)
    )


def greedy_plus_local(
    objective: NationalObjective,
    assembly: ModelAssembly,
    k: int,
    candidates: tuple[str, ...] | None = None,
    max_iters: int = LOCAL_SEARCH_MAX_ITERS,
) -> LocalSearchResult:
    """Greedy, then add / drop / swap until no single move improves the objective.

    Drop moves matter here in a way they would not for a monotone objective: because coverage can
    fall when a venom is added, the best portfolio of size at most k is sometimes smaller than k.
    """
    pool = list(candidates or candidate_populations(assembly))
    start = greedy_sites(objective, assembly, k, tuple(pool))
    current = list(start.sites)
    best = start.coverage
    moves: list[str] = [f"greedy -> {start.sites} ({best:.5f})"]

    for iteration in range(max_iters):
        improved = False

        # swap
        for position, incumbent in enumerate(current):
            for candidate in pool:
                if candidate in current:
                    continue
                trial = list(current)
                trial[position] = candidate
                value = evaluate_subset(objective, assembly, tuple(trial))
                if value > best + 1e-12:
                    best, current, improved = value, trial, True
                    moves.append(f"swap {incumbent} -> {candidate} ({value:.5f})")
                    break
            if improved:
                break

        # add, if there is room
        if not improved and len(current) < k:
            for candidate in pool:
                if candidate in current:
                    continue
                trial = [*current, candidate]
                value = evaluate_subset(objective, assembly, tuple(trial))
                if value > best + 1e-12:
                    best, current, improved = value, trial, True
                    moves.append(f"add {candidate} ({value:.5f})")
                    break

        # drop: only useful because the objective is non-monotone
        if not improved and len(current) > 1:
            for position, incumbent in enumerate(current):
                trial = [s for i, s in enumerate(current) if i != position]
                value = evaluate_subset(objective, assembly, tuple(trial))
                if value > best + 1e-12:
                    best, current, improved = value, trial, True
                    moves.append(f"drop {incumbent} ({value:.5f}) -- dilution")
                    break

        if not improved:
            return LocalSearchResult(
                sites=tuple(current), coverage=best, iterations=iteration, moves=tuple(moves)
            )
    logger.warning("local search hit the %d-iteration cap without converging", max_iters)
    return LocalSearchResult(
        sites=tuple(current), coverage=best, iterations=max_iters, moves=tuple(moves)
    )
