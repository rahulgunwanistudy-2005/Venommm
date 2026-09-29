"""Mixture-weight optimisation on the simplex for a fixed immunogen set.

For fixed `S` the objective is smooth in `w`, so SLSQP on the simplex is appropriate. This is where
the dilution effect becomes actionable rather than merely descriptive: given the venoms you have,
how much of each should go into the mixture?
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import minimize

from venomgap.errors import OptimisationError
from venomgap.model.assemble import ImmunogenSpec
from venomgap.optimize.objective import NationalObjective

logger = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class WeightResult:
    weights: NDArray[np.float64]
    coverage: float
    uniform_coverage: float
    improved: bool
    converged: bool


def optimise_weights(
    objective: NationalObjective,
    immunogen: ImmunogenSpec,
    min_weight: float = 0.0,
) -> WeightResult:
    """Maximise national coverage over the mixture weights, holding the venom set fixed."""
    n = immunogen.size
    if n == 1:
        value = objective.national_coverage(immunogen.with_weights(np.ones(1)))
        return WeightResult(
            weights=np.ones(1), coverage=value, uniform_coverage=value,
            improved=False, converged=True,
        )

    uniform = np.full(n, 1.0 / n)
    uniform_coverage = objective.national_coverage(immunogen.with_weights(uniform))

    def negative(w: NDArray[np.float64]) -> float:
        # Project onto the simplex defensively: SLSQP can step marginally outside the constraint set
        # and the coverage model validates its weights strictly.
        clipped = np.clip(w, 0.0, None)
        total = clipped.sum()
        if total <= 0.0:
            return 0.0
        return -objective.national_coverage(immunogen.with_weights(clipped / total))

    constraints = [{"type": "eq", "fun": lambda w: float(np.sum(w) - 1.0)}]
    bounds = [(min_weight, 1.0)] * n
    result = minimize(
        negative, uniform, method="SLSQP", bounds=bounds, constraints=constraints,
        options={"maxiter": 300, "ftol": 1e-10},
    )
    weights = np.clip(np.asarray(result.x, dtype=np.float64), 0.0, None)
    total = float(weights.sum())
    if total <= 0.0:
        raise OptimisationError("weight optimisation collapsed to an all-zero mixture")
    weights = weights / total
    coverage = objective.national_coverage(immunogen.with_weights(weights))

    if coverage < uniform_coverage:
        # Never return a worse answer than the uniform mixture we started from.
        logger.info("weight optimisation did not beat the uniform mixture; keeping uniform")
        return WeightResult(
            weights=uniform, coverage=uniform_coverage, uniform_coverage=uniform_coverage,
            improved=False, converged=bool(result.success),
        )
    return WeightResult(
        weights=weights, coverage=coverage, uniform_coverage=uniform_coverage,
        improved=coverage > uniform_coverage + 1e-9, converged=bool(result.success),
    )
