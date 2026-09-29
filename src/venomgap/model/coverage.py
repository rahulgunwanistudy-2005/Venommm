"""The coverage model.

An antivenom is an estimator trained on the immunogen mixture. Neutralisation is a limiting-reagent
problem, not a geometric one: a fixed mass of antibody is partitioned across toxin families in
proportion to the immunising mixture, and each family of the incoming venom is neutralised only to
the extent that enough cross-reactive antibody was raised against it.

    supply_f(S, w)   = B * sum_q w_q a_qf                      antibody mass allocated to family f
    eff_f(S, w, p)   = B * sum_q w_q a_qf x_f(q -> p)           of which is usable against p
    demand_f(p, D)   = D * a_pf * kappa_f                       antibody mass family f of p requires
    n_f              = min(1, eff_f / demand_f)                 saturating neutralised fraction
    C(p | S, w, D)   = sum_f sigma_f a_pf n_f / sum_f sigma_f a_pf

Three properties follow from this form and none of them are fitted:

* **Non-monotone in |S|.** Adding a venom to the mixture moves antibody supply away from the
  families a given bite actually delivers. Cosine similarity between composition vectors cannot
  produce this, which is the ablation in `validate/baselines.py`.
* **Dose-dependent.** Demand is linear in D while supply is not, so coverage falls as dose rises.
* **Saturating.** Antibody raised in excess of demand for a family is wasted, encoded by min(1, .).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

from venomgap.config import (
    DEFAULT_B_MG,
    DEFAULT_D_MG,
    EPS,
    FAMILIES,
    LIMITING_FAMILY_MIN_WEIGHTED_SHARE,
    SIGMA,
)
from venomgap.errors import DataValidationError
from venomgap.types import CoverageResult, ResultFlag

N_FAMILIES = len(FAMILIES)

SIGMA_VECTOR: NDArray[np.float64] = np.array([SIGMA[f] for f in FAMILIES], dtype=np.float64)
SIGMA_VECTOR.setflags(write=False)  # frozen a priori; a writeable global would be a tuning vector


@dataclass(frozen=True, slots=True)
class CoverageInputs:
    """Dense, index-aligned inputs to the coverage calculation.

    Keeping this as arrays rather than dicts is what lets the optimiser evaluate hundreds of
    thousands of candidate mixtures. `immunogen_compositions` is (n_immunogen, n_families) and
    `crossreact` is (n_immunogen, n_families): entry (q, f) is x_f(q -> p) for the single target p.
    """

    target_composition: NDArray[np.float64]
    immunogen_compositions: NDArray[np.float64]
    crossreact: NDArray[np.float64]
    kappa: NDArray[np.float64]
    sigma: NDArray[np.float64] = field(default_factory=lambda: SIGMA_VECTOR)

    def __post_init__(self) -> None:
        if self.target_composition.shape != (N_FAMILIES,):
            raise DataValidationError(
                f"target composition has shape {self.target_composition.shape}, "
                f"expected ({N_FAMILIES},)"
            )
        shape = self.immunogen_compositions.shape
        if self.immunogen_compositions.ndim != 2 or shape[1] != N_FAMILIES:
            raise DataValidationError(
                f"immunogen compositions have shape {self.immunogen_compositions.shape}, "
                f"expected (n, {N_FAMILIES})"
            )
        if self.crossreact.shape != self.immunogen_compositions.shape:
            raise DataValidationError(
                f"crossreact shape {self.crossreact.shape} does not match immunogen "
                f"compositions {self.immunogen_compositions.shape}"
            )
        if self.kappa.shape != (N_FAMILIES,):
            raise DataValidationError(
                f"kappa has shape {self.kappa.shape}, expected ({N_FAMILIES},)"
            )
        if np.any(self.kappa <= 0.0):
            raise DataValidationError("kappa must be strictly positive for every family")
        if np.any(self.crossreact < 0.0) or np.any(self.crossreact > 1.0):
            raise DataValidationError("cross-reactivity must lie in [0, 1]")


def antibody_supply(
    immunogen_compositions: NDArray[np.float64],
    weights: NDArray[np.float64],
    B: float = DEFAULT_B_MG,
) -> NDArray[np.float64]:
    """supply_f(S, w) = B * sum_q w_q * a_qf.

    The finite total `B` is the whole mechanism: the returned vector always sums to B (up to the
    weights summing to 1), so antibody given to one family is antibody taken from another.
    """
    _check_weights(weights, immunogen_compositions.shape[0])
    return B * (weights @ immunogen_compositions)


def effective_supply(
    immunogen_compositions: NDArray[np.float64],
    weights: NDArray[np.float64],
    crossreact: NDArray[np.float64],
    B: float = DEFAULT_B_MG,
) -> NDArray[np.float64]:
    """eff_f(S, w, p) = B * sum_q w_q * a_qf * x_f(q -> p).

    Cross-reactivity discounts supply but does not return it to the budget: antibody raised against
    a divergent orthologue is spent, not recycled. That asymmetry is why the dilution penalty bites.
    """
    _check_weights(weights, immunogen_compositions.shape[0])
    return B * (weights @ (immunogen_compositions * crossreact))


def antibody_demand(
    target_composition: NDArray[np.float64],
    kappa: NDArray[np.float64],
    D: float = DEFAULT_D_MG,
) -> NDArray[np.float64]:
    """demand_f(p) = D * a_pf * kappa_f."""
    return D * target_composition * kappa


def neutralised_fractions(
    inputs: CoverageInputs,
    weights: NDArray[np.float64],
    B: float = DEFAULT_B_MG,
    D: float = DEFAULT_D_MG,
) -> NDArray[np.float64]:
    """n_f = min(1, eff_f / demand_f), with n_f = 1 where the venom contains none of family f.

    A family the snake does not deliver cannot be under-neutralised, so its demand is zero and its
    contribution to the weighted mean is zero as well; setting n_f = 1 there keeps the ratio
    well-defined without affecting the result.
    """
    eff = effective_supply(inputs.immunogen_compositions, weights, inputs.crossreact, B=B)
    demand = antibody_demand(inputs.target_composition, inputs.kappa, D=D)
    n = np.ones(N_FAMILIES, dtype=np.float64)
    present = demand > EPS
    n[present] = np.minimum(1.0, eff[present] / demand[present])
    return n


def coverage(
    inputs: CoverageInputs,
    weights: NDArray[np.float64],
    B: float = DEFAULT_B_MG,
    D: float = DEFAULT_D_MG,
) -> float:
    """C(p | S, w, D): the medically weighted fraction of the delivered venom dose neutralised."""
    n = neutralised_fractions(inputs, weights, B=B, D=D)
    weight = inputs.sigma * inputs.target_composition
    denom = float(weight.sum())
    if denom <= EPS:
        raise DataValidationError(
            "target composition has no medically weighted mass; coverage is undefined"
        )
    return float(np.clip((weight * n).sum() / denom, 0.0, 1.0))


def limiting_families(
    inputs: CoverageInputs,
    neutralised: NDArray[np.float64],
    min_weighted_share: float = LIMITING_FAMILY_MIN_WEIGHTED_SHARE,
) -> tuple[str, ...]:
    """Families with n_f < 1, worst first, restricted to those that matter medically.

    The filter stops a trace family with a rounding-level abundance from being reported as the
    reason a district is red. Without it the answer to "why" is noise.
    """
    weight = inputs.sigma * inputs.target_composition
    denom = float(weight.sum())
    if denom <= EPS:
        return ()
    share = weight / denom
    # Sort by neutralised fraction, then by how much of the deficit the family actually accounts
    # for. The tie-break matters: where cross-reactivity is zero many families sit at exactly
    # n_f = 0, and without it the "worst limiting family" would be decided alphabetically.
    candidates = [
        (float(neutralised[i]), -float(share[i] * (1.0 - neutralised[i])), FAMILIES[i])
        for i in range(N_FAMILIES)
        if neutralised[i] < 1.0 - 1e-9 and share[i] >= min_weighted_share
    ]
    candidates.sort()
    return tuple(name for _, _, name in candidates)


def evaluate(
    pop_id: str,
    inputs: CoverageInputs,
    weights: NDArray[np.float64],
    B: float = DEFAULT_B_MG,
    D: float = DEFAULT_D_MG,
    uncertainty: float = 0.0,
    flags: tuple[ResultFlag, ...] = (),
) -> CoverageResult:
    """Full per-population result, including why coverage is short where it is."""
    n = neutralised_fractions(inputs, weights, B=B, D=D)
    c = coverage(inputs, weights, B=B, D=D)
    return CoverageResult(
        pop_id=pop_id,
        coverage=c,
        deficit=1.0 - c,
        per_family_neutralised={FAMILIES[i]: float(n[i]) for i in range(N_FAMILIES)},
        limiting_families=limiting_families(inputs, n),
        uncertainty=uncertainty,
        flags=flags,
    )


def coverage_many(
    target_compositions: NDArray[np.float64],
    immunogen_compositions: NDArray[np.float64],
    crossreact: NDArray[np.float64],
    weights: NDArray[np.float64],
    kappa: NDArray[np.float64],
    sigma: NDArray[np.float64] | None = None,
    B: float = DEFAULT_B_MG,
    D: float = DEFAULT_D_MG,
) -> NDArray[np.float64]:
    """Vectorised coverage for many targets at once.

    `crossreact` is (n_targets, n_immunogen, n_families). This is the hot path: the optimiser calls
    it once per candidate mixture over every district, so the loop stays in numpy.
    """
    sigma = SIGMA_VECTOR if sigma is None else sigma
    n_targets = target_compositions.shape[0]
    n_immunogen = immunogen_compositions.shape[0]
    if target_compositions.shape[1] != N_FAMILIES:
        raise DataValidationError("target compositions must be (n_targets, n_families)")
    if crossreact.shape != (n_targets, n_immunogen, N_FAMILIES):
        raise DataValidationError(
            f"crossreact shape {crossreact.shape} does not match "
            f"({n_targets}, {n_immunogen}, {N_FAMILIES})"
        )
    _check_weights(weights, n_immunogen)

    # eff[t, f] = B * sum_q w_q * a_qf * x[t, q, f]
    weighted = immunogen_compositions[None, :, :] * crossreact
    eff = B * np.einsum("q,tqf->tf", weights, weighted)
    demand = D * target_compositions * kappa[None, :]

    n = np.ones_like(demand)
    present = demand > EPS
    np.divide(eff, demand, out=n, where=present)
    np.minimum(n, 1.0, out=n)
    n[~present] = 1.0

    weight = sigma[None, :] * target_compositions
    denom = weight.sum(axis=1)
    if np.any(denom <= EPS):
        raise DataValidationError("a target composition has no medically weighted mass")
    result: NDArray[np.float64] = np.clip((weight * n).sum(axis=1) / denom, 0.0, 1.0)
    return result


def _check_weights(weights: NDArray[np.float64], expected: int) -> None:
    if weights.shape != (expected,):
        raise DataValidationError(
            f"weights have shape {weights.shape}, expected ({expected},)"
        )
    if np.any(weights < -EPS):
        raise DataValidationError("mixture weights must be non-negative")
    total = float(weights.sum())
    if abs(total - 1.0) > 1e-6:
        raise DataValidationError(f"mixture weights sum to {total:.9f}, expected 1")
