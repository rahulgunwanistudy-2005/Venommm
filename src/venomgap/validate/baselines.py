"""The cosine-similarity baseline, present only to be shown failing.

The existing literature scores antivenom-venom match by geometric similarity between composition
vectors. That is what this module implements, faithfully, so that the comparison is fair:

    C_cosine(p | S, w) = sum_q w_q * cos(a_q, a_p)

and its natural paraspecific variant, the best single match in the mixture:

    C_cosine_max(p | S) = max_q cos(a_q, a_p)

Both are **monotone non-decreasing in the mixture size**, and the reason is structural rather than
empirical. The weighted-mean form is an average of fixed per-venom scores, so adding a venom can
only
move it toward that venom's score -- and if you are free to choose the best subset at each size, a
larger subset can always retain the previous best and so never scores lower. The max form is
monotone by inspection. Neither has a budget, so neither can represent antibody raised against one
family being antibody not raised against another.

That is criterion R4c. A metric with no finite supply cannot produce a dilution penalty, so it
cannot
reproduce an experimentally observed phenomenon that the stoichiometric model gets for free.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from itertools import combinations
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from venomgap.config import EPS, FAMILIES, MAX_K, SIGMA
from venomgap.model.assemble import ModelAssembly
from venomgap.optimize.objective import NationalObjective
from venomgap.types import FittedParameters

if TYPE_CHECKING:
    from venomgap.validate.retrodiction import MixtureSizeCurve

logger = logging.getLogger(__name__)

SIGMA_VECTOR: NDArray[np.float64] = np.array([SIGMA[f] for f in FAMILIES], dtype=np.float64)


def cosine_similarity(a: NDArray[np.float64], b: NDArray[np.float64]) -> float:
    """Plain cosine between two composition vectors."""
    na, nb = float(np.linalg.norm(a)), float(np.linalg.norm(b))
    if na <= EPS or nb <= EPS:
        return 0.0
    return float(np.clip(np.dot(a, b) / (na * nb), 0.0, 1.0))


def cosine_coverage(
    target: NDArray[np.float64],
    immunogens: NDArray[np.float64],
    weights: NDArray[np.float64],
    use_max: bool = False,
) -> float:
    """The baseline score. No budget, no dose, no stoichiometry -- by design."""
    scores = np.array([cosine_similarity(immunogens[q], target) for q in range(len(immunogens))])
    if use_max:
        return float(scores.max())
    return float(np.dot(weights, scores))


@dataclass(frozen=True, slots=True)
class CosineObjective:
    """National coverage under the cosine baseline, on exactly the same district grid."""

    objective: NationalObjective
    assembly: ModelAssembly
    use_max: bool = False

    def national_coverage(self, pop_ids: tuple[str, ...]) -> float:
        if not pop_ids:
            return 0.0
        immunogen = self.assembly.immunogen_from_pop_ids(pop_ids)
        cells = np.array(
            [
                cosine_coverage(
                    self.objective.cell_compositions[i],
                    immunogen.compositions,
                    immunogen.weights,
                    use_max=self.use_max,
                )
                for i in range(self.objective.n_cells)
            ]
        )
        numerator = np.bincount(
            self.objective.cell_district_index,
            weights=cells * self.objective.cell_weights,
            minlength=len(self.objective.district_ids),
        )
        denominator = np.bincount(
            self.objective.cell_district_index,
            weights=self.objective.cell_weights,
            minlength=len(self.objective.district_ids),
        )
        district = np.zeros_like(numerator)
        present = denominator > 0.0
        district[present] = numerator[present] / denominator[present]
        return float(np.dot(self.objective.burden, district))


def cosine_mixture_size_curve(
    objective: NationalObjective,
    assembly: ModelAssembly,
    fitted: FittedParameters,  # noqa: ARG001 - the baseline has no fitted parameters, by design
    max_size: int = MAX_K,
    use_max: bool = False,
) -> MixtureSizeCurve:
    """The R4c curve: the cosine baseline's best achievable score at each mixture size."""
    from venomgap.optimize.greedy import _immunogen_candidate_pool
    from venomgap.validate.retrodiction import MixtureSizeCurve

    pool = list(_immunogen_candidate_pool(assembly))
    baseline = CosineObjective(objective=objective, assembly=assembly, use_max=use_max)

    sizes: list[int] = []
    coverages: list[float] = []
    best_sets: list[tuple[str, ...]] = []
    for size in range(1, min(max_size, len(pool)) + 1):
        best: tuple[float, tuple[str, ...]] | None = None
        for subset in combinations(pool, size):
            value = baseline.national_coverage(subset)
            if best is None or value > best[0]:
                best = (value, subset)
        assert best is not None
        sizes.append(size)
        coverages.append(best[0])
        best_sets.append(best[1])
    logger.info(
        "cosine baseline curve: %s",
        ", ".join(f"{s}:{c:.4f}" for s, c in zip(sizes, coverages, strict=True)),
    )
    return MixtureSizeCurve(
        sizes=tuple(sizes), coverages=tuple(coverages), best_sets=tuple(best_sets)
    )


def cosine_deficit_ranking(
    objective: NationalObjective,  # noqa: ARG001 - kept for signature parity with the model path
    assembly: ModelAssembly,
    use_max: bool = False,
) -> dict[str, float]:
    """Per-population deficit under the baseline, for the R1-R3 ablation comparison."""
    immunogen = assembly.big_four_immunogen()
    out: dict[str, float] = {}
    for population in assembly.populations:
        if population.country != "India":
            continue
        score = cosine_coverage(
            population.vector(), immunogen.compositions, immunogen.weights, use_max=use_max
        )
        out[population.pop_id] = 1.0 - score
    return out


def sequence_free_deficit_ranking(
    assembly: ModelAssembly, fitted: FittedParameters
) -> dict[str, float]:
    """Pre-registered negative control: the stoichiometric model with x_f = 1 for every pair.

    If cross-reactivity carried no signal, this would rank the holdout targets as well as the full
    model does. It is the check that the sequence layer is doing work rather than decorating.
    """
    from venomgap.config import DEFAULT_B_MG, DEFAULT_D_MG
    from venomgap.model.coverage import CoverageInputs, coverage

    immunogen = assembly.big_four_immunogen()
    kappa = assembly.kappa(fitted.kappa_scale)
    ones = np.ones((immunogen.size, len(FAMILIES)))
    out: dict[str, float] = {}
    for population in assembly.populations:
        if population.country != "India":
            continue
        inputs = CoverageInputs(
            target_composition=population.vector(),
            immunogen_compositions=immunogen.compositions,
            crossreact=ones,
            kappa=kappa,
            sigma=SIGMA_VECTOR,
        )
        value = coverage(inputs, immunogen.weights, B=DEFAULT_B_MG, D=DEFAULT_D_MG)
        out[population.pop_id] = 1.0 - value
    return out
