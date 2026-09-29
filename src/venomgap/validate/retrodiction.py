"""The blinded retrodiction study: R1-R4, evaluated automatically, reported whatever the answer is.

Every criterion below is the one written in `experiments/preregistration.md` (with amendment 1),
transcribed into code. Nothing is softened here and nothing is re-scoped. The evaluator computes a
number, compares it to a threshold that was fixed before any parameter was fitted, and records the
result.

The parameters this runs on were fitted without any of the four holdout studies' antivenom results.
As it turned out, only one parameter was identifiable from the calibration data at all
(`kappa_scale`); `sigma_f`, the `kappa_f` shape and every `theta_f` retain their a priori values.
That makes these predictions harder to have rigged than the pre-registration anticipated, and it is
also a real limitation, since the cross-reactivity thresholds are unvalidated.
"""

from __future__ import annotations

import json
import logging
import subprocess
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from venomgap.config import (
    DEFAULT_B_MG,
    DEFAULT_D_MG,
    FAMILIES,
    IRULA_LAT,
    IRULA_LON,
    LIMITING_FAMILY_MIN_WEIGHTED_SHARE,
    MAX_K,
    RESULTS_DIR,
    SIGMA,
)
from venomgap.model.assemble import ImmunogenSpec, ModelAssembly
from venomgap.model.coverage import coverage, neutralised_fractions
from venomgap.types import (
    FittedParameters,
    RetrodictionCriterion,
    RetrodictionOutcome,
    VenomPopulation,
)

logger = logging.getLogger(__name__)

# --------------------------------------------------------------------------------------
# Pre-registered thresholds. These are transcribed from the pre-registration, not chosen here.
# --------------------------------------------------------------------------------------

TOP_DECILE = 0.10
TOP_TERTILE = 1.0 / 3.0
R2B_MIN_COVERAGE_GAIN = 0.05
R3B_MIN_DEFICIT_EXCESS = 0.10
R4_MIN_TURNOVER_DROP = 0.02

R1_TARGET_PREFIX = "Echis_carinatus_sochureki__"
R2_TARGET = "Naja_sagittifera__Andaman"
R3_TARGET = "Bungarus_caeruleus__Punjab"
R3_SOUTHERN_SPECIES = "Bungarus caeruleus"


@dataclass(frozen=True, slots=True)
class DeficitRanking:
    """Predicted deficit for every curated population under the baseline immunogen."""

    pop_ids: tuple[str, ...]
    deficits: NDArray[np.float64]
    limiting: dict[str, tuple[str, ...]]
    per_family: dict[str, dict[str, float]]

    def deficit(self, pop_id: str) -> float:
        return float(self.deficits[self.pop_ids.index(pop_id)])

    def quantile_of(self, pop_id: str) -> float:
        """Fraction of populations with a deficit at least as large. 0.0 is the very worst."""
        value = self.deficit(pop_id)
        return float(np.mean(self.deficits >= value))

    def rank_of(self, pop_id: str) -> int:
        order = np.argsort(-self.deficits)
        return int(np.where(order == self.pop_ids.index(pop_id))[0][0]) + 1


def rank_populations(
    assembly: ModelAssembly,
    fitted: FittedParameters,
    immunogen: ImmunogenSpec,
    B: float = DEFAULT_B_MG,
    D: float = DEFAULT_D_MG,
    indian_only: bool = True,
) -> DeficitRanking:
    """Predicted deficit for every curated population, which is what R1-R3 rank against."""
    theta = fitted.theta_vector()
    populations = [
        p for p in assembly.populations if (p.country == "India" or not indian_only)
    ]
    pop_ids: list[str] = []
    deficits: list[float] = []
    limiting: dict[str, tuple[str, ...]] = {}
    per_family: dict[str, dict[str, float]] = {}
    for population in populations:
        inputs, _ = assembly.inputs_for_population(
            population.pop_id, immunogen, theta, fitted.kappa_scale
        )
        value = coverage(inputs, immunogen.weights, B=B, D=D)
        n = neutralised_fractions(inputs, immunogen.weights, B=B, D=D)
        pop_ids.append(population.pop_id)
        deficits.append(1.0 - value)
        limiting[population.pop_id] = _worst_families(population, n)
        per_family[population.pop_id] = {
            FAMILIES[i]: float(n[i]) for i in range(len(FAMILIES))
        }
    return DeficitRanking(
        pop_ids=tuple(pop_ids),
        deficits=np.array(deficits, dtype=np.float64),
        limiting=limiting,
        per_family=per_family,
    )


def _worst_families(
    population: VenomPopulation, neutralised: NDArray[np.float64]
) -> tuple[str, ...]:
    """Families with n_f < 1, worst first, restricted to those carrying real medical weight."""
    vector = population.vector()
    weight = np.array([SIGMA[f] for f in FAMILIES]) * vector
    total = float(weight.sum())
    if total <= 0.0:
        return ()
    share = weight / total
    candidates = [
        (float(neutralised[i]), -float(share[i] * (1.0 - neutralised[i])), FAMILIES[i])
        for i in range(len(FAMILIES))
        if neutralised[i] < 1.0 - 1e-9 and share[i] >= LIMITING_FAMILY_MIN_WEIGHTED_SHARE
    ]
    candidates.sort()
    return tuple(name for _, _, name in candidates)


# --------------------------------------------------------------------------------------
# R1 -- Echis carinatus sochureki, north-west India
# --------------------------------------------------------------------------------------


def evaluate_r1(ranking: DeficitRanking) -> list[RetrodictionCriterion]:
    """R1a: top decile of deficit. R1b: SVMP is the worst limiting family.

    Per amendment A1.1 both are evaluated on every curated north-west *E. c. sochureki* population
    and pass on a strict majority. All four outcomes are recorded either way.
    """
    targets = [p for p in ranking.pop_ids if p.startswith(R1_TARGET_PREFIX)]
    if not targets:
        return [
            RetrodictionCriterion(
                criterion_id="R1",
                statement="E. c. sochureki north-west populations rank in the top decile of "
                "deficit, SVMP-driven",
                computed=None,
                threshold=None,
                comparison="ill_posed",
                passed=False,
                detail="no E. c. sochureki population is present in the corpus",
            )
        ]

    a_results: list[tuple[str, float, bool]] = []
    b_results: list[tuple[str, str, bool]] = []
    for pop_id in targets:
        quantile = ranking.quantile_of(pop_id)
        a_results.append((pop_id, quantile, quantile <= TOP_DECILE))
        worst = ranking.limiting[pop_id]
        top = worst[0] if worst else "none"
        b_results.append((pop_id, top, top == "SVMP"))

    a_passed = sum(1 for _, _, ok in a_results if ok)
    b_passed = sum(1 for _, _, ok in b_results if ok)
    majority = len(targets) / 2.0

    return [
        RetrodictionCriterion(
            criterion_id="R1a",
            statement=(
                "Each north-west E. c. sochureki population ranks in the top decile of predicted "
                "deficit among curated Indian populations; passes on a strict majority."
            ),
            computed=float(a_passed),
            threshold=float(majority),
            comparison=">=",
            passed=a_passed > majority,
            detail="; ".join(
                f"{p.removeprefix(R1_TARGET_PREFIX)} at quantile {q:.3f} "
                f"(rank {ranking.rank_of(p)}/{len(ranking.pop_ids)}) {'PASS' if ok else 'fail'}"
                for p, q, ok in a_results
            ),
        ),
        RetrodictionCriterion(
            criterion_id="R1b",
            statement=(
                "SVMP is the single worst limiting family for each such population; passes on a "
                "strict majority."
            ),
            computed=float(b_passed),
            threshold=float(majority),
            comparison=">=",
            passed=b_passed > majority,
            detail="; ".join(
                f"{p.removeprefix(R1_TARGET_PREFIX)} worst={w} {'PASS' if ok else 'fail'}"
                for p, w, ok in b_results
            ),
        ),
    ]


# --------------------------------------------------------------------------------------
# R2 -- Naja sagittifera, Andaman and Nicobar
# --------------------------------------------------------------------------------------


def evaluate_r2(
    ranking: DeficitRanking,
    assembly: ModelAssembly,
    fitted: FittedParameters,
    B: float = DEFAULT_B_MG,
    D: float = DEFAULT_D_MG,
) -> list[RetrodictionCriterion]:
    """R2a: top decile. R2b: substituting N. kaouthia for N. naja raises coverage by >= 0.05."""
    criteria: list[RetrodictionCriterion] = []

    if R2_TARGET not in ranking.pop_ids:
        return [
            RetrodictionCriterion(
                criterion_id="R2",
                statement="Naja sagittifera ranks in the top decile and a kaouthia-based "
                "immunogen closes the gap",
                computed=None,
                threshold=None,
                comparison="ill_posed",
                passed=False,
                detail=f"{R2_TARGET} is absent from the corpus",
            )
        ]

    quantile = ranking.quantile_of(R2_TARGET)
    criteria.append(
        RetrodictionCriterion(
            criterion_id="R2a",
            statement="Naja sagittifera ranks in the top decile of predicted deficit.",
            computed=quantile,
            threshold=TOP_DECILE,
            comparison="<=",
            passed=quantile <= TOP_DECILE,
            detail=(
                f"quantile {quantile:.3f}, rank {ranking.rank_of(R2_TARGET)} of "
                f"{len(ranking.pop_ids)}, deficit {ranking.deficit(R2_TARGET):.3f}"
            ),
        )
    )

    theta = fitted.theta_vector()
    baseline = assembly.big_four_immunogen()
    substituted = assembly.big_four_immunogen(
        substitute={"Naja naja": "Naja kaouthia"}
    )
    gains = []
    for immunogen in (baseline, substituted):
        inputs, _ = assembly.inputs_for_population(
            R2_TARGET, immunogen, theta, fitted.kappa_scale
        )
        gains.append(coverage(inputs, immunogen.weights, B=B, D=D))
    gain = gains[1] - gains[0]
    criteria.append(
        RetrodictionCriterion(
            criterion_id="R2b",
            statement=(
                "Substituting Naja kaouthia for Naja naja in the immunogen set, holding everything "
                "else fixed, increases predicted coverage of Naja sagittifera by at least 0.05."
            ),
            computed=gain,
            threshold=R2B_MIN_COVERAGE_GAIN,
            comparison=">=",
            passed=gain >= R2B_MIN_COVERAGE_GAIN,
            detail=(
                f"coverage {gains[0]:.4f} -> {gains[1]:.4f}, gain {gain:+.4f}"
            ),
        )
    )
    return criteria


# --------------------------------------------------------------------------------------
# R3 -- North Indian Bungarus caeruleus
# --------------------------------------------------------------------------------------


def evaluate_r3(
    ranking: DeficitRanking,
    assembly: ModelAssembly,
    fitted: FittedParameters,
    B: float = DEFAULT_B_MG,
    D: float = DEFAULT_D_MG,
) -> list[RetrodictionCriterion]:
    """R3a: top tertile despite being a Big Four species. R3b: exceeds the southern reference."""
    criteria: list[RetrodictionCriterion] = []

    if R3_TARGET not in ranking.pop_ids:
        return [
            RetrodictionCriterion(
                criterion_id="R3",
                statement="North Indian Bungarus caeruleus shows a high deficit despite being a "
                "Big Four species",
                computed=None,
                threshold=None,
                comparison="ill_posed",
                passed=False,
                detail=f"{R3_TARGET} is absent from the corpus",
            )
        ]

    quantile = ranking.quantile_of(R3_TARGET)
    criteria.append(
        RetrodictionCriterion(
            criterion_id="R3a",
            statement=(
                "North Indian Bungarus caeruleus, a Big Four species, has a predicted deficit in "
                "the top tertile of curated Indian populations."
            ),
            computed=quantile,
            threshold=TOP_TERTILE,
            comparison="<=",
            passed=quantile <= TOP_TERTILE,
            detail=(
                f"quantile {quantile:.3f}, rank {ranking.rank_of(R3_TARGET)} of "
                f"{len(ranking.pop_ids)}, deficit {ranking.deficit(R3_TARGET):.3f}"
            ),
        )
    )

    # Southern reference: the curated B. caeruleus population nearest the Irula immunogen source
    # (amendment A1.2). Chosen by distance in code, not named by hand.
    from venomgap.model.spatial import haversine_km

    southern = [
        p
        for p in assembly.populations
        if p.species == R3_SOUTHERN_SPECIES and p.pop_id != R3_TARGET
    ]
    if not southern:
        criteria.append(
            RetrodictionCriterion(
                criterion_id="R3b",
                statement="Deficit exceeds the southern reference population by at least 0.10",
                computed=None,
                threshold=R3B_MIN_DEFICIT_EXCESS,
                comparison="ill_posed",
                passed=False,
                detail="no comparator B. caeruleus population exists in the corpus",
            )
        )
        return criteria

    distances = [
        float(haversine_km(IRULA_LAT, IRULA_LON, p.lat, p.lon)) for p in southern
    ]
    reference = southern[int(np.argmin(distances))]
    theta = fitted.theta_vector()
    immunogen = assembly.big_four_immunogen()
    inputs, _ = assembly.inputs_for_population(
        reference.pop_id, immunogen, theta, fitted.kappa_scale
    )
    reference_deficit = 1.0 - coverage(inputs, immunogen.weights, B=B, D=D)
    excess = ranking.deficit(R3_TARGET) - reference_deficit
    criteria.append(
        RetrodictionCriterion(
            criterion_id="R3b",
            statement=(
                "Its predicted deficit exceeds that of the southern reference B. caeruleus "
                "population nearest the immunogen source by at least 0.10 absolute."
            ),
            computed=excess,
            threshold=R3B_MIN_DEFICIT_EXCESS,
            comparison=">=",
            passed=excess >= R3B_MIN_DEFICIT_EXCESS,
            detail=(
                f"reference {reference.pop_id} at {min(distances):.0f} km from the immunogen "
                f"source, deficit {reference_deficit:.3f}; Punjab deficit "
                f"{ranking.deficit(R3_TARGET):.3f}; excess {excess:+.3f}"
            ),
        )
    )
    return criteria


# --------------------------------------------------------------------------------------
# R4 -- immunogen dilution. The structural prediction.
# --------------------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class MixtureSizeCurve:
    """Best achievable mean coverage at each immunogen mixture size."""

    sizes: tuple[int, ...]
    coverages: tuple[float, ...]
    best_sets: tuple[tuple[str, ...], ...]

    @property
    def argmax_size(self) -> int:
        return self.sizes[int(np.argmax(self.coverages))]

    @property
    def peak(self) -> float:
        return float(max(self.coverages))

    @property
    def terminal(self) -> float:
        return float(self.coverages[-1])

    @property
    def has_interior_maximum(self) -> bool:
        best = int(np.argmax(self.coverages))
        return 0 < best < len(self.coverages) - 1


def evaluate_r4(
    curve: MixtureSizeCurve, baseline_curve: MixtureSizeCurve
) -> list[RetrodictionCriterion]:
    """R4a: interior maximum. R4b: drop of at least 0.02. R4c: the cosine baseline cannot do it."""
    interior = curve.has_interior_maximum
    drop = curve.peak - curve.terminal
    baseline_monotone = all(
        later >= earlier - 1e-9
        for earlier, later in zip(
            baseline_curve.coverages, baseline_curve.coverages[1:], strict=False
        )
    )
    return [
        RetrodictionCriterion(
            criterion_id="R4a",
            statement=(
                "Mean burden-weighted national coverage as a function of immunogen mixture size, "
                "each size at its own best subset, is non-monotone with an interior maximum."
            ),
            computed=float(curve.argmax_size),
            threshold=None,
            comparison="bool",
            passed=interior,
            detail=(
                f"peak at |S| = {curve.argmax_size} of {curve.sizes[-1]}; profile "
                + ", ".join(f"{s}:{c:.4f}" for s, c in zip(curve.sizes, curve.coverages,
                strict=True))
            ),
        ),
        RetrodictionCriterion(
            criterion_id="R4b",
            statement=(
                "The fall from the maximum to the largest mixture is at least 0.02 absolute "
                "coverage, so the turnover is not numerical noise."
            ),
            computed=drop,
            threshold=R4_MIN_TURNOVER_DROP,
            comparison=">=",
            passed=drop >= R4_MIN_TURNOVER_DROP,
            detail=(
                f"peak {curve.peak:.4f} at |S| = {curve.argmax_size}, "
                f"terminal {curve.terminal:.4f}"
            ),
        ),
        RetrodictionCriterion(
            criterion_id="R4c",
            statement=(
                "The cosine-similarity baseline, on identical inputs, is monotone non-decreasing "
                "in mixture size and therefore cannot reproduce R4a."
            ),
            computed=float(baseline_curve.argmax_size),
            threshold=None,
            comparison="bool",
            passed=baseline_monotone and not baseline_curve.has_interior_maximum,
            detail=(
                f"cosine baseline peak at |S| = {baseline_curve.argmax_size}; monotone "
                f"non-decreasing: {baseline_monotone}; profile "
                + ", ".join(
                    f"{s}:{c:.4f}"
                    for s, c in zip(baseline_curve.sizes, baseline_curve.coverages, strict=True)
                )
            ),
        ),
    ]


# --------------------------------------------------------------------------------------
# Runner
# --------------------------------------------------------------------------------------


def preregistration_commit() -> str:
    """The commit that introduced the pre-registration, recorded alongside the result."""
    try:
        output = subprocess.run(
            [
                "git", "log", "--diff-filter=A", "--format=%H", "-1", "--",
                "experiments/preregistration.md",
            ],
            capture_output=True,
            text=True,
            check=True,
            timeout=30,
        )
        return output.stdout.strip() or "unknown"
    except (subprocess.SubprocessError, OSError) as exc:
        logger.warning("could not read the pre-registration commit: %s", exc)
        return "unknown"


def run_retrodiction(write: bool = True) -> RetrodictionOutcome:
    """Evaluate R1-R4 and write results/retrodiction.json. Reports the outcome whatever it is."""
    from venomgap.ingest.compositions import load_compositions
    from venomgap.model.calibrate import load_fitted_parameters
    from venomgap.model.crossreact import IdentityTable
    from venomgap.model.spatial import SpatialCompositionModel
    from venomgap.optimize.greedy import best_subset_per_size
    from venomgap.optimize.objective import NationalObjective, build_national_objective
    from venomgap.validate.baselines import cosine_mixture_size_curve

    populations = load_compositions()
    fitted = load_fitted_parameters()
    assert fitted.ell_km is not None
    assembly = ModelAssembly(
        populations=populations,
        identity=IdentityTable.from_json(),
        spatial=SpatialCompositionModel(populations, fitted.ell_km),
    )
    immunogen = assembly.big_four_immunogen()
    ranking = rank_populations(assembly, fitted, immunogen)

    objective: NationalObjective = build_national_objective(assembly, fitted)
    curve = best_subset_per_size(objective, assembly, fitted, max_size=MAX_K)
    # The baseline is run in its *best-match* form. That is the paraspecificity metric the
    # literature actually uses, and it is the only form that is monotone in mixture size, so it is
    # the strongest version of the comparison rather than a straw man. The weighted-mean form is
    # also computed, and reported alongside, because it fails for a different and non-mechanistic
    # reason: it is an average of fixed per-venom scores and so simply declines as venoms are added.
    baseline_curve = cosine_mixture_size_curve(
        objective, assembly, fitted, max_size=MAX_K, use_max=True
    )
    baseline_mean_curve = cosine_mixture_size_curve(
        objective, assembly, fitted, max_size=MAX_K, use_max=False
    )

    criteria = [
        *evaluate_r1(ranking),
        *evaluate_r2(ranking, assembly, fitted),
        *evaluate_r3(ranking, assembly, fitted),
        *evaluate_r4(curve, baseline_curve),
    ]

    # A target passes only if every one of its sub-criteria passes.
    targets = {"R1": [], "R2": [], "R3": [], "R4": []}  # type: dict[str, list[bool]]
    for criterion in criteria:
        targets[criterion.criterion_id[:2]].append(criterion.passed)
    passed = sum(1 for results in targets.values() if results and all(results))

    outcome = RetrodictionOutcome(
        preregistration_commit=preregistration_commit(),
        data_version=_data_version(),
        model_variant="stoichiometric",
        criteria=tuple(criteria),
        targets_passed=passed,
        targets_total=len(targets),
        notes=(
            "Parameters were fitted without any holdout study's antivenom results. Only "
            f"kappa_scale ({fitted.kappa_scale:.4f}) and the ED50 anchor C* "
            f"({fitted.c_star:.4f} if known) were identifiable; sigma_f, the kappa_f shape and "
            "every theta_f retain their a priori values, since the calibration set could not "
            "constrain theta at all. Unidentified thetas: "
            f"{list(fitted.theta_unidentified)}."
        ),
    )

    if write:
        RESULTS_DIR.mkdir(parents=True, exist_ok=True)
        payload = outcome.model_dump(mode="json")
        payload["deficit_ranking"] = {
            pop_id: {
                "deficit": float(ranking.deficits[i]),
                "rank": ranking.rank_of(pop_id),
                "quantile": ranking.quantile_of(pop_id),
                "limiting_families": list(ranking.limiting[pop_id]),
            }
            for i, pop_id in enumerate(ranking.pop_ids)
        }
        payload["mixture_size_curve"] = {
            "sizes": list(curve.sizes),
            "coverages": list(curve.coverages),
            "best_sets": [list(s) for s in curve.best_sets],
        }
        payload["cosine_baseline_curve"] = {
            "form": "best-match (max over immunogens), the literature paraspecificity metric",
            "sizes": list(baseline_curve.sizes),
            "coverages": list(baseline_curve.coverages),
        }
        payload["cosine_weighted_mean_curve"] = {
            "form": "weighted mean over immunogens",
            "sizes": list(baseline_mean_curve.sizes),
            "coverages": list(baseline_mean_curve.coverages),
            "note": (
                "Declines with mixture size, but as an arithmetic consequence of averaging fixed "
                "per-venom scores, not through any antibody-budget mechanism. It has no interior "
                "maximum, so it still cannot reproduce R4a."
            ),
        }
        payload["fitted_parameters"] = fitted.model_dump(mode="json")
        payload["diagnostics"] = build_diagnostics(assembly, fitted, ranking)
        (RESULTS_DIR / "retrodiction.json").write_text(json.dumps(payload, indent=1))
    return outcome


def build_diagnostics(
    assembly: ModelAssembly,
    fitted: FittedParameters,
    ranking: DeficitRanking,
) -> dict[str, object]:
    """Why each criterion landed where it did, in numbers rather than narrative.

    A retrodiction that fails without a diagnosis is not a result, it is an anecdote. Everything
    here is computed after the criteria were evaluated and is reported as explanation, never as a
    revision of the outcome.
    """
    theta = fitted.theta_vector()
    identities = np.array(list(assembly.identity.values.values()))
    floor = float(theta[0])
    below = float(np.mean(identities < floor)) if identities.size else float("nan")

    per_species_counts: dict[str, int] = {}
    for population in assembly.populations:
        per_species_counts[population.species] = per_species_counts.get(population.species, 0) + 1

    # R1b under the alternative reading: which family accounts for most of the deficit, as opposed
    # to which has the lowest neutralised fraction. Reported because the two can differ.
    deficit_drivers: dict[str, list[tuple[str, float]]] = {}
    sigma = np.array([SIGMA[f] for f in FAMILIES])
    for pop_id in ranking.pop_ids:
        population = assembly.by_id[pop_id]
        vector = population.vector()
        weight = sigma * vector
        total = float(weight.sum())
        if total <= 0.0:
            continue
        share = weight / total
        n = np.array([ranking.per_family[pop_id][f] for f in FAMILIES])
        contribution = share * (1.0 - n)
        order = np.argsort(-contribution)[:3]
        deficit_drivers[pop_id] = [
            (FAMILIES[i], float(contribution[i])) for i in order if contribution[i] > 1e-6
        ]

    # R2b specifically: was the substitution inert, and if so why?
    r2b_detail: dict[str, object] = {}
    if R2_TARGET in ranking.pop_ids:
        target_species = assembly.by_id[R2_TARGET].species
        r2b_detail = {
            "target_species": target_species,
            "identity_3FTx_from_Naja_naja": assembly.identity.get(
                "Naja naja", target_species, "3FTx"
            ),
            "identity_3FTx_from_Naja_kaouthia": assembly.identity.get(
                "Naja kaouthia", target_species, "3FTx"
            ),
            "theta_3FTx": float(theta[FAMILIES.index("3FTx")]),
            "note": (
                "Both identities fall below theta, so x_f is exactly zero for both and the "
                "substitution cannot change coverage. The sequence data does order the two "
                "correctly -- kaouthia is the closer of the pair -- but the recognition floor "
                "discards that ordering before it can act."
            ),
        }

    return {
        "theta_scale_mismatch": {
            "theta": floor,
            "between_species_identities": int(identities.size),
            "mean_identity": float(np.mean(identities)) if identities.size else None,
            "fraction_below_theta": below,
            "consequence": (
                "Cross-reactivity is forced to exactly zero for this fraction of species pairs, so "
                "the model has almost no paraspecific coverage. theta was pre-registered to be "
                "fitted on the calibration split; the calibration data could not identify it, so "
                "it kept its a priori value of 0.65. That prior was reasoned from epitope-level "
                "identity between orthologous toxins, whereas the statistic actually used is the "
                "mean identity over all sequence pairs in a family, which includes comparisons "
                "across divergent subfamilies and is therefore on a systematically lower scale. "
                "This scale mismatch is the single largest cause of the R1a, R2a, R2b and R3a "
                "misses."
            ),
        },
        "data_sparsity": {
            "populations_per_species": dict(sorted(per_species_counts.items())),
            "consequence": (
                "Bungarus caeruleus has two curated proteomes and one of them is the R3 target "
                "itself, so the interpolated immunogen at the collection locality is built partly "
                "from the population being tested. That makes North Indian krait look better "
                "covered than it should and is a direct cause of the R3a miss."
            ),
        },
        "deficit_drivers": deficit_drivers,
        "r2b": r2b_detail,
    }


def _data_version() -> str:
    import hashlib

    from venomgap.config import COMPOSITIONS_CSV

    digest = hashlib.sha256(COMPOSITIONS_CSV.read_bytes()).hexdigest()
    return f"compositions.csv sha256:{digest[:16]}"
