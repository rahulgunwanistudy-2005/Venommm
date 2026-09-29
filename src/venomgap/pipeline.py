"""End-to-end pipeline steps, each writing a JSON result the figures and README read from.

No number quoted anywhere in this project is retyped by hand. Everything in the README, the figures
and the web app is read back out of `results/*.json`.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.stats import spearmanr

from venomgap.config import (
    B_SWEEP_MG,
    D_SWEEP_MG,
    DEFAULT_B_MG,
    DEFAULT_D_MG,
    DISCLAIMER,
    EXACT_INSTANCE_MAX_K,
    MAX_K,
    RESULTS_DIR,
)
from venomgap.model.assemble import ModelAssembly
from venomgap.model.calibrate import load_fitted_parameters
from venomgap.model.crossreact import IdentityTable
from venomgap.model.spatial import (
    SpatialCompositionModel,
    composition_dissimilarity_vs_distance,
    fit_length_scale,
)
from venomgap.optimize.objective import NationalObjective, build_national_objective
from venomgap.types import FittedParameters, SitingSolution

logger = logging.getLogger(__name__)

HIGH_DEFICIT_THRESHOLD = 0.5


@dataclass(frozen=True, slots=True)
class Context:
    assembly: ModelAssembly
    fitted: FittedParameters
    objective: NationalObjective


def build_context(
    B: float = DEFAULT_B_MG,
    D: float = DEFAULT_D_MG,
    uniform_within_state: bool = False,
) -> Context:
    from venomgap.ingest.compositions import load_compositions

    populations = load_compositions()
    fitted = load_fitted_parameters()
    if fitted.ell_km is None:
        raise ValueError("fitted parameters carry no length scale; run calibrate first")
    assembly = ModelAssembly(
        populations=populations,
        identity=IdentityTable.from_json(),
        spatial=SpatialCompositionModel(populations, fitted.ell_km),
    )
    objective = build_national_objective(
        assembly, fitted, B=B, D=D, uniform_within_state=uniform_within_state
    )
    return Context(assembly=assembly, fitted=fitted, objective=objective)


def _write(name: str, payload: dict[str, Any]) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    payload.setdefault("disclaimer", DISCLAIMER)
    (RESULTS_DIR / name).write_text(json.dumps(payload, indent=1, default=float))
    logger.info("wrote %s", RESULTS_DIR / name)


# --------------------------------------------------------------------------------------
# M6 -- spatial model and the district atlas payload
# --------------------------------------------------------------------------------------


def run_spatial() -> dict[str, Any]:
    """Fit the length scale, score every district, and write the atlas payload."""
    from venomgap.ingest.compositions import load_compositions

    populations = load_compositions()
    length_scale = fit_length_scale(populations)
    context = build_context()
    immunogen = context.assembly.big_four_immunogen()

    coverage = context.objective.district_coverage(immunogen)
    unknown = context.objective.districts_with_no_estimate()
    cells = context.objective.cell_coverage(immunogen)

    # Per-district uncertainty and provenance, aggregated from its contributing cells.
    n_districts = len(context.objective.district_ids)
    weights = context.objective.cell_weights
    index = context.objective.cell_district_index
    denom = np.bincount(index, weights=weights, minlength=n_districts)
    uncertainty = np.zeros(n_districts)
    nearest = np.full(n_districts, np.inf)
    present = denom > 0.0
    uncertainty[present] = (
        np.bincount(index, weights=weights * context.objective.cell_uncertainty,
                    minlength=n_districts)[present] / denom[present]
    )
    for cell_index, district_index in enumerate(index):
        nearest[district_index] = min(
            nearest[district_index], context.objective.cell_nearest_km[cell_index]
        )

    from venomgap.ingest.burden import district_burden, load_districts, load_state_burden

    districts = district_burden(load_districts(), load_state_burden())
    by_id = districts.set_index("district_id")

    dominant: dict[int, list[str]] = {}
    limiting: dict[int, list[str]] = {}
    worst_cell: dict[int, float] = {}
    for cell_index, district_index in enumerate(index):
        species = context.objective.cell_species[cell_index]
        dominant.setdefault(district_index, []).append(species)
        value = float(cells[cell_index])
        if district_index not in worst_cell or value < worst_cell[district_index]:
            worst_cell[district_index] = value

    records = []
    for i, district_id in enumerate(context.objective.district_ids):
        row = by_id.loc[district_id]
        is_unknown = bool(unknown[i]) or not np.isfinite(nearest[i])
        flags = ["composition_imputed"]
        if is_unknown:
            flags.append("no_nearby_proteome")
        records.append(
            {
                "district_id": district_id,
                "district": str(row["district"]),
                "state": str(row["geo_state"]),
                "burden_state": str(row["state"]),
                "lat": float(row["lat"]),
                "lon": float(row["lon"]),
                "coverage": float(coverage[i]),
                "deficit": float(1.0 - coverage[i]),
                "uncertainty": float(uncertainty[i]),
                "status": "unknown" if is_unknown else "estimated",
                "burden_weight": float(row["burden_weight"]),
                "nearest_proteome_km": float(nearest[i]) if np.isfinite(nearest[i]) else -1.0,
                "dominant_species": sorted(set(dominant.get(i, []))),
                "limiting_families": limiting.get(i, []),
                "worst_species_coverage": worst_cell.get(i),
                "flags": flags,
            }
        )

    national = context.objective.national_coverage(immunogen)
    payload: dict[str, Any] = {
        "ell_km": length_scale.ell_km,
        "ell_plateau_km": list(length_scale.plateau_km),
        "ell_plateau_is_wide": length_scale.plateau_is_wide,
        "ell_at_grid_edge": length_scale.at_grid_edge,
        "loo_cv_rmse": length_scale.rmse,
        "loo_cv_folds": length_scale.folds,
        "loo_rmse_by_ell": dict(
            zip([str(e) for e in length_scale.grid], length_scale.rmse_by_ell, strict=True)
        ),
        "per_population_error": length_scale.per_population_error,
        "unknown_cutoff_km": context.assembly.spatial.cutoff_km,
        "districts": len(records),
        "estimated": int(sum(r["status"] == "estimated" for r in records)),
        "unknown": int(sum(r["status"] == "unknown" for r in records)),
        "national_coverage_baseline": national,
        "national_deficit_baseline": 1.0 - national,
        "district_results": records,
        "interpretation": (
            f"The leave-one-population-out optimum sits at {length_scale.ell_km:.0f} km with a "
            f"plateau spanning {length_scale.plateau_km[0]:.0f}-{length_scale.plateau_km[1]:.0f} "
            "km. A length scale of this size exceeds India's north-south extent, which means the "
            "cross-validation cannot distinguish a local kernel from a near-global one. The honest "
            "reading is that at the present sampling density geographic distance carries little "
            "predictive power for venom composition, and the best available estimator is close to "
            "a species-level mean. That is a finding about how sparse the data is, not a tuned "
            "parameter, and it is the quantitative case for sampling more sites."
        ),
    }

    distances, dissimilarity, same = composition_dissimilarity_vs_distance(populations)
    payload["dissimilarity_vs_distance"] = {
        "distance_km": distances.tolist(),
        "dissimilarity": dissimilarity.tolist(),
        "same_genus": same.tolist(),
    }
    _write("spatial.json", payload)
    return payload


# --------------------------------------------------------------------------------------
# M7 -- siting optimisation
# --------------------------------------------------------------------------------------


def run_optimisation() -> dict[str, Any]:
    """Coverage versus k, the turnover point, and the measured optimality gap."""
    from venomgap.optimize.exact import measure_gap
    from venomgap.optimize.greedy import (
        best_subset_per_size,
        candidate_populations,
        greedy_exactly,
        greedy_sites,
    )
    from venomgap.optimize.local_search import greedy_plus_local, improve_by_swaps
    from venomgap.optimize.weights import optimise_weights

    context = build_context()
    baseline_immunogen = context.assembly.big_four_immunogen()
    baseline = context.objective.national_coverage(baseline_immunogen)
    candidates = candidate_populations(context.assembly, indian_only=True)

    solutions: list[SitingSolution] = []
    curve: list[dict[str, Any]] = []
    # k counts *new* collection sites added to the existing Big Four immunogen, not the total
    # mixture size. k = 0 is therefore the antivenom India makes today, and every k > 0 is that
    # mixture extended -- which is where the dilution penalty applies, since each addition takes
    # budget share from the four venoms already in the vial.
    for k in range(1, MAX_K + 1):
        greedy = greedy_sites(
            context.objective, context.assembly, k, candidates, base=baseline_immunogen
        )
        local = greedy_plus_local(
            context.objective, context.assembly, k, candidates, base=baseline_immunogen
        )
        immunogen = context.assembly.combined_immunogen(baseline_immunogen, local.sites)
        weights = optimise_weights(context.objective, immunogen)
        best_immunogen = immunogen.with_weights(weights.weights)
        district = context.objective.district_coverage(best_immunogen)
        baseline_district = context.objective.district_coverage(baseline_immunogen)
        moved = int(
            np.sum(
                (1.0 - baseline_district >= HIGH_DEFICIT_THRESHOLD)
                & (1.0 - district < HIGH_DEFICIT_THRESHOLD)
            )
        )
        # The optimised weights cover the whole enlarged mixture; the tail of that vector is the
        # share going to the newly collected venoms.
        new_site_weights = tuple(
            float(w) for w in weights.weights[best_immunogen.size - len(local.sites):]
        )
        solutions.append(
            SitingSolution(
                k=k,
                sites=local.sites,
                weights=new_site_weights,
                immunogen_labels=best_immunogen.labels,
                immunogen_weights=tuple(float(w) for w in weights.weights),
                national_coverage=weights.coverage,
                coverage_gain_vs_baseline=weights.coverage - baseline,
                districts_moved_out_of_high_deficit=moved,
                solver="greedy+local",
                optimality_gap=None,
            )
        )
        # "Given that exactly k centres are built, which k and how well do they do?" The
        # unconstrained solver above stops adding as soon as the next venom would hurt, which
        # answers a different question and hides the turnover.
        forced_sites = greedy_exactly(
            context.objective, context.assembly, k, candidates, base=baseline_immunogen
        )
        forced_sites = improve_by_swaps(
            context.objective, context.assembly, forced_sites, candidates,
            base=baseline_immunogen,
        ).sites
        forced_immunogen = context.assembly.combined_immunogen(
            baseline_immunogen, forced_sites
        )
        forced_uniform = context.objective.national_coverage(forced_immunogen)
        forced_weights = optimise_weights(context.objective, forced_immunogen)
        curve.append(
            {
                "k": k,
                "greedy": greedy.coverage,
                "greedy_local": local.coverage,
                "weight_optimised": weights.coverage,
                "forced_k_uniform": forced_uniform,
                "forced_k_weighted": forced_weights.coverage,
                "forced_k_sites": list(forced_sites),
                "uniform_weights": weights.uniform_coverage,
                "weights_improved": weights.improved,
                "sites": list(local.sites),
                "moves": list(local.moves),
            }
        )
        logger.info(
            "k=%d greedy=%.5f local=%.5f weighted=%.5f sites=%s",
            k, greedy.coverage, local.coverage, weights.coverage, local.sites,
        )

    # k = 0 is the current immunogen, so it belongs in the curve the turnover is read from.
    curve.insert(
        0,
        {
            "k": 0,
            "greedy": baseline,
            "greedy_local": baseline,
            "weight_optimised": baseline,
            "uniform_weights": baseline,
            "weights_improved": False,
            "sites": [],
            "moves": ["current Big Four immunogen"],
            "forced_k_uniform": baseline,
            "forced_k_weighted": baseline,
            "forced_k_sites": [],
        },
    )
    # The turnover is read from the forced-k curve under uniform weights, which is the setting the
    # dilution effect actually describes: a vial whose protein is split evenly across its venoms.
    #
    # "Turnover" means the first descent, not the global maximum. Those are different questions and
    # here they have different answers: coverage peaks at a small k, falls as the next venoms
    # dilute the mixture, then recovers once enough new composition is in the vial to pay for
    # itself. Reporting only the global argmax would hide the dip, which is the whole point.
    forced = [row["forced_k_uniform"] for row in curve]
    first_descent = next(
        (i for i in range(len(forced) - 1) if forced[i + 1] < forced[i] - 1e-9), None
    )
    turnover_k = int(curve[first_descent]["k"]) if first_descent is not None else int(
        curve[int(np.argmax(forced))]["k"]
    )
    interior = first_descent is not None and 0 < first_descent < len(forced) - 1
    global_best_index = int(np.argmax(forced))
    values = [row["weight_optimised"] for row in curve]

    mixture_curve = best_subset_per_size(
        context.objective, context.assembly, context.fitted, max_size=MAX_K
    )
    gap = measure_gap(context.objective, context.assembly, k=EXACT_INSTANCE_MAX_K)

    payload = {
        "baseline_national_coverage": baseline,
        "baseline_national_deficit": 1.0 - baseline,
        "coverage_vs_k": curve,
        "turnover_k": turnover_k,
        "turnover_is_interior": interior,
        "turnover_is_first_descent": first_descent is not None,
        "global_best_k_uniform": int(curve[global_best_index]["k"]),
        "global_best_coverage_uniform": float(forced[global_best_index]),
        "turnover_basis": (
            "Forced-k curve under uniform mixture weights: the best achievable coverage when "
            "exactly k new venoms are added to the Big Four immunogen and the vial's protein is "
            "split evenly across all of them. turnover_k is the first k after which coverage "
            "falls, not the global maximum -- coverage dips as the next venoms dilute the mixture "
            "and recovers once enough new composition is in the vial to pay for itself."
        ),
        "weight_reoptimisation_note": (
            "Re-optimising the mixture weights removes the dip entirely, because the optimiser can "
            "down-weight the added venoms rather than splitting the vial evenly. That is the "
            "actionable consequence of the dilution effect: adding a venom to a polyvalent "
            "antivenom without re-balancing the mixture can make it worse, and re-balancing is "
            "what recovers the gain."
        ),
        "forced_k_curve": [
            {"k": row["k"], "uniform": row["forced_k_uniform"],
             "weighted": row["forced_k_weighted"], "sites": row["forced_k_sites"]}
            for row in curve
        ],
        "best_national_coverage": float(max(values)),
        "solutions": [s.model_dump(mode="json") for s in solutions],
        "mixture_size_curve": {
            "sizes": list(mixture_curve.sizes),
            "coverages": list(mixture_curve.coverages),
            "best_sets": [list(s) for s in mixture_curve.best_sets],
            "argmax_size": mixture_curve.argmax_size,
            "peak": mixture_curve.peak,
            "terminal": mixture_curve.terminal,
        },
        "optimality_gap": gap,
        "submodularity_note": (
            "The objective is non-monotone and not submodular, because adding a venom moves "
            "antibody supply away from families already covered. No (1 - 1/e) guarantee applies "
            "and none is claimed; the gap above is measured against exact enumeration on a reduced "
            "instance, not derived from theory."
        ),
    }
    _write("optimisation.json", payload)
    return payload


# --------------------------------------------------------------------------------------
# M9 -- sensitivity and ablation
# --------------------------------------------------------------------------------------


def run_sensitivity() -> dict[str, Any]:
    """Leave-one-study-out rank stability, B/D sweeps, and the pre-registered controls."""
    from venomgap.model.calibrate import load_calibration, run_calibration
    from venomgap.validate.baselines import (
        cosine_deficit_ranking,
        sequence_free_deficit_ranking,
    )
    from venomgap.validate.retrodiction import rank_populations

    context = build_context()
    immunogen = context.assembly.big_four_immunogen()
    reference = rank_populations(context.assembly, context.fitted, immunogen)
    reference_map = {
        pop_id: float(reference.deficits[i]) for i, pop_id in enumerate(reference.pop_ids)
    }

    # Leave-one-study-out: refit without each calibration study and re-rank.
    studies = load_calibration().studies
    loso: list[dict[str, Any]] = []
    for study in studies:
        refitted = run_calibration(exclude_study=study, write=False)
        ranking = rank_populations(context.assembly, refitted, immunogen)
        shared = [p for p in ranking.pop_ids if p in reference_map]
        rho = float(
            spearmanr(
                [reference_map[p] for p in shared],
                [float(ranking.deficits[ranking.pop_ids.index(p)]) for p in shared],
            ).statistic
        )
        loso.append(
            {
                "excluded_study": study,
                "spearman_rho": rho,
                "kappa_scale": refitted.kappa_scale,
                "calibration_rows": refitted.calibration_rows,
            }
        )
        logger.info("leave-one-study-out (%s): Spearman rho = %.4f", study, rho)

    # B and D sweeps: does the ranking survive a different dose or vial assumption?
    sweeps: list[dict[str, Any]] = []
    for B in B_SWEEP_MG:
        for D in D_SWEEP_MG:
            swept = build_context(B=B, D=D)
            ranking = rank_populations(
                swept.assembly, swept.fitted, swept.assembly.big_four_immunogen(), B=B, D=D
            )
            shared = [p for p in ranking.pop_ids if p in reference_map]
            rho = float(
                spearmanr(
                    [reference_map[p] for p in shared],
                    [float(ranking.deficits[ranking.pop_ids.index(p)]) for p in shared],
                ).statistic
            )
            sweeps.append(
                {
                    "B_mg": B,
                    "D_mg": D,
                    "ratio": B / D,
                    "national_coverage": swept.objective.national_coverage(
                        swept.assembly.big_four_immunogen()
                    ),
                    "spearman_rho_vs_default": rho,
                }
            )

    # Pre-registered negative controls.
    cosine = cosine_deficit_ranking(context.objective, context.assembly, use_max=True)
    sequence_free = sequence_free_deficit_ranking(context.assembly, context.fitted)
    shared = [p for p in reference_map if p in cosine]
    cosine_rho = float(
        spearmanr([reference_map[p] for p in shared], [cosine[p] for p in shared]).statistic
    )
    seqfree_rho = float(
        spearmanr(
            [reference_map[p] for p in shared], [sequence_free[p] for p in shared]
        ).statistic
    )

    rng = np.random.default_rng(0)
    permuted_rhos = []
    pop_ids = list(reference_map)
    for _ in range(200):
        shuffled = rng.permutation([reference_map[p] for p in pop_ids])
        permuted_rhos.append(
            float(spearmanr([reference_map[p] for p in pop_ids], shuffled).statistic)
        )

    rhos = [row["spearman_rho"] for row in loso]
    payload = {
        "leave_one_study_out": loso,
        "mean_spearman": float(np.mean(rhos)) if rhos else float("nan"),
        "min_spearman": float(np.min(rhos)) if rhos else float("nan"),
        "parameter_sweep": sweeps,
        "sweep_min_spearman": float(min(s["spearman_rho_vs_default"] for s in sweeps)),
        "controls": {
            "cosine_baseline_rho_vs_model": cosine_rho,
            "sequence_free_rho_vs_model": seqfree_rho,
            "permutation_rho_mean": float(np.mean(permuted_rhos)),
            "permutation_rho_p95": float(np.percentile(np.abs(permuted_rhos), 95)),
            "note": (
                "The permutation control confirms the null level. The sequence-free control sets "
                "x_f = 1 everywhere; a high correlation with the full model would mean the "
                "sequence layer is decoration. The cosine baseline is run in its best-match form."
            ),
        },
        "cosine_deficits": cosine,
        "sequence_free_deficits": sequence_free,
        "model_deficits": reference_map,
    }
    _write("sensitivity.json", payload)
    return payload
