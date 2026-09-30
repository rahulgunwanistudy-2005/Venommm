"""Export everything the browser needs to run the coverage model with no backend.

The model is arithmetic over small fixed tables, so it ports to JavaScript exactly. What the
browser cannot do is the sequence alignment and the spatial interpolation, so both are precomputed
here: cross-reactivity collapses to a (species x species x family) table because x_f depends only
on the species pair, and district compositions are already interpolated.

This makes the hosted atlas a static page that reproduces the Python results bit for bit.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from venomgap.config import DEFAULT_B_MG, DEFAULT_D_MG, FAMILIES, IRULA_LAT, IRULA_LON
from venomgap.model.coverage import SIGMA_VECTOR
from venomgap.model.crossreact import crossreact_matrix
from venomgap.pipeline import build_context

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "web" / "public" / "venomgap-bundle.json"

R = 5  # decimal places; composition values are relative abundances, 5 dp is far past precision


def r(x: float) -> float:
    return round(float(x), R)


def main() -> None:
    context = build_context()
    assembly, fitted, objective = context.assembly, context.fitted, context.objective
    theta = fitted.theta_vector()
    kappa = fitted.kappa_vector()

    species = sorted({p.species for p in assembly.populations} | set(objective.cell_species))

    # x_f(q -> p) depends only on the species pair, so a 9 x 9 x 16 table is complete.
    cross: dict[str, list[float]] = {}
    for q in species:
        for p in species:
            matrix, _ = crossreact_matrix([q], p, assembly.identity, theta, assembly.genus_of)
            cross[f"{q}|{p}"] = [r(v) for v in matrix[0]]

    base = assembly.big_four_immunogen()
    spatial = json.loads((ROOT / "results" / "spatial.json").read_text())
    opt = json.loads((ROOT / "results" / "optimisation.json").read_text())
    retro = json.loads((ROOT / "results" / "retrodiction.json").read_text())
    sens = json.loads((ROOT / "results" / "sensitivity.json").read_text())

    district_by_id = {d["district_id"]: d for d in spatial["district_results"]}

    bundle = {
        "meta": {
            "generated_from": "results/*.json + the live model",
            "preregistration_commit": retro["preregistration_commit"],
            "data_version": retro["data_version"],
            "disclaimer": (
                "Research model. Not clinical guidance. Predictions are computational and "
                "require experimental validation."
            ),
        },
        "families": list(FAMILIES),
        "sigma": [r(v) for v in SIGMA_VECTOR],
        "kappa": [r(v) for v in kappa],
        "B": DEFAULT_B_MG,
        "D": DEFAULT_D_MG,
        "kappaScale": r(fitted.kappaScale) if hasattr(fitted, "kappaScale") else r(fitted.kappa_scale),
        "crossreact": cross,
        "immunogenSource": {"lat": IRULA_LAT, "lon": IRULA_LON},
        "baseImmunogen": {
            "labels": list(base.labels),
            "species": list(base.species),
            "compositions": [[r(v) for v in row] for row in base.compositions],
            "weights": [r(v) for v in base.weights],
        },
        "populations": [
            {
                "popId": p.pop_id,
                "species": p.species,
                "locality": p.locality,
                "state": p.state,
                "country": p.country,
                "lat": p.lat,
                "lon": p.lon,
                "composition": [r(v) for v in p.vector()],
                "dominant": max(p.composition, key=lambda k: p.composition[k]),
                "doi": p.provenance.doi,
                "table": p.provenance.table,
                "method": p.provenance.method,
                "accessed": p.provenance.accessed.isoformat(),
                "flags": list(p.flags),
                "holdout": p.holdout,
                "note": p.provenance.note,
            }
            for p in assembly.populations
        ],
        "cells": {
            "districtIndex": [int(i) for i in objective.cell_district_index],
            "species": list(objective.cell_species),
            "compositions": [[r(v) for v in row] for row in objective.cell_compositions],
            "biteShare": [r(v) for v in objective.cell_weights],
            "uncertainty": [r(v) for v in objective.cell_uncertainty],
            "nearestKm": [round(float(v), 1) for v in objective.cell_nearest_km],
            "status": list(objective.cell_status),
        },
        "districts": [
            {
                "id": district_id,
                "name": district_by_id[district_id]["district"],
                "state": district_by_id[district_id]["state"],
                "lat": district_by_id[district_id]["lat"],
                "lon": district_by_id[district_id]["lon"],
                "burden": r(float(objective.burden[i])),
                "coverage": r(district_by_id[district_id]["coverage"]),
                "deficit": r(district_by_id[district_id]["deficit"]),
                "uncertainty": r(district_by_id[district_id]["uncertainty"]),
                "status": district_by_id[district_id]["status"],
                "nearestKm": round(district_by_id[district_id]["nearest_proteome_km"], 1),
                "species": district_by_id[district_id]["dominant_species"],
                "flags": district_by_id[district_id]["flags"],
            }
            for i, district_id in enumerate(objective.district_ids)
        ],
        "results": {
            "spatial": {
                k: spatial[k]
                for k in (
                    "ell_km", "ell_plateau_km", "loo_cv_rmse", "loo_cv_folds",
                    "loo_rmse_by_ell", "unknown_cutoff_km", "estimated", "unknown",
                    "national_coverage_baseline", "national_deficit_baseline",
                    "interpretation",
                )
            },
            "optimisation": {
                k: opt[k]
                for k in (
                    "baseline_national_coverage", "forced_k_curve", "coverage_vs_k",
                    "turnover_k", "turnover_basis", "weight_reoptimisation_note",
                    "global_best_k_uniform", "optimality_gap", "submodularity_note",
                )
            },
            "retrodiction": {
                "targets_passed": retro["targets_passed"],
                "targets_total": retro["targets_total"],
                "criteria": retro["criteria"],
                "notes": retro["notes"],
                "diagnostics": retro["diagnostics"],
                "mixture_size_curve": retro["mixture_size_curve"],
                "cosine_baseline_curve": retro["cosine_baseline_curve"],
                "deficit_ranking": retro["deficit_ranking"],
                "fitted_parameters": {
                    k: retro["fitted_parameters"][k]
                    for k in (
                        "kappa_scale", "c_star", "theta", "calibration_rows",
                        "calibration_rmse", "calibration_studies", "theta_unidentified",
                        "ordinal_checks", "ordinal_violations", "ell_km",
                    )
                },
            },
            "sensitivity": {
                "mean_spearman": sens["mean_spearman"],
                "min_spearman": sens["min_spearman"],
                "sweep_min_spearman": sens["sweep_min_spearman"],
                "leave_one_study_out": sens["leave_one_study_out"],
                "controls": sens["controls"],
                "model_deficits": {k: r(v) for k, v in sens["model_deficits"].items()},
                "cosine_deficits": {k: r(v) for k, v in sens["cosine_deficits"].items()},
                "sequence_free_deficits": {
                    k: r(v) for k, v in sens["sequence_free_deficits"].items()
                },
            },
        },
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(bundle, separators=(",", ":")))
    size_kb = OUT.stat().st_size / 1024
    print(f"wrote {OUT} ({size_kb:.0f} KB)")
    print(f"  districts {len(bundle['districts'])}, cells {len(bundle['cells']['species'])}, "
          f"populations {len(bundle['populations'])}, crossreact pairs {len(cross)}")

    # Verify the exported tables reproduce the Python baseline exactly.
    comps = np.array(bundle["cells"]["compositions"])
    sigma = np.array(bundle["sigma"])
    kap = np.array(bundle["kappa"])
    w = np.array(bundle["baseImmunogen"]["weights"])
    imm = np.array(bundle["baseImmunogen"]["compositions"])
    burden = np.array([d["burden"] for d in bundle["districts"]])
    idx = np.array(bundle["cells"]["districtIndex"])
    share = np.array(bundle["cells"]["biteShare"])
    x = np.array([
        [cross[f"{s}|{t}"] for s in bundle["baseImmunogen"]["species"]]
        for t in bundle["cells"]["species"]
    ])
    eff = bundle["B"] * np.einsum("q,tqf->tf", w, imm[None, :, :] * x)
    demand = bundle["D"] * comps * kap[None, :]
    n = np.ones_like(demand)
    present = demand > 1e-12
    np.divide(eff, demand, out=n, where=present)
    np.minimum(n, 1.0, out=n)
    n[~present] = 1.0
    weight = sigma[None, :] * comps
    cell_cov = (weight * n).sum(axis=1) / weight.sum(axis=1)
    num = np.bincount(idx, weights=cell_cov * share, minlength=len(burden))
    den = np.bincount(idx, weights=share, minlength=len(burden))
    dist = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
    national = float(np.dot(burden / burden.sum(), dist))
    expected = spatial["national_coverage_baseline"]
    print(f"  round-trip national coverage {national:.6f} vs python {expected:.6f} "
          f"(delta {abs(national - expected):.2e})")
    assert abs(national - expected) < 2e-4, "exported tables do not reproduce the Python result"


if __name__ == "__main__":
    main()
