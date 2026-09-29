"""API routes. Three screens' worth of data plus a live mixture evaluator."""

from __future__ import annotations

import json
import logging
from functools import lru_cache
from typing import Any

import numpy as np
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field

from venomgap.config import DISCLAIMER, FAMILIES, RESULTS_DIR
from venomgap.errors import VenomGapError

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api")


def _result(name: str) -> dict[str, Any]:
    path = RESULTS_DIR / name
    if not path.exists():
        raise HTTPException(
            status_code=503,
            detail=f"{name} has not been generated; run `python -m venomgap.cli pipeline`",
        )
    return json.loads(path.read_text())  # type: ignore[no-any-return]


@lru_cache(maxsize=1)
def _context():  # type: ignore[no-untyped-def]
    """The live model, built once. Only the mixture explorer needs it."""
    from venomgap.pipeline import build_context

    return build_context()


# --------------------------------------------------------------------------------------
# Screen 1 -- the atlas
# --------------------------------------------------------------------------------------


@router.get("/atlas")
def atlas() -> dict[str, Any]:
    """District choropleth payload: deficit, uncertainty, burden, status."""
    spatial = _result("spatial.json")
    return {
        "districts": spatial["district_results"],
        "national_coverage": spatial["national_coverage_baseline"],
        "national_deficit": spatial["national_deficit_baseline"],
        "ell_km": spatial["ell_km"],
        "ell_plateau_km": spatial["ell_plateau_km"],
        "unknown_cutoff_km": spatial["unknown_cutoff_km"],
        "estimated": spatial["estimated"],
        "unknown": spatial["unknown"],
        "interpretation": spatial["interpretation"],
        "disclaimer": DISCLAIMER,
    }


@router.get("/siting")
def siting() -> dict[str, Any]:
    """Cached siting solutions for k = 0..6, which the atlas slider moves through."""
    opt = _result("optimisation.json")
    context = _context()
    populations = {p.pop_id: p for p in context.assembly.populations}

    def pin(pop_id: str, weight: float) -> dict[str, Any] | None:
        population = populations.get(pop_id)
        if population is None:
            return None
        return {
            "pop_id": pop_id,
            "species": population.species,
            "locality": population.locality,
            "state": population.state,
            "lat": population.lat,
            "lon": population.lon,
            "weight": weight,
        }

    forced = []
    for row in opt.get("forced_k_curve", []):
        sites = row["sites"]
        share = 1.0 / (len(sites) + 4) if sites else 0.0
        forced.append(
            {
                **row,
                "pins": [p for p in (pin(s, share) for s in sites) if p is not None],
            }
        )

    solutions = []
    for solution in opt["solutions"]:
        pins = [
            p
            for p in (
                pin(pop_id, weight)
                for pop_id, weight in zip(
                    solution["sites"], solution["weights"], strict=False
                )
            )
            if p is not None
        ]
        solutions.append({**solution, "pins": pins})
    return {
        "baseline_national_coverage": opt["baseline_national_coverage"],
        "solutions": solutions,
        "coverage_vs_k": opt["coverage_vs_k"],
        "forced_k_curve": forced,
        "turnover_k": opt["turnover_k"],
        "turnover_basis": opt.get("turnover_basis", ""),
        "weight_reoptimisation_note": opt.get("weight_reoptimisation_note", ""),
        "global_best_k_uniform": opt.get("global_best_k_uniform"),
        "optimality_gap": opt["optimality_gap"],
        "submodularity_note": opt["submodularity_note"],
        "disclaimer": DISCLAIMER,
    }


# --------------------------------------------------------------------------------------
# Screen 2 -- district detail
# --------------------------------------------------------------------------------------


@router.get("/district/{district_id}")
def district(district_id: str) -> dict[str, Any]:
    """Everything the detail panel shows, including provenance for every contributing proteome."""
    spatial = _result("spatial.json")
    match = next(
        (d for d in spatial["district_results"] if d["district_id"] == district_id), None
    )
    if match is None:
        raise HTTPException(status_code=404, detail=f"unknown district {district_id}")

    context = _context()
    immunogen = context.assembly.big_four_immunogen()
    theta = context.fitted.theta_vector()

    per_species: list[dict[str, Any]] = []
    provenance: list[dict[str, Any]] = []
    seen: set[str] = set()
    for index, cell_district in enumerate(context.objective.cell_district_index):
        if context.objective.district_ids[cell_district] != district_id:
            continue
        species = context.objective.cell_species[index]
        composition = context.objective.cell_compositions[index]
        inputs, imputed = context.assembly.inputs_for(
            composition, species, immunogen, theta, context.fitted.kappa_scale
        )
        from venomgap.model.coverage import coverage, limiting_families, neutralised_fractions

        B, D = context.objective.B, context.objective.D
        n = neutralised_fractions(inputs, immunogen.weights, B=B, D=D)
        value = coverage(inputs, immunogen.weights, B=B, D=D)
        per_species.append(
            {
                "species": species,
                "coverage": value,
                "deficit": 1.0 - value,
                "bite_share": float(context.objective.cell_weights[index]),
                "uncertainty": float(context.objective.cell_uncertainty[index]),
                "nearest_proteome_km": float(context.objective.cell_nearest_km[index]),
                "status": context.objective.cell_status[index],
                "composition": {FAMILIES[i]: float(composition[i]) for i in range(len(FAMILIES))},
                "per_family_neutralised": {FAMILIES[i]: float(n[i]) for i in range(len(FAMILIES))},
                "limiting_families": list(limiting_families(inputs, n)),
                "sequence_imputed": bool(imputed),
            }
        )
        estimate = context.assembly.spatial.estimate(match["lat"], match["lon"], species)
        for contributor in estimate.contributors:
            if contributor in seen:
                continue
            seen.add(contributor)
            population = context.assembly.by_id[contributor]
            provenance.append(
                {
                    "pop_id": contributor,
                    "species": population.species,
                    "locality": population.locality,
                    "state": population.state,
                    "doi": population.provenance.doi,
                    "table": population.provenance.table,
                    "method": population.provenance.method,
                    "accessed": population.provenance.accessed.isoformat(),
                    "flags": list(population.flags),
                    "note": population.provenance.note,
                }
            )

    per_species.sort(key=lambda row: -row["deficit"])
    return {
        **match,
        "per_species": per_species,
        "provenance": provenance,
        "flag_explanations": FLAG_EXPLANATIONS,
        "disclaimer": DISCLAIMER,
    }


FLAG_EXPLANATIONS: dict[str, str] = {
    "composition_imputed": (
        "No venom proteome has been published for this exact location. The composition is "
        "interpolated from the nearest sampled populations of the same species."
    ),
    "no_nearby_proteome": (
        "The nearest published proteome is further away than this model is willing to extrapolate. "
        "The district is reported as unknown, not as covered."
    ),
    "sequence_imputed": (
        "Toxin sequences for this species and family are missing or very few, so cross-reactivity "
        "falls back to a genus-level or global estimate."
    ),
    "single_study_basis": "Everything known about this population comes from one study.",
    "partial_table_renormalised": (
        "The source study characterised only part of whole venom. The unreported remainder is "
        "assigned to 'other' rather than redistributed among the families that were measured."
    ),
    "genus_consensus_sequences": "Sequences were borrowed from congeners of the same genus.",
}


# --------------------------------------------------------------------------------------
# Screen 3 -- the mixture explorer
# --------------------------------------------------------------------------------------


class MixtureRequest(BaseModel):
    pop_ids: list[str] = Field(min_length=1, max_length=12)
    weights: list[float] | None = None
    B_mg: float | None = Field(default=None, gt=0.0)
    D_mg: float | None = Field(default=None, gt=0.0)


@router.get("/populations")
def populations() -> dict[str, Any]:
    """The candidate venoms the mixture explorer offers."""
    context = _context()
    rows = []
    for population in context.assembly.populations:
        dominant = max(population.composition, key=lambda k: population.composition[k])
        rows.append(
            {
                "pop_id": population.pop_id,
                "species": population.species,
                "locality": population.locality,
                "state": population.state,
                "country": population.country,
                "lat": population.lat,
                "lon": population.lon,
                "dominant_family": dominant,
                "composition": population.composition,
                "doi": population.provenance.doi,
                "table": population.provenance.table,
                "flags": list(population.flags),
                "holdout": population.holdout,
            }
        )
    rows.sort(key=lambda r: (r["species"], r["state"]))
    return {"populations": rows, "families": list(FAMILIES), "disclaimer": DISCLAIMER}


@router.post("/mixture")
def mixture(request: MixtureRequest) -> dict[str, Any]:
    """Evaluate an arbitrary immunogen mixture live. This is the interactive proof of R4."""
    context = _context()
    unknown = [p for p in request.pop_ids if p not in context.assembly.by_id]
    if unknown:
        raise HTTPException(status_code=400, detail=f"unknown populations: {unknown}")

    weights = None
    if request.weights is not None:
        if len(request.weights) != len(request.pop_ids):
            raise HTTPException(status_code=400, detail="weights and pop_ids differ in length")
        total = sum(request.weights)
        if total <= 0:
            raise HTTPException(status_code=400, detail="weights must sum to a positive number")
        weights = np.array([w / total for w in request.weights], dtype=float)

    try:
        immunogen = context.assembly.immunogen_from_pop_ids(tuple(request.pop_ids), weights)
    except VenomGapError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc

    objective = context.objective
    if request.B_mg is not None or request.D_mg is not None:
        from venomgap.pipeline import build_context

        objective = build_context(
            B=request.B_mg or objective.B, D=request.D_mg or objective.D
        ).objective

    district = objective.district_coverage(immunogen)
    national = float(np.dot(objective.burden, district))
    baseline = objective.national_coverage(context.assembly.big_four_immunogen())

    return {
        "pop_ids": request.pop_ids,
        "weights": [float(w) for w in immunogen.weights],
        "national_coverage": national,
        "national_deficit": 1.0 - national,
        "baseline_national_coverage": baseline,
        "change_vs_baseline": national - baseline,
        "mixture_size": immunogen.size,
        "district_coverage": {
            objective.district_ids[i]: float(district[i]) for i in range(len(district))
        },
        "disclaimer": DISCLAIMER,
    }


@router.get("/retrodiction")
def retrodiction() -> dict[str, Any]:
    payload = _result("retrodiction.json")
    return {
        "targets_passed": payload["targets_passed"],
        "targets_total": payload["targets_total"],
        "criteria": payload["criteria"],
        "preregistration_commit": payload["preregistration_commit"],
        "notes": payload["notes"],
        "diagnostics": payload["diagnostics"],
        "mixture_size_curve": payload["mixture_size_curve"],
        "cosine_baseline_curve": payload["cosine_baseline_curve"],
        "deficit_ranking": payload["deficit_ranking"],
        "disclaimer": DISCLAIMER,
    }
