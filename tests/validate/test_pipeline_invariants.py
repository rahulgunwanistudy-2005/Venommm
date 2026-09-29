"""Invariants the whole pipeline must hold, checked against the committed results.

These are the tests that catch a regression in the *claims*, not in the code. They read
`results/*.json` and assert the things the README, the figures and the demo all depend on.
"""

from __future__ import annotations

import json
from itertools import pairwise

import pytest

from venomgap.config import RESULTS_DIR


def _load(name: str) -> dict:
    path = RESULTS_DIR / name
    if not path.exists():
        pytest.skip(f"{name} not generated; run `python -m venomgap.cli pipeline`")
    return json.loads(path.read_text())


@pytest.fixture(scope="module")
def retro() -> dict:
    return _load("retrodiction.json")


@pytest.fixture(scope="module")
def spatial() -> dict:
    return _load("spatial.json")


@pytest.fixture(scope="module")
def optimisation() -> dict:
    return _load("optimisation.json")


@pytest.fixture(scope="module")
def sensitivity() -> dict:
    return _load("sensitivity.json")


# --------------------------------------------------------------------------------------
# The headline claim: the model produces the dilution effect and the baseline cannot
# --------------------------------------------------------------------------------------


def test_the_model_turns_over_in_mixture_size(retro: dict) -> None:
    coverages = retro["mixture_size_curve"]["coverages"]
    peak = coverages.index(max(coverages))
    assert 0 < peak < len(coverages) - 1, (
        f"the mixture-size curve has no interior maximum: {coverages}"
    )


def test_the_cosine_baseline_cannot_turn_over(retro: dict) -> None:
    """R4c. If this ever fails, the ablation no longer distinguishes the two models."""
    coverages = retro["cosine_baseline_curve"]["coverages"]
    peak = coverages.index(max(coverages))
    assert peak == len(coverages) - 1, (
        f"the cosine baseline produced an interior maximum: {coverages}"
    )
    assert all(
        later >= earlier - 1e-9
        for earlier, later in pairwise(coverages)
    ), "the best-match cosine baseline must be monotone non-decreasing in mixture size"


def test_baseline_is_run_in_its_strongest_form(retro: dict) -> None:
    """A straw-man ablation proves nothing. The comparison must use the best-match form."""
    assert "best-match" in retro["cosine_baseline_curve"]["form"]


# --------------------------------------------------------------------------------------
# Reporting integrity
# --------------------------------------------------------------------------------------


def test_every_criterion_records_its_threshold_and_outcome(retro: dict) -> None:
    for criterion in retro["criteria"]:
        assert criterion["statement"]
        assert isinstance(criterion["passed"], bool)
        assert criterion["comparison"] in {">=", "<=", "==", "bool", "ill_posed"}
        if criterion["comparison"] in {">=", "<="}:
            assert criterion["threshold"] is not None
            assert criterion["computed"] is not None


def test_a_target_passes_only_if_all_its_subcriteria_pass(retro: dict) -> None:
    groups: dict[str, list[bool]] = {}
    for criterion in retro["criteria"]:
        groups.setdefault(criterion["criterion_id"][:2], []).append(criterion["passed"])
    expected = sum(1 for results in groups.values() if all(results))
    assert retro["targets_passed"] == expected


def test_the_preregistration_commit_is_recorded(retro: dict) -> None:
    commit = retro["preregistration_commit"]
    assert commit != "unknown" and len(commit) == 40, (
        "the retrodiction must record the commit that introduced the pre-registration"
    )


def test_failures_are_diagnosed_not_just_reported(retro: dict) -> None:
    """A miss without a diagnosis is an anecdote. Every failure needs computed evidence."""
    if retro["targets_passed"] < retro["targets_total"]:
        diagnostics = retro["diagnostics"]
        assert diagnostics["theta_scale_mismatch"]["fraction_below_theta"] is not None
        assert diagnostics["theta_scale_mismatch"]["consequence"]
        assert diagnostics["data_sparsity"]["populations_per_species"]


# --------------------------------------------------------------------------------------
# Unknown districts are never coloured as covered
# --------------------------------------------------------------------------------------


def test_unknown_districts_carry_the_flag(spatial: dict) -> None:
    for district in spatial["district_results"]:
        if district["status"] == "unknown":
            assert "no_nearby_proteome" in district["flags"], (
                f"{district['district_id']} is unknown but does not say why"
            )


def test_district_coverage_and_deficit_are_complementary(spatial: dict) -> None:
    for district in spatial["district_results"]:
        assert district["coverage"] + district["deficit"] == pytest.approx(1.0, abs=1e-9)
        assert 0.0 <= district["coverage"] <= 1.0


def test_the_length_scale_plateau_is_reported(spatial: dict) -> None:
    """A flat-bottomed CV curve quoted as a point estimate would overstate what was fitted."""
    low, high = spatial["ell_plateau_km"]
    assert low <= spatial["ell_km"] <= high
    assert "plateau" in spatial["interpretation"].lower()


# --------------------------------------------------------------------------------------
# The optimiser makes no guarantee it does not have
# --------------------------------------------------------------------------------------


def test_no_submodularity_guarantee_is_claimed(optimisation: dict) -> None:
    note = optimisation["submodularity_note"]
    assert "not submodular" in note
    assert "1 - 1/e" in note or "1 − 1/e" in note


def test_the_optimality_gap_is_measured(optimisation: dict) -> None:
    gap = optimisation["optimality_gap"]
    assert gap["subsets_evaluated"] > 0
    assert 0.0 <= gap["optimality_gap"] <= 1.0
    assert gap["heuristic_coverage"] <= gap["exact_coverage"] + 1e-9


def test_k_zero_is_the_current_antivenom(optimisation: dict) -> None:
    """k counts *new* sites, so k = 0 must reproduce the baseline exactly."""
    zero = next(row for row in optimisation["coverage_vs_k"] if row["k"] == 0)
    assert zero["weight_optimised"] == pytest.approx(
        optimisation["baseline_national_coverage"]
    )
    assert zero["sites"] == []


def test_adding_sites_never_reports_a_worse_unconstrained_solution(
    optimisation: dict,
) -> None:
    """The solver is free to add fewer than k, so its curve must not decrease."""
    values = [row["weight_optimised"] for row in optimisation["coverage_vs_k"]]
    assert all(
        later >= earlier - 1e-9 for earlier, later in pairwise(values)
    ), f"the unconstrained solver returned a worse portfolio at a larger budget: {values}"


# --------------------------------------------------------------------------------------
# Sensitivity: the pre-registered falsification thresholds
# --------------------------------------------------------------------------------------


def test_rank_stability_clears_the_preregistered_threshold(sensitivity: dict) -> None:
    """Section 7 of the pre-registration: rho below 0.5 falsifies the ranking."""
    assert sensitivity["min_spearman"] >= 0.5, (
        f"leave-one-study-out rho fell to {sensitivity['min_spearman']:.3f}, which the "
        f"pre-registration lists as a falsification condition"
    )


def test_the_sequence_layer_is_doing_work(sensitivity: dict) -> None:
    """The pre-registered sequence-free control. If x_f = 1 ranked as well as the real model,
    the cross-reactivity layer would be decoration."""
    control = sensitivity["controls"]
    assert abs(control["sequence_free_rho_vs_model"]) < 0.9


def test_the_cosine_baseline_ranks_differently(sensitivity: dict) -> None:
    control = sensitivity["controls"]
    assert control["cosine_baseline_rho_vs_model"] < 0.9


def test_the_permutation_null_is_near_zero(sensitivity: dict) -> None:
    control = sensitivity["controls"]
    assert abs(control["permutation_rho_mean"]) < 0.2
