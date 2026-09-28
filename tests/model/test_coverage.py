"""M1 gate: analytic properties of the coverage model on synthetic populations.

Every test here checks a closed-form consequence of the model equations, not that the code runs.
If any of these fail, nothing downstream in VenomGap is meaningful.
"""

from __future__ import annotations

from itertools import pairwise

import numpy as np
import pytest

from venomgap.config import FAMILIES, FAMILY_INDEX, SIGMA
from venomgap.errors import DataValidationError
from venomgap.model.coverage import (
    CoverageInputs,
    antibody_demand,
    antibody_supply,
    coverage,
    coverage_many,
    effective_supply,
    evaluate,
    neutralised_fractions,
)

N = len(FAMILIES)


def comp(**kwargs: float) -> np.ndarray:
    """Build a composition vector from family=abundance keywords, asserting it sums to 1."""
    v = np.zeros(N)
    for family, value in kwargs.items():
        v[FAMILY_INDEX[family]] = value
    assert abs(v.sum() - 1.0) < 1e-12, f"synthetic composition sums to {v.sum()}"
    return v


def unit_kappa() -> np.ndarray:
    return np.ones(N)


def perfect_crossreact(n_immunogen: int) -> np.ndarray:
    return np.ones((n_immunogen, N))


# --------------------------------------------------------------------------------------
# Property 1 — a perfectly matched, amply dosed antivenom reaches coverage 1
# --------------------------------------------------------------------------------------


def test_identical_population_and_immunogen_gives_full_coverage() -> None:
    p = comp(**{"3FTx": 0.6, "PLA2": 0.3, "SVMP": 0.1})
    inputs = CoverageInputs(
        target_composition=p,
        immunogen_compositions=p[None, :],
        crossreact=perfect_crossreact(1),
        kappa=unit_kappa(),
    )
    # B == D with unit kappa makes supply exactly equal demand in every family.
    assert coverage(inputs, np.array([1.0]), B=100.0, D=100.0) == pytest.approx(1.0)


def test_exactly_sufficient_budget_is_the_boundary_case() -> None:
    """With B = D and unit kappa, supply equals demand family by family: n_f == 1 everywhere."""
    p = comp(**{"3FTx": 0.5, "SVMP": 0.5})
    inputs = CoverageInputs(
        target_composition=p,
        immunogen_compositions=p[None, :],
        crossreact=perfect_crossreact(1),
        kappa=unit_kappa(),
    )
    n = neutralised_fractions(inputs, np.array([1.0]), B=40.0, D=40.0)
    assert np.allclose(n, 1.0)
    # One percent less budget and both families fall short by exactly one percent.
    n_short = neutralised_fractions(inputs, np.array([1.0]), B=39.6, D=40.0)
    assert n_short[FAMILY_INDEX["3FTx"]] == pytest.approx(0.99)
    assert n_short[FAMILY_INDEX["SVMP"]] == pytest.approx(0.99)


# --------------------------------------------------------------------------------------
# Property 2 — disjoint family support gives coverage ~ 0
# --------------------------------------------------------------------------------------


def test_disjoint_families_give_zero_coverage() -> None:
    target = comp(**{"3FTx": 1.0})
    immunogen = comp(**{"SVMP": 1.0})
    inputs = CoverageInputs(
        target_composition=target,
        immunogen_compositions=immunogen[None, :],
        crossreact=perfect_crossreact(1),
        kappa=unit_kappa(),
    )
    assert coverage(inputs, np.array([1.0]), B=10_000.0, D=1.0) == pytest.approx(0.0)


def test_zero_crossreactivity_gives_zero_coverage_even_with_identical_composition() -> None:
    """Composition match is not recognition. This is what separates the model from cosine."""
    p = comp(**{"3FTx": 0.7, "PLA2": 0.3})
    inputs = CoverageInputs(
        target_composition=p,
        immunogen_compositions=p[None, :],
        crossreact=np.zeros((1, N)),
        kappa=unit_kappa(),
    )
    assert coverage(inputs, np.array([1.0]), B=10_000.0, D=1.0) == pytest.approx(0.0)


# --------------------------------------------------------------------------------------
# Property 3 — THE DILUTION EFFECT. Adding a venom can lower coverage. (R4 mechanism.)
# --------------------------------------------------------------------------------------


def test_adding_a_dilutant_venom_lowers_coverage() -> None:
    """The R4 mechanism, in its purest synthetic form.

    The target is pure 3FTx. Immunogen A is the same. Immunogen B is pure SVMP and shares nothing
    with the target. Adding B to the mixture halves the antibody raised against 3FTx, so coverage
    must fall. No parameter was fitted to produce this.
    """
    target = comp(**{"3FTx": 1.0})
    a = comp(**{"3FTx": 1.0})
    b = comp(**{"SVMP": 1.0})

    single = CoverageInputs(
        target_composition=target,
        immunogen_compositions=a[None, :],
        crossreact=perfect_crossreact(1),
        kappa=unit_kappa(),
    )
    pair = CoverageInputs(
        target_composition=target,
        immunogen_compositions=np.vstack([a, b]),
        crossreact=perfect_crossreact(2),
        kappa=unit_kappa(),
    )

    # Budget deliberately short of ample, so the model is in its limiting-reagent regime.
    c_single = coverage(single, np.array([1.0]), B=40.0, D=40.0)
    c_pair = coverage(pair, np.array([0.5, 0.5]), B=40.0, D=40.0)

    assert c_single == pytest.approx(1.0)
    assert c_pair == pytest.approx(0.5)
    assert c_pair < c_single


def test_coverage_is_non_monotone_in_mixture_size() -> None:
    """Mean coverage over a population set rises then falls as venoms are added.

    Two target populations share 3FTx and PLA2. Adding the second immunogen helps because it
    covers PLA2. Adding four further immunogens that contribute only families neither target
    delivers must hurt, because supply is finite.
    """
    t1 = comp(**{"3FTx": 0.8, "PLA2": 0.2})
    t2 = comp(**{"3FTx": 0.2, "PLA2": 0.8})
    targets = np.vstack([t1, t2])

    pool = [
        comp(**{"3FTx": 1.0}),
        comp(**{"PLA2": 1.0}),
        comp(**{"LAAO": 1.0}),
        comp(**{"NUC": 1.0}),
        comp(**{"PDE": 1.0}),
        comp(**{"AChE": 1.0}),
    ]

    means = []
    for size in range(1, len(pool) + 1):
        immunogens = np.vstack(pool[:size])
        weights = np.full(size, 1.0 / size)
        crossreact = np.ones((2, size, N))
        c = coverage_many(
            targets,
            immunogens,
            crossreact,
            weights,
            kappa=unit_kappa(),
            B=40.0,
            D=40.0,
        )
        means.append(float(c.mean()))

    best = int(np.argmax(means))
    assert 0 < best < len(means) - 1, f"expected an interior optimum, got profile {means}"
    assert means[best] > means[0], "adding a complementary venom should help"
    assert means[-1] < means[best], "adding pure dilutants should hurt"


def test_dilution_vanishes_when_budget_is_ample() -> None:
    """The penalty is a budget effect, not an artefact: with enough antibody it disappears.

    This is the check that the dilution term is mechanistic. If coverage fell when B is effectively
    unlimited, the model would be penalising mixture size for its own sake.
    """
    target = comp(**{"3FTx": 1.0})
    a = comp(**{"3FTx": 1.0})
    b = comp(**{"SVMP": 1.0})
    pair = CoverageInputs(
        target_composition=target,
        immunogen_compositions=np.vstack([a, b]),
        crossreact=perfect_crossreact(2),
        kappa=unit_kappa(),
    )
    assert coverage(pair, np.array([0.5, 0.5]), B=1_000_000.0, D=40.0) == pytest.approx(1.0)


# --------------------------------------------------------------------------------------
# Property 4 — dose response is monotone decreasing
# --------------------------------------------------------------------------------------


def test_coverage_is_monotone_decreasing_in_dose() -> None:
    p = comp(**{"3FTx": 0.5, "SVMP": 0.3, "PLA2": 0.2})
    inputs = CoverageInputs(
        target_composition=p,
        immunogen_compositions=p[None, :],
        crossreact=perfect_crossreact(1),
        kappa=unit_kappa(),
    )
    # B = 100 with unit kappa means the budget is exhausted once D exceeds 100 mg, so the
    # grid deliberately straddles the saturation boundary.
    doses = [10.0, 50.0, 100.0, 200.0, 400.0, 800.0]
    values = [coverage(inputs, np.array([1.0]), B=100.0, D=d) for d in doses]
    assert all(later <= earlier + 1e-12 for earlier, later in pairwise(values))
    assert values[0] == pytest.approx(1.0), "a small dose is fully covered"
    assert values[-1] == pytest.approx(0.125), "at D = 8B coverage is exactly B/D"
    assert values[0] > values[-1], "a large enough dose must reduce coverage"


def test_coverage_is_monotone_increasing_in_budget() -> None:
    p = comp(**{"3FTx": 0.5, "SVMP": 0.5})
    inputs = CoverageInputs(
        target_composition=p,
        immunogen_compositions=p[None, :],
        crossreact=perfect_crossreact(1),
        kappa=unit_kappa(),
    )
    values = [coverage(inputs, np.array([1.0]), B=b, D=100.0) for b in [10.0, 50.0, 100.0, 500.0]]
    assert all(later >= earlier - 1e-12 for earlier, later in pairwise(values))


# --------------------------------------------------------------------------------------
# Property 5 — saturation: excess antibody against a family is wasted
# --------------------------------------------------------------------------------------


def test_excess_antibody_against_one_family_does_not_rescue_another() -> None:
    """Antibody is family-specific. A surplus against SVMP cannot neutralise 3FTx.

    The target is 3FTx-dominant; the immunogen is SVMP-dominant. Raising B enormously saturates
    SVMP long before it fixes the 3FTx shortfall, so coverage plateaus below 1 at the 3FTx ceiling.
    """
    target = comp(**{"3FTx": 0.9, "SVMP": 0.1})
    immunogen = comp(**{"3FTx": 0.05, "SVMP": 0.95})
    inputs = CoverageInputs(
        target_composition=target,
        immunogen_compositions=immunogen[None, :],
        crossreact=perfect_crossreact(1),
        kappa=unit_kappa(),
    )
    low = coverage(inputs, np.array([1.0]), B=40.0, D=40.0)
    high = coverage(inputs, np.array([1.0]), B=200.0, D=40.0)
    assert high > low
    n = neutralised_fractions(inputs, np.array([1.0]), B=200.0, D=40.0)
    assert n[FAMILY_INDEX["SVMP"]] == pytest.approx(1.0), "SVMP should be saturated"
    assert n[FAMILY_INDEX["3FTx"]] < 1.0, "3FTx should still be short"


def test_neutralised_fraction_never_exceeds_one() -> None:
    p = comp(**{"3FTx": 1.0})
    inputs = CoverageInputs(
        target_composition=p,
        immunogen_compositions=p[None, :],
        crossreact=perfect_crossreact(1),
        kappa=unit_kappa(),
    )
    n = neutralised_fractions(inputs, np.array([1.0]), B=1e9, D=1.0)
    assert np.all(n <= 1.0 + 1e-12)


# --------------------------------------------------------------------------------------
# Property 6 — scale invariance in B/D, and the supply budget really is finite
# --------------------------------------------------------------------------------------


def test_coverage_depends_only_on_the_ratio_of_budget_to_dose() -> None:
    p = comp(**{"3FTx": 0.4, "PLA2": 0.35, "SVMP": 0.25})
    inputs = CoverageInputs(
        target_composition=p,
        immunogen_compositions=p[None, :],
        crossreact=np.full((1, N), 0.7),
        kappa=unit_kappa(),
    )
    base = coverage(inputs, np.array([1.0]), B=100.0, D=250.0)
    for scale in [0.1, 3.0, 17.5]:
        scaled = coverage(inputs, np.array([1.0]), B=100.0 * scale, D=250.0 * scale)
        assert scaled == pytest.approx(base)


def test_supply_sums_to_the_total_budget() -> None:
    """The finite-budget invariant. If this breaks, the dilution mechanism is gone."""
    rng = np.random.default_rng(0)
    immunogens = rng.dirichlet(np.ones(N), size=4)
    weights = rng.dirichlet(np.ones(4))
    supply = antibody_supply(immunogens, weights, B=1234.5)
    assert supply.sum() == pytest.approx(1234.5)


def test_effective_supply_never_exceeds_raw_supply() -> None:
    rng = np.random.default_rng(1)
    immunogens = rng.dirichlet(np.ones(N), size=3)
    weights = rng.dirichlet(np.ones(3))
    x = rng.uniform(0.0, 1.0, size=(3, N))
    raw = antibody_supply(immunogens, weights, B=500.0)
    eff = effective_supply(immunogens, weights, x, B=500.0)
    assert np.all(eff <= raw + 1e-12)


def test_demand_is_linear_in_dose_and_kappa() -> None:
    p = comp(**{"3FTx": 0.5, "SVMP": 0.5})
    k = np.full(N, 3.0)
    d1 = antibody_demand(p, k, D=10.0)
    d2 = antibody_demand(p, k, D=20.0)
    assert np.allclose(d2, 2.0 * d1)
    assert d1[FAMILY_INDEX["3FTx"]] == pytest.approx(10.0 * 0.5 * 3.0)


# --------------------------------------------------------------------------------------
# Limiting families: the model must answer "why", not only "where"
# --------------------------------------------------------------------------------------


def test_limiting_families_names_the_worst_family_first() -> None:
    target = comp(**{"3FTx": 0.4, "SVMP": 0.4, "PLA2": 0.2})
    immunogen = comp(**{"3FTx": 0.05, "SVMP": 0.60, "PLA2": 0.35})
    inputs = CoverageInputs(
        target_composition=target,
        immunogen_compositions=immunogen[None, :],
        crossreact=perfect_crossreact(1),
        kappa=unit_kappa(),
    )
    result = evaluate("synthetic__test", inputs, np.array([1.0]), B=40.0, D=40.0)
    assert result.limiting_families[0] == "3FTx"
    assert "SVMP" not in result.limiting_families, "SVMP is over-supplied here"
    assert result.coverage + result.deficit == pytest.approx(1.0)


def test_trace_families_are_not_reported_as_limiting() -> None:
    """A family at 0.1% abundance must not be blamed for a district's deficit."""
    target = comp(**{"3FTx": 0.999, "PDE": 0.001})
    immunogen = comp(**{"3FTx": 1.0})
    inputs = CoverageInputs(
        target_composition=target,
        immunogen_compositions=immunogen[None, :],
        crossreact=perfect_crossreact(1),
        kappa=unit_kappa(),
    )
    result = evaluate("synthetic__trace", inputs, np.array([1.0]), B=100.0, D=40.0)
    assert "PDE" not in result.limiting_families


# --------------------------------------------------------------------------------------
# Severity weighting behaves as specified
# --------------------------------------------------------------------------------------


def test_severity_weighting_favours_neutralising_the_dangerous_family() -> None:
    """Two mixtures with equal total shortfall differ in coverage by severity.

    3FTx carries sigma 1.00 and NUC carries 0.10, so failing on 3FTx must cost more coverage than
    failing on NUC by the same mass.
    """
    target = comp(**{"3FTx": 0.5, "NUC": 0.5})
    kappa = unit_kappa()

    fails_on_3ftx = CoverageInputs(
        target_composition=target,
        immunogen_compositions=comp(**{"NUC": 1.0})[None, :],
        crossreact=perfect_crossreact(1),
        kappa=kappa,
    )
    fails_on_nuc = CoverageInputs(
        target_composition=target,
        immunogen_compositions=comp(**{"3FTx": 1.0})[None, :],
        crossreact=perfect_crossreact(1),
        kappa=kappa,
    )
    c_3ftx_missing = coverage(fails_on_3ftx, np.array([1.0]), B=40.0, D=40.0)
    c_nuc_missing = coverage(fails_on_nuc, np.array([1.0]), B=40.0, D=40.0)
    assert c_nuc_missing > c_3ftx_missing

    # Closed form: sigma-weighted shares are 1.0*0.5 and 0.1*0.5, normalised.
    expected_nuc_missing = (SIGMA["3FTx"] * 0.5) / (SIGMA["3FTx"] * 0.5 + SIGMA["NUC"] * 0.5)
    assert c_nuc_missing == pytest.approx(expected_nuc_missing)


# --------------------------------------------------------------------------------------
# kappa: small toxins are antibody-expensive
# --------------------------------------------------------------------------------------


def test_high_kappa_family_is_harder_to_cover() -> None:
    """kappa encodes antibody mass per unit venom mass. Doubling it halves n_f."""
    p = comp(**{"3FTx": 1.0})
    base = CoverageInputs(
        target_composition=p,
        immunogen_compositions=p[None, :],
        crossreact=perfect_crossreact(1),
        kappa=np.full(N, 1.0),
    )
    doubled = CoverageInputs(
        target_composition=p,
        immunogen_compositions=p[None, :],
        crossreact=perfect_crossreact(1),
        kappa=np.full(N, 2.0),
    )
    assert coverage(base, np.array([1.0]), B=40.0, D=40.0) == pytest.approx(1.0)
    assert coverage(doubled, np.array([1.0]), B=40.0, D=40.0) == pytest.approx(0.5)


# --------------------------------------------------------------------------------------
# coverage_many agrees with the scalar path
# --------------------------------------------------------------------------------------


def test_vectorised_and_scalar_paths_agree() -> None:
    rng = np.random.default_rng(7)
    targets = rng.dirichlet(np.ones(N) * 0.4, size=5)
    immunogens = rng.dirichlet(np.ones(N) * 0.4, size=3)
    x = rng.uniform(0.3, 1.0, size=(5, 3, N))
    weights = rng.dirichlet(np.ones(3))
    kappa = rng.uniform(0.5, 15.0, size=N)

    batch = coverage_many(targets, immunogens, x, weights, kappa, B=900.0, D=45.0)
    for t in range(5):
        one = CoverageInputs(
            target_composition=targets[t],
            immunogen_compositions=immunogens,
            crossreact=x[t],
            kappa=kappa,
        )
        assert coverage(one, weights, B=900.0, D=45.0) == pytest.approx(batch[t])


# --------------------------------------------------------------------------------------
# Input validation: bad inputs raise typed errors, they do not silently produce numbers
# --------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"kappa": np.zeros(N)}, "strictly positive"),
        ({"crossreact": np.full((1, N), 1.5)}, r"\[0, 1\]"),
    ],
)
def test_invalid_inputs_raise(kwargs: dict[str, np.ndarray], match: str) -> None:
    p = comp(**{"3FTx": 1.0})
    base = {
        "target_composition": p,
        "immunogen_compositions": p[None, :],
        "crossreact": perfect_crossreact(1),
        "kappa": np.ones(N),
    }
    base.update(kwargs)
    with pytest.raises(DataValidationError, match=match):
        CoverageInputs(**base)  # type: ignore[arg-type]


def test_weights_must_sum_to_one() -> None:
    p = comp(**{"3FTx": 1.0})
    inputs = CoverageInputs(
        target_composition=p,
        immunogen_compositions=np.vstack([p, p]),
        crossreact=perfect_crossreact(2),
        kappa=np.ones(N),
    )
    with pytest.raises(DataValidationError, match="sum to"):
        coverage(inputs, np.array([0.5, 0.9]))


def test_wrong_weight_length_raises() -> None:
    p = comp(**{"3FTx": 1.0})
    inputs = CoverageInputs(
        target_composition=p,
        immunogen_compositions=np.vstack([p, p]),
        crossreact=perfect_crossreact(2),
        kappa=np.ones(N),
    )
    with pytest.raises(DataValidationError, match="expected"):
        coverage(inputs, np.array([1.0]))
