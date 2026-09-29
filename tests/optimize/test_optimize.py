"""M7 gate: the solvers, and the properties that make the dilution effect real rather than a bug."""

from __future__ import annotations

import numpy as np
import pytest

from venomgap.errors import OptimisationError
from venomgap.optimize.exact import optimality_gap


def test_optimality_gap_is_zero_when_the_heuristic_matches() -> None:
    assert optimality_gap(0.65, 0.65) == pytest.approx(0.0)


def test_optimality_gap_is_relative_and_non_negative() -> None:
    assert optimality_gap(0.60, 0.80) == pytest.approx(0.25)
    # A heuristic cannot beat the exact optimum on the same instance; if it appears to,
    # report zero rather than a negative gap.
    assert optimality_gap(0.90, 0.80) == pytest.approx(0.0)


def test_optimality_gap_is_nan_for_a_degenerate_instance() -> None:
    assert np.isnan(optimality_gap(0.0, 0.0))


@pytest.mark.slow
def test_greedy_stops_rather_than_accepting_a_worse_portfolio() -> None:
    """The early stop is the dilution effect inside the solver, and it must not be silently
    overridden: a larger budget can never return worse coverage."""
    from venomgap.optimize.greedy import candidate_populations, greedy_sites
    from venomgap.pipeline import build_context

    context = build_context()
    base = context.assembly.big_four_immunogen()
    candidates = candidate_populations(context.assembly, indian_only=True)
    previous = context.objective.national_coverage(base)
    for k in range(1, 5):
        result = greedy_sites(context.objective, context.assembly, k, candidates, base=base)
        assert result.coverage >= previous - 1e-9, (
            f"greedy returned worse coverage at k={k}: {result.coverage} < {previous}"
        )
        assert len(result.sites) <= k
        previous = result.coverage


@pytest.mark.slow
def test_weight_optimisation_never_returns_worse_than_uniform() -> None:
    from venomgap.optimize.weights import optimise_weights
    from venomgap.pipeline import build_context

    context = build_context()
    immunogen = context.assembly.big_four_immunogen()
    result = optimise_weights(context.objective, immunogen)
    assert result.coverage >= result.uniform_coverage - 1e-12
    assert result.weights.sum() == pytest.approx(1.0)
    assert np.all(result.weights >= -1e-12)


@pytest.mark.slow
def test_exact_enumeration_refuses_an_instance_it_cannot_enumerate() -> None:
    from venomgap.optimize.exact import exact_best
    from venomgap.pipeline import build_context

    context = build_context()
    with pytest.raises(OptimisationError, match="limited to"):
        exact_best(context.objective, context.assembly, k=9)


@pytest.mark.slow
def test_combined_immunogen_extends_rather_than_replaces() -> None:
    """k new sites means the Big Four plus k, not k on their own."""
    from venomgap.pipeline import build_context

    context = build_context()
    base = context.assembly.big_four_immunogen()
    added = context.assembly.combined_immunogen(base, ("Naja_naja__Rajasthan",))
    assert added.size == base.size + 1
    assert added.labels[: base.size] == base.labels
    assert added.weights.sum() == pytest.approx(1.0)
    # An empty addition is a no-op.
    assert context.assembly.combined_immunogen(base, ()).size == base.size
