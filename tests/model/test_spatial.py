"""M6 gate: kernel interpolation, the fitted length scale, and honest `unknown` labelling."""

from __future__ import annotations

from itertools import pairwise

import numpy as np
import pytest

from venomgap.config import FAMILIES
from venomgap.errors import SpatialModelError
from venomgap.ingest.compositions import load_compositions
from venomgap.model.spatial import (
    SpatialCompositionModel,
    fit_length_scale,
    gaussian_kernel,
    haversine_km,
)

N = len(FAMILIES)


@pytest.fixture(scope="module")
def populations() -> list:
    return load_compositions()


# --------------------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------------------


def test_haversine_against_known_distances() -> None:
    # Chennai to Delhi is about 1760 km great-circle.
    distance = float(haversine_km(13.08, 80.27, 28.61, 77.21))
    assert 1700 < distance < 1820
    assert float(haversine_km(0.0, 0.0, 0.0, 0.0)) == pytest.approx(0.0)


def test_haversine_is_symmetric() -> None:
    a = float(haversine_km(12.0, 78.0, 26.0, 73.0))
    b = float(haversine_km(26.0, 73.0, 12.0, 78.0))
    assert a == pytest.approx(b)


def test_kernel_decays_and_is_one_at_zero() -> None:
    assert float(gaussian_kernel(np.array([0.0]), 500.0)[0]) == pytest.approx(1.0)
    values = gaussian_kernel(np.array([0.0, 250.0, 500.0, 1000.0]), 500.0)
    assert all(later < earlier for earlier, later in pairwise(values))


def test_non_positive_length_scale_raises() -> None:
    with pytest.raises(SpatialModelError, match="positive"):
        gaussian_kernel(np.array([1.0]), 0.0)


# --------------------------------------------------------------------------------------
# Interpolation
# --------------------------------------------------------------------------------------


def test_estimate_at_a_sampled_point_is_close_to_that_population(populations: list) -> None:
    model = SpatialCompositionModel(populations, ell_km=50.0)
    target = next(p for p in populations if p.species == "Naja naja")
    estimate = model.estimate(target.lat, target.lon, target.species)
    # With a tight kernel the nearest population dominates.
    assert float(np.linalg.norm(estimate.composition - target.vector())) < 0.5
    assert estimate.nearest_km == pytest.approx(0.0, abs=1.0)


def test_estimates_always_lie_on_the_simplex(populations: list) -> None:
    model = SpatialCompositionModel(populations, ell_km=800.0)
    for lat, lon in [(28.6, 77.2), (11.6, 92.7), (19.0, 73.0), (26.9, 71.9)]:
        for species in {p.species for p in populations}:
            estimate = model.estimate(lat, lon, species)
            assert estimate.composition.sum() == pytest.approx(1.0, abs=1e-9)
            assert np.all(estimate.composition >= -1e-12)


def test_distant_points_are_reported_unknown(populations: list) -> None:
    """A point far from every proteome must be unknown, not confidently covered."""
    model = SpatialCompositionModel(populations, ell_km=100.0)
    estimate = model.estimate(-40.0, 10.0, "Naja naja")  # the South Atlantic
    assert estimate.status == "unknown"
    assert not estimate.is_known


def test_uncertainty_rises_with_distance(populations: list) -> None:
    model = SpatialCompositionModel(populations, ell_km=600.0)
    target = next(p for p in populations if p.species == "Daboia russelii")
    near = model.estimate(target.lat, target.lon, target.species)
    far = model.estimate(target.lat + 12.0, target.lon + 12.0, target.species)
    assert far.uncertainty >= near.uncertainty


def test_leave_one_out_exclusion_actually_excludes(populations: list) -> None:
    model = SpatialCompositionModel(populations, ell_km=300.0)
    target = next(p for p in populations if p.species == "Naja naja")
    with_self = model.estimate(target.lat, target.lon, target.species)
    without_self = model.estimate(target.lat, target.lon, target.species, exclude=target.pop_id)
    assert target.pop_id in with_self.contributors
    assert target.pop_id not in without_self.contributors


def test_subspecies_borrows_from_its_nominate_species(populations: list) -> None:
    """`Echis carinatus sochureki` must be able to use `Echis carinatus` proteomes as neighbours."""
    model = SpatialCompositionModel(populations, ell_km=600.0)
    neighbours = {p.pop_id for p in model.neighbours("Echis carinatus sochureki")}
    assert any(pop_id.startswith("Echis_carinatus__") for pop_id in neighbours)


def test_species_with_no_proteome_falls_back_to_its_genus(populations: list) -> None:
    model = SpatialCompositionModel(populations, ell_km=600.0)
    neighbours = model.neighbours("Bungarus sindanus")
    assert neighbours, "a species with sparse data must still get genus-level neighbours"
    assert all(p.genus == "Bungarus" for p in neighbours)


# --------------------------------------------------------------------------------------
# The fit
# --------------------------------------------------------------------------------------


def test_length_scale_fit_is_deterministic_and_reports_its_plateau(populations: list) -> None:
    first = fit_length_scale(populations)
    second = fit_length_scale(populations)
    assert first.ell_km == second.ell_km
    assert first.rmse == pytest.approx(second.rmse)
    low, high = first.plateau_km
    assert low <= first.ell_km <= high
    assert first.folds > 0


def test_the_fit_uses_only_composition_information(populations: list) -> None:
    """The length scale must not depend on any antivenom quantity, which is what makes it safe to
    fit on the whole corpus including holdout populations."""
    holdout = [p for p in populations if p.holdout]
    assert holdout, "the corpus should contain holdout populations"
    fit = fit_length_scale(populations)
    # Every holdout population takes part in the CV; if it did not, the fit would be using
    # holdout status, which is antivenom information.
    assert any(p.pop_id in fit.per_population_error for p in holdout)


def test_too_few_populations_raises() -> None:
    with pytest.raises(SpatialModelError, match="at least three"):
        fit_length_scale(load_compositions()[:2])
