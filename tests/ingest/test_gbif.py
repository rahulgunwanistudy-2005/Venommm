"""Occurrence handling: presence only, never abundance."""

from __future__ import annotations

import inspect

import numpy as np

from venomgap.ingest import gbif
from venomgap.ingest.gbif import Occurrence, range_mask


def occ(lat: float, lon: float) -> Occurrence:
    return Occurrence(lat=lat, lon=lon, state="", year=None, basis="HUMAN_OBSERVATION")


def test_range_mask_includes_nearby_and_excludes_distant() -> None:
    lats = np.array([13.0, 28.6])
    lons = np.array([80.2, 77.2])
    mask = range_mask(lats, lons, [occ(13.1, 80.3)], radius_km=250.0)
    assert bool(mask[0]) is True
    assert bool(mask[1]) is False


def test_no_occurrences_means_no_range() -> None:
    mask = range_mask(np.array([13.0]), np.array([80.2]), [], radius_km=250.0)
    assert not mask.any()


def test_mask_is_boolean_and_carries_no_count() -> None:
    """The central guarantee: occurrence *counts* must not survive into the model.

    A district with one nearby record and a district with a thousand must be indistinguishable,
    because GBIF is effort-biased and a count would import that bias as abundance.
    """
    lats, lons = np.array([13.0]), np.array([80.2])
    one = range_mask(lats, lons, [occ(13.05, 80.25)], radius_km=250.0)
    many = range_mask(lats, lons, [occ(13.0 + i * 0.001, 80.2) for i in range(500)],
                      radius_km=250.0)
    assert one.dtype == bool
    assert np.array_equal(one, many)


def test_range_mask_returns_one_entry_per_district() -> None:
    lats = np.array([13.0, 20.0, 28.6])
    lons = np.array([80.2, 77.0, 77.2])
    mask = range_mask(lats, lons, [occ(13.1, 80.3)], radius_km=250.0)
    assert mask.shape == lats.shape


def test_the_module_documents_that_counts_are_not_abundance() -> None:
    """Belt and braces: the prohibition is stated where the next person will read it."""
    doc = inspect.getdoc(gbif) or ""
    assert "not abundance" in doc.lower()


def test_occurrence_summary_says_what_counts_are_not_for() -> None:
    summary = gbif.occurrence_summary({"Naja naja": [occ(13.0, 80.2)]})
    assert "never used as abundance" in str(summary["note"])
    assert summary["records"] == {"Naja naja": 1}
