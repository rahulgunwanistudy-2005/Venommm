"""GBIF occurrence ingest: which snakes are plausibly where.

**Occurrence counts are not abundance.** GBIF records are presence-only and heavily effort-biased:
a district near a university has more records than a district that is harder to reach, and that
says nothing about how many snakes live there. So occurrences are used in exactly one direction --
to *exclude* a district from a species' range when the species is not plausibly present -- and
never to infer how common a species is. A test asserts that no abundance is derived from counts.

Range membership is decided by proximity: a district is in range if a georeferenced occurrence of
that species falls within `RANGE_RADIUS_KM` of its centroid. Species whose records are too sparse
to support that fall back to a published range description held in `data/curated/species_range.csv`
and are flagged.
"""

from __future__ import annotations

import json
import logging
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from urllib.parse import urlencode

import numpy as np
import requests
from numpy.typing import NDArray

from venomgap.config import OCCURRENCE_CACHE, OFFLINE
from venomgap.errors import MissingCacheError, VenomGapError

logger = logging.getLogger(__name__)

GBIF_SEARCH = "https://api.gbif.org/v1/occurrence/search"
REQUEST_TIMEOUT_S = 60
MAX_RETRIES = 3
RETRY_BACKOFF_S = 3.0
PAGE_LIMIT = 300
MAX_RECORDS = 3000

# A district centroid within this distance of a georeferenced occurrence is treated as in range.
# Generous on purpose: the cost of wrongly excluding a district is a false "unknown", the cost of
# wrongly including one is a diluted species mixture, and the first is the more honest failure.
RANGE_RADIUS_KM = 250.0

# Below this many georeferenced records a species' GBIF range is not trustworthy on its own.
MIN_RECORDS_FOR_RANGE = 25


@dataclass(frozen=True, slots=True)
class Occurrence:
    lat: float
    lon: float
    state: str
    year: int | None
    basis: str


def _cache_path(species: str) -> Path:
    return OCCURRENCE_CACHE / f"{species.replace(' ', '_')}.json"


def fetch_species_occurrences(
    species: str, refresh: bool = False, country: str = "IN"
) -> list[Occurrence]:
    """Georeferenced occurrences for a species in one country, cached to disk."""
    path = _cache_path(species)
    if path.exists() and not refresh:
        payload = json.loads(path.read_text())
        return [Occurrence(**record) for record in payload["occurrences"]]

    if OFFLINE:
        raise MissingCacheError(
            f"offline mode is set and {path} is not cached; "
            f"run `python -m venomgap.cli fetch-occurrences` with network access first"
        )

    records: list[Occurrence] = []
    offset = 0
    while offset < MAX_RECORDS:
        params = {
            "scientificName": species,
            "country": country,
            "hasCoordinate": "true",
            "hasGeospatialIssue": "false",
            "limit": str(PAGE_LIMIT),
            "offset": str(offset),
        }
        payload = _get_json(f"{GBIF_SEARCH}?{urlencode(params)}")
        results = payload.get("results", [])
        for row in results:
            lat, lon = row.get("decimalLatitude"), row.get("decimalLongitude")
            if lat is None or lon is None:
                continue
            records.append(
                Occurrence(
                    lat=float(lat),
                    lon=float(lon),
                    state=str(row.get("stateProvince") or ""),
                    year=int(row["year"]) if row.get("year") else None,
                    basis=str(row.get("basisOfRecord") or ""),
                )
            )
        if payload.get("endOfRecords", True) or not results:
            break
        offset += PAGE_LIMIT

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {
                "species": species,
                "country": country,
                "query": "presence-only occurrences; counts are NOT used as abundance",
                "records": len(records),
                "occurrences": [asdict(r) for r in records],
            },
            indent=1,
        )
    )
    logger.info("%s: %d georeferenced occurrences in %s", species, len(records), country)
    return records


def _get_json(url: str) -> dict[str, object]:
    last_error: Exception | None = None
    for attempt in range(MAX_RETRIES):
        try:
            response = requests.get(url, timeout=REQUEST_TIMEOUT_S)
            response.raise_for_status()
            data: dict[str, object] = response.json()
        except (requests.RequestException, ValueError) as exc:
            last_error = exc
            logger.warning("GBIF request failed (%d/%d): %s", attempt + 1, MAX_RETRIES, exc)
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_BACKOFF_S * (attempt + 1))
        else:
            return data
    raise VenomGapError(f"GBIF request failed after {MAX_RETRIES} attempts: {last_error}")


def range_mask(
    district_lats: NDArray[np.float64],
    district_lons: NDArray[np.float64],
    occurrences: list[Occurrence],
    radius_km: float = RANGE_RADIUS_KM,
) -> NDArray[np.bool_]:
    """True where a district centroid is within `radius_km` of any occurrence.

    This is a presence mask and nothing more. The number of nearby occurrences is deliberately
    discarded so that it cannot leak into a weighting anywhere downstream.
    """
    from venomgap.model.spatial import haversine_km

    if not occurrences:
        return np.zeros(district_lats.shape, dtype=bool)
    occ_lats = np.array([o.lat for o in occurrences])
    occ_lons = np.array([o.lon for o in occurrences])
    mask = np.zeros(district_lats.shape, dtype=bool)
    # Chunked so a large occurrence set does not build a huge distance matrix at once.
    for start in range(0, len(occ_lats), 500):
        chunk_lat = occ_lats[start : start + 500]
        chunk_lon = occ_lons[start : start + 500]
        distances = haversine_km(
            district_lats[:, None], district_lons[:, None], chunk_lat[None, :], chunk_lon[None, :]
        )
        mask |= np.any(distances <= radius_km, axis=1)
    return mask


def occurrence_summary(species_occurrences: dict[str, list[Occurrence]]) -> dict[str, object]:
    """Record counts for SOURCES.md, with an explicit note about what they are not used for."""
    return {
        "note": (
            "Record counts are reported for provenance only. They are presence-only and "
            "effort-biased, and are never used as abundance or to weight a species."
        ),
        "records": {
            species: len(records) for species, records in sorted(species_occurrences.items())
        },
        "sparse_species": sorted(
            species
            for species, records in species_occurrences.items()
            if len(records) < MIN_RECORDS_FOR_RANGE
        ),
    }
