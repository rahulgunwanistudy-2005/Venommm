"""Snakebite burden: state death rates apportioned to districts.

The weighting that turns per-district coverage into a national objective comes from Suraweera et
al. 2020 (eLife 9:e54076), the Million Death Study analysis of snakebite mortality, Table 3.

The apportionment is deliberately crude and deliberately explicit: a state's estimated deaths are
divided among its districts in proportion to **rural** population, because snakebite mortality in
India is overwhelmingly rural. Nothing here claims district-level mortality is known. It is not.
The burden weight is a stated allocation assumption, and the sensitivity analysis re-runs the
optimiser under a uniform-within-state allocation so that the effect of the assumption is measured
rather than asserted.
"""

from __future__ import annotations

import logging

import pandas as pd

from venomgap.config import BURDEN_CSV, DISTRICTS_CSV
from venomgap.errors import DataValidationError

logger = logging.getLogger(__name__)

# States whose Suraweera row is an aggregate rather than a single state. Districts in these states
# inherit the aggregate rate, and the fact is carried through to the output.
AGGREGATE_ROWS: frozenset[str] = frozenset({"Northeastern states", "All other states"})


def load_state_burden(path=BURDEN_CSV) -> pd.DataFrame:  # type: ignore[no-untyped-def]
    """State-level snakebite mortality, validated."""
    if not path.exists():
        raise DataValidationError(f"state burden file not found: {path}")
    frame = pd.read_csv(path)
    required = {"state", "death_rate_2010_2014", "estimated_deaths_thousands_2001_14", "source_doi"}
    missing = required - set(frame.columns)
    if missing:
        raise DataValidationError(f"{path} is missing columns: {sorted(missing)}")
    if frame["death_rate_2010_2014"].isna().any():
        raise DataValidationError(f"{path} has a missing death rate")
    if (frame["death_rate_2010_2014"] < 0).any():
        raise DataValidationError(f"{path} has a negative death rate")
    return frame


def load_districts(path=DISTRICTS_CSV) -> pd.DataFrame:  # type: ignore[no-untyped-def]
    """District table with centroid, state and rural population."""
    if not path.exists():
        raise DataValidationError(f"district file not found: {path}")
    frame = pd.read_csv(path)
    required = {"district_id", "district", "state", "lat", "lon", "rural_population"}
    missing = required - set(frame.columns)
    if missing:
        raise DataValidationError(f"{path} is missing columns: {sorted(missing)}")
    if frame["district_id"].duplicated().any():
        duplicated = frame.loc[frame["district_id"].duplicated(), "district_id"].tolist()
        raise DataValidationError(f"{path} has duplicate district_id values: {duplicated[:5]}")
    return frame


def district_burden(
    districts: pd.DataFrame,
    state_burden: pd.DataFrame,
    uniform_within_state: bool = False,
) -> pd.DataFrame:
    """Attach a burden weight to every district.

    `burden_weight` is the district's share of estimated national annual snakebite deaths, so the
    weights sum to 1 across India and the national objective is a true burden-weighted mean.

    With `uniform_within_state=True` the within-state split is by district count instead of rural
    population. That is the sensitivity variant, not an alternative truth.
    """
    rates = state_burden.set_index("state")
    merged = districts.copy()

    unmatched = sorted(set(merged["state"]) - set(rates.index))
    if unmatched:
        raise DataValidationError(
            f"districts reference states with no burden row: {unmatched}. Every state must map to "
            f"a Suraweera row or to an explicit aggregate."
        )

    merged["state_deaths"] = merged["state"].map(
        rates["estimated_deaths_thousands_2001_14"].to_dict()
    )
    merged["state_death_rate"] = merged["state"].map(rates["death_rate_2010_2014"].to_dict())
    merged["burden_is_aggregate"] = merged["state"].map(
        rates.index.to_series().isin(AGGREGATE_ROWS).to_dict()
    )

    if uniform_within_state:
        share = merged.groupby("state")["district_id"].transform("count")
        merged["within_state_share"] = 1.0 / share
    else:
        totals = merged.groupby("state")["rural_population"].transform("sum")
        # A state whose rural population is entirely missing falls back to a uniform split rather
        # than producing NaN weights; the fallback is logged, not hidden.
        safe = totals.replace(0.0, pd.NA)
        merged["within_state_share"] = merged["rural_population"] / safe
        degenerate = merged["within_state_share"].isna()
        if degenerate.any():
            counts = merged.groupby("state")["district_id"].transform("count")
            merged.loc[degenerate, "within_state_share"] = 1.0 / counts[degenerate]
            logger.warning(
                "%d districts had no rural population and fell back to a uniform within-state "
                "share", int(degenerate.sum())
            )

    raw = merged["state_deaths"] * merged["within_state_share"]
    total = float(raw.sum())
    if total <= 0.0:
        raise DataValidationError("total burden is zero; cannot normalise")
    merged["burden_weight"] = raw / total

    logger.info(
        "burden attached to %d districts across %d states; top 5 states carry %.1f%% of weight",
        len(merged),
        merged["state"].nunique(),
        100.0
        * merged.groupby("state")["burden_weight"].sum().nlargest(5).sum(),
    )
    return merged
