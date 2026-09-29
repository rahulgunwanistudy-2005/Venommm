"""Parameter fitting, and the guard that makes the blinded retrodiction meaningful.

Two things are fitted and nothing else:

* `kappa_scale` -- one global number. The *shape* of `kappa_f` is fixed a priori by molar
  equivalence (`kappa_f = kappa_scale * M_ab / M_f`), so the fit can rescale how antibody-expensive
  venom is in total but cannot decide that one family is cheaper than the mass argument says.
* `theta_f` -- the sequence-identity recognition floor, for the six families the calibration set has
  evidence about. Every other family keeps the a priori floor exactly.

Plus one nuisance parameter, `C_star`: the coverage level the model should predict at a murine ED50.
It is not a claim about clinical neutralisation, only the level the fit anchors on.

**The guard.** `fit_parameters` raises `HoldoutLeakError` if any input row is labelled
`holdout=True` **or** cites one of the four holdout studies. That second check matters: the 2019
"beyond the big four" paper supplies compositions for several non-holdout populations, and its
compositions are legitimate model inputs, but none of its antivenom results may enter a fit.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from scipy.optimize import minimize

from venomgap.config import (
    CALIBRATION_CSV,
    FAMILIES,
    FITTED_PARAMS_JSON,
    HOLDOUT_STUDY_DOIS,
    KAPPA_L2_PENALTY,
    KAPPA_PRIOR,
    KAPPA_SCALE_BOUNDS,
    RANDOM_SEED,
    THETA_BOUNDS,
    THETA_IDENTIFIABILITY_TOL,
    THETA_L2_PENALTY,
    THETA_PRIOR,
    THETA_PROBE_DELTA,
)
from venomgap.errors import CalibrationError, DataValidationError, HoldoutLeakError
from venomgap.model.assemble import ModelAssembly, theta_from_parameters
from venomgap.model.coverage import coverage, neutralised_fractions
from venomgap.types import FittedParameters

logger = logging.getLogger(__name__)

# Families the calibration set carries direct evidence about. Everything else keeps THETA_PRIOR.
FITTED_THETA_FAMILIES: tuple[str, ...] = ("3FTx", "PLA2", "SVMP", "SVSP", "CTL", "KSPI")

# Minimum gap the ordinal constraint asks for between a poorly and a well recognised family.
ORDINAL_MARGIN = 0.10

# Relative weights of the three evidence terms. Potency data is quantitative so it dominates; the
# ordinal term is qualitative and is weighted to inform without driving the fit.
WEIGHT_POTENCY = 1.0
WEIGHT_CENSORED = 1.0
WEIGHT_ORDINAL = 0.25

C_STAR_BOUNDS = (0.20, 0.999)


@dataclass(frozen=True, slots=True)
class CalibrationSet:
    """The validated, guard-checked calibration data, split by evidence type."""

    potency: pd.DataFrame
    censored: pd.DataFrame
    poorly_recognised: pd.DataFrame
    well_recognised: pd.DataFrame
    studies: tuple[str, ...]

    @property
    def rows(self) -> int:
        return (
            len(self.potency)
            + len(self.censored)
            + len(self.poorly_recognised)
            + len(self.well_recognised)
        )


def load_calibration(path=CALIBRATION_CSV, exclude_study: str | None = None) -> CalibrationSet:
    """Load the calibration CSV and enforce the holdout split. Raises rather than filtering."""
    if not path.exists():
        raise DataValidationError(f"calibration file not found: {path}")
    frame = pd.read_csv(path)

    required = {
        "record_id", "study", "doi", "table", "antivenom", "pop_id",
        "observable", "family", "observed_value", "censored", "holdout",
    }
    missing = required - set(frame.columns)
    if missing:
        raise DataValidationError(f"{path} is missing columns: {sorted(missing)}")

    _assert_no_holdout(frame, str(path))

    if exclude_study is not None:
        before = len(frame)
        frame = frame[frame["study"] != exclude_study]
        logger.info(
            "leave-one-study-out: excluded %r, %d of %d rows remain",
            exclude_study, len(frame), before,
        )
        if frame.empty:
            raise CalibrationError(f"excluding study {exclude_study!r} leaves no calibration data")

    potency = frame[
        (frame["observable"] == "neutralisation_potency_mg_per_ml") & (frame["censored"] == "none")
    ].copy()
    censored = frame[
        (frame["observable"] == "neutralisation_potency_mg_per_ml") & (frame["censored"] != "none")
    ].copy()
    poorly = frame[frame["observable"] == "poorly_recognised_family"].copy()
    well = frame[frame["observable"] == "well_recognised_family"].copy()

    recognised = len(potency) + len(censored) + len(poorly) + len(well)
    if recognised != len(frame):
        unknown = sorted(set(frame["observable"]) - {
            "neutralisation_potency_mg_per_ml",
            "poorly_recognised_family",
            "well_recognised_family",
        })
        raise DataValidationError(f"{path} has unhandled observables: {unknown}")

    return CalibrationSet(
        potency=potency,
        censored=censored,
        poorly_recognised=poorly,
        well_recognised=well,
        studies=tuple(sorted(frame["study"].unique())),
    )


def _assert_no_holdout(frame: pd.DataFrame, source: str) -> None:
    """The split is enforced here, in code. Not a convention -- an assertion."""
    flagged = frame[frame["holdout"].astype(str).str.lower().isin({"true", "1", "yes"})]
    if not flagged.empty:
        raise HoldoutLeakError(
            f"{source} contains {len(flagged)} row(s) labelled holdout=True: "
            f"{sorted(flagged['record_id'])}. Calibration must never see a holdout row."
        )
    leaked = frame[frame["doi"].isin(HOLDOUT_STUDY_DOIS)]
    if not leaked.empty:
        raise HoldoutLeakError(
            f"{source} cites holdout studies {sorted(leaked['doi'].unique())} in rows "
            f"{sorted(leaked['record_id'])}. A holdout study's compositions are legitimate model "
            f"inputs, but none of its antivenom results may enter a parameter fit."
        )


def fit_parameters(
    assembly: ModelAssembly,
    calibration: CalibrationSet,
    seed: int = RANDOM_SEED,
) -> FittedParameters:
    """Fit `kappa_scale` and the six supported `theta_f` against the calibration split.

    The potency term is the heart of it. At a murine ED50 every population has, by definition,
    reached the same level of protection, so the model evaluated at each population's own observed
    antivenom-to-venom ratio should return the same coverage. That ratio is
    `B/D = protein_concentration / potency` in mg antibody per mg venom, using the protein
    concentration the source study measured for the batch it tested. What remains unknown -- the
    venom-specific fraction of that protein -- is what `kappa_scale` absorbs.
    """
    immunogen = assembly.big_four_immunogen()
    rng = np.random.default_rng(seed)

    n_theta = len(FITTED_THETA_FAMILIES)
    theta_low, theta_high = THETA_BOUNDS

    def unpack(x: NDArray[np.float64]) -> tuple[float, NDArray[np.float64], float]:
        kappa_scale = float(np.exp(x[0]))
        theta_global = float(x[1])
        deviations = dict(zip(FITTED_THETA_FAMILIES, x[2 : 2 + n_theta], strict=True))
        theta = np.clip(
            theta_from_parameters(theta_global, deviations), theta_low, theta_high
        )
        c_star = float(x[2 + n_theta])
        return kappa_scale, theta, c_star

    def data_loss(x: NDArray[np.float64]) -> float:
        """The evidence terms only. Probing this is what reveals whether a parameter is pinned by
        data or merely by its prior."""
        kappa_scale, theta, c_star = unpack(x)
        total = 0.0

        # 1. Potency consistency.
        for _, row in calibration.potency.iterrows():
            ratio = _ratio(row)
            if ratio is None:
                continue
            inputs, _ = assembly.inputs_for_population(
                str(row["pop_id"]), immunogen, theta, kappa_scale
            )
            predicted = coverage(inputs, immunogen.weights, B=ratio, D=1.0)
            total += WEIGHT_POTENCY * (predicted - c_star) ** 2

        # 2. Censored failures: coverage must stay below the anchor even at a generous ratio.
        for _, row in calibration.censored.iterrows():
            ratio = _ratio(row)
            if ratio is None:
                continue
            inputs, _ = assembly.inputs_for_population(
                str(row["pop_id"]), immunogen, theta, kappa_scale
            )
            predicted = coverage(inputs, immunogen.weights, B=ratio, D=1.0)
            total += WEIGHT_CENSORED * max(0.0, predicted - c_star) ** 2

        # 3. Ordinal per-family constraints from published antivenomics.
        for pop_id, poor_families, good_families in _ordinal_groups(calibration):
            inputs, _ = assembly.inputs_for_population(
                pop_id, immunogen, theta, kappa_scale
            )
            n = neutralised_fractions(inputs, immunogen.weights, B=1.0, D=1.0)
            for poor in poor_families:
                for good in good_families:
                    gap = n[FAMILIES.index(poor)] - n[FAMILIES.index(good)] + ORDINAL_MARGIN
                    total += WEIGHT_ORDINAL * max(0.0, float(gap)) ** 2

        return total

    def penalty(x: NDArray[np.float64]) -> float:
        """Regularisation toward the a priori values."""
        return (
            KAPPA_L2_PENALTY * x[0] ** 2
            + THETA_L2_PENALTY * (float(x[1]) - THETA_PRIOR) ** 2
            + THETA_L2_PENALTY * float(np.sum(np.square(x[2 : 2 + n_theta])))
        )

    def objective(x: NDArray[np.float64]) -> float:
        return data_loss(x) + penalty(x)

    bounds = [
        (np.log(KAPPA_SCALE_BOUNDS[0]), np.log(KAPPA_SCALE_BOUNDS[1])),
        THETA_BOUNDS,
        *[(-0.20, 0.20)] * n_theta,
        C_STAR_BOUNDS,
    ]

    best: tuple[float, NDArray[np.float64]] | None = None
    # Multi-start: the objective has hinge terms and is not convex, so a single local solve from one
    # point is not enough. The seed is fixed, so the restarts are reproducible.
    starts = [np.array([0.0, THETA_PRIOR, *([0.0] * n_theta), 0.5])]
    for _ in range(7):
        starts.append(
            np.array([
                rng.uniform(-1.5, 1.5),
                rng.uniform(theta_low, theta_high),
                *rng.uniform(-0.15, 0.15, size=n_theta),
                rng.uniform(*C_STAR_BOUNDS),
            ])
        )
    for start in starts:
        result = minimize(objective, start, method="L-BFGS-B", bounds=bounds)
        if best is None or result.fun < best[0]:
            best = (float(result.fun), np.asarray(result.x, dtype=np.float64))

    if best is None:
        raise CalibrationError("no calibration start converged")

    loss, x = best
    kappa_scale, theta, c_star = unpack(x)

    residuals = []
    for _, row in calibration.potency.iterrows():
        ratio = _ratio(row)
        if ratio is None:
            continue
        inputs, _ = assembly.inputs_for_population(
            str(row["pop_id"]), immunogen, theta, kappa_scale
        )
        predicted = coverage(inputs, immunogen.weights, B=ratio, D=1.0)
        residuals.append(predicted - c_star)
    rmse = float(np.sqrt(np.mean(np.square(residuals)))) if residuals else float("nan")

    at_bound = [
        family
        for family, value in zip(FAMILIES, theta, strict=True)
        if abs(value - theta_low) < 1e-6 or abs(value - theta_high) < 1e-6
    ]
    if len(at_bound) > len(FAMILIES) // 2:
        # Pre-registered falsification signal, section 7. Reported, not silently accepted.
        logger.warning(
            "theta is at a bound for %d of %d families (%s); the pre-registration lists this as "
            "an indication that the identity signal carries little information",
            len(at_bound), len(FAMILIES), at_bound,
        )

    # How much does the calibration objective actually care about each fitted theta? A family whose
    # objective barely moves is not constrained by this data, and reporting it as "fitted" would
    # overstate what the calibration did. The dominant reason here is structural: the calibration
    # targets are populations of species that are themselves in the immunogen set, where
    # x_f(q -> p) = 1 by construction and theta has no purchase.
    identifiability: dict[str, float] = {}
    base_data_loss = data_loss(x)
    for index, family in enumerate(FITTED_THETA_FAMILIES):
        probe = x.copy()
        probe[2 + index] = float(np.clip(x[2 + index] + THETA_PROBE_DELTA, -0.20, 0.20))
        up = data_loss(probe)
        probe[2 + index] = float(np.clip(x[2 + index] - THETA_PROBE_DELTA, -0.20, 0.20))
        down = data_loss(probe)
        identifiability[family] = float(
            max(abs(up - base_data_loss), abs(down - base_data_loss))
        )
    unidentified = tuple(
        family
        for family, sensitivity in identifiability.items()
        if sensitivity < THETA_IDENTIFIABILITY_TOL
    )
    if unidentified:
        logger.warning(
            "theta_f is not identified by the calibration data for %s; these keep their a priori "
            "value of %.2f. The calibration targets share species with the immunogen set, where "
            "cross-reactivity is 1 by construction.",
            list(unidentified), THETA_PRIOR,
        )

    # Internal consistency check: does the a priori kappa shape reproduce the published ordering of
    # well- against poorly-recognised families? This is not fitted, so a violation is information.
    violations: list[str] = []
    checks = 0
    for pop_id, poor_families, good_families in _ordinal_groups(calibration):
        inputs, _ = assembly.inputs_for_population(pop_id, immunogen, theta, kappa_scale)
        n = neutralised_fractions(inputs, immunogen.weights, B=1.0, D=1.0)
        for poor in poor_families:
            for good in good_families:
                checks += 1
                if n[FAMILIES.index(poor)] > n[FAMILIES.index(good)]:
                    violations.append(f"{pop_id}:{poor}>{good}")

    kappa = {f: kappa_scale * KAPPA_PRIOR[f] for f in FAMILIES}
    fitted = FittedParameters(
        c_star=c_star,
        theta_identifiability=identifiability,
        theta_unidentified=unidentified,
        ordinal_violations=tuple(violations),
        ordinal_checks=checks,
        theta={f: float(v) for f, v in zip(FAMILIES, theta, strict=True)},
        kappa=kappa,
        kappa_scale=kappa_scale,
        calibration_rows=calibration.rows,
        calibration_rmse=rmse,
        calibration_studies=calibration.studies,
        excluded_study=None,
        converged=bool(np.isfinite(loss)),
    )
    logger.info(
        "calibrated on %d rows from %d studies: kappa_scale=%.4f, C*=%.3f, RMSE=%.4f, loss=%.5f",
        calibration.rows, len(calibration.studies), kappa_scale, c_star, rmse, loss,
    )
    return fitted


def _ratio(row: pd.Series[object]) -> float | None:
    """Antibody-to-venom mass ratio implied by a measured neutralisation potency.

    potency is mg of venom neutralised per mL of antivenom, and the study also measured the
    antivenom's protein concentration in mg/mL, so their quotient is mg antibody protein per mg
    venom -- exactly the `B/D` the coverage model takes.
    """
    potency = float(row["observed_value"])
    if potency <= 0.0:
        return None
    concentration = float(row["antivenom_protein_mg_per_ml"])
    return concentration / potency


def _ordinal_groups(
    calibration: CalibrationSet,
) -> list[tuple[str, list[str], list[str]]]:
    """Group the ordinal rows into (pop_id, poorly recognised, well recognised) triples."""
    groups: list[tuple[str, list[str], list[str]]] = []
    populations = set(calibration.poorly_recognised["pop_id"]) & set(
        calibration.well_recognised["pop_id"]
    )
    for pop_id in sorted(populations):
        poor = sorted(
            calibration.poorly_recognised.loc[
                calibration.poorly_recognised["pop_id"] == pop_id, "family"
            ]
        )
        good = sorted(
            calibration.well_recognised.loc[
                calibration.well_recognised["pop_id"] == pop_id, "family"
            ]
        )
        if poor and good:
            groups.append((pop_id, poor, good))
    return groups


def run_calibration(
    exclude_study: str | None = None, write: bool = True
) -> FittedParameters:
    """Load everything, fit, and persist the parameters."""
    from venomgap.ingest.compositions import load_compositions
    from venomgap.model.crossreact import IdentityTable
    from venomgap.model.spatial import SpatialCompositionModel, fit_length_scale

    populations = load_compositions()
    identity = IdentityTable.from_json()
    length_scale = fit_length_scale(populations)
    assembly = ModelAssembly(
        populations=populations,
        identity=identity,
        spatial=SpatialCompositionModel(populations, length_scale.ell_km),
    )
    calibration = load_calibration(exclude_study=exclude_study)
    fitted = fit_parameters(assembly, calibration)
    fitted = FittedParameters(
        theta=fitted.theta,
        kappa=fitted.kappa,
        kappa_scale=fitted.kappa_scale,
        c_star=fitted.c_star,
        theta_identifiability=fitted.theta_identifiability,
        theta_unidentified=fitted.theta_unidentified,
        ordinal_violations=fitted.ordinal_violations,
        ordinal_checks=fitted.ordinal_checks,
        ell_km=length_scale.ell_km,
        calibration_rows=fitted.calibration_rows,
        calibration_rmse=fitted.calibration_rmse,
        calibration_studies=fitted.calibration_studies,
        excluded_study=exclude_study,
        converged=fitted.converged,
    )
    if write:
        FITTED_PARAMS_JSON.parent.mkdir(parents=True, exist_ok=True)
        FITTED_PARAMS_JSON.write_text(json.dumps(fitted.model_dump(mode="json"), indent=1))
    return fitted


def load_fitted_parameters(path=FITTED_PARAMS_JSON) -> FittedParameters:
    if not path.exists():
        raise CalibrationError(
            f"{path} does not exist; run `python -m venomgap.cli calibrate` first"
        )
    return FittedParameters.model_validate_json(path.read_text())
