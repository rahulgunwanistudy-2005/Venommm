"""M4 gate: the holdout guard is an assertion, and the fit reports what it could not identify."""

from __future__ import annotations

import csv

import pytest

from venomgap.config import CALIBRATION_CSV, FAMILIES, HOLDOUT_STUDY_DOIS, THETA_PRIOR
from venomgap.errors import HoldoutLeakError
from venomgap.model.calibrate import load_calibration, load_fitted_parameters


@pytest.fixture(scope="module")
def calibration_rows() -> tuple[list[dict[str, str]], list[str]]:
    with CALIBRATION_CSV.open() as fh:
        rows = list(csv.DictReader(fh))
    return rows, list(rows[0])


def _write(tmp_path, rows, header):  # type: ignore[no-untyped-def]
    path = tmp_path / "calibration.csv"
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=header)
        writer.writeheader()
        writer.writerows(rows)
    return path


# --------------------------------------------------------------------------------------
# The guard. This is the test the whole blinded design rests on.
# --------------------------------------------------------------------------------------


def test_calibration_refuses_a_row_labelled_holdout(tmp_path, calibration_rows) -> None:  # type: ignore[no-untyped-def]
    rows, header = calibration_rows
    poisoned = [dict(r) for r in rows]
    poisoned[0]["holdout"] = "true"
    with pytest.raises(HoldoutLeakError, match="holdout=True"):
        load_calibration(_write(tmp_path, poisoned, header))


def test_calibration_refuses_a_row_citing_a_holdout_study(tmp_path, calibration_rows) -> None:  # type: ignore[no-untyped-def]
    """A holdout study's compositions are inputs; its antivenom results are not."""
    rows, header = calibration_rows
    poisoned = [dict(r) for r in rows]
    poisoned[0]["doi"] = "10.1371/journal.pntd.0007899"  # the R3 source
    with pytest.raises(HoldoutLeakError, match="holdout studies"):
        load_calibration(_write(tmp_path, poisoned, header))


@pytest.mark.parametrize("doi", sorted(HOLDOUT_STUDY_DOIS))
def test_every_holdout_study_is_refused(tmp_path, calibration_rows, doi: str) -> None:  # type: ignore[no-untyped-def]
    rows, header = calibration_rows
    poisoned = [dict(r) for r in rows]
    poisoned[0]["doi"] = doi
    with pytest.raises(HoldoutLeakError):
        load_calibration(_write(tmp_path, poisoned, header))


def test_the_committed_calibration_file_is_clean() -> None:
    calibration = load_calibration()
    assert calibration.rows > 0
    for frame in (
        calibration.potency,
        calibration.censored,
        calibration.poorly_recognised,
        calibration.well_recognised,
    ):
        assert not frame["doi"].isin(HOLDOUT_STUDY_DOIS).any()
        assert not frame["holdout"].astype(str).str.lower().isin({"true", "1", "yes"}).any()


def test_leave_one_study_out_actually_removes_a_study() -> None:
    full = load_calibration()
    assert len(full.studies) >= 2
    reduced = load_calibration(exclude_study=full.studies[0])
    assert full.studies[0] not in reduced.studies
    assert reduced.rows < full.rows


# --------------------------------------------------------------------------------------
# What the fit produced, and what it honestly could not
# --------------------------------------------------------------------------------------


def test_fitted_parameters_cover_the_whole_vocabulary() -> None:
    fitted = load_fitted_parameters()
    assert set(fitted.theta) == set(FAMILIES)
    assert set(fitted.kappa) == set(FAMILIES)
    assert all(v > 0.0 for v in fitted.kappa.values())


def test_kappa_shape_is_a_priori_and_only_its_scale_is_fitted() -> None:
    """kappa_f must stay proportional to 1 / M_family. If the shape ever becomes fitted, the
    mechanistic molar-equivalence argument is gone and this test should fail loudly."""
    from venomgap.config import KAPPA_PRIOR

    fitted = load_fitted_parameters()
    ratios = [fitted.kappa[f] / KAPPA_PRIOR[f] for f in FAMILIES]
    assert max(ratios) - min(ratios) < 1e-9
    assert ratios[0] == pytest.approx(fitted.kappa_scale)


def test_unidentified_thetas_keep_their_a_priori_value() -> None:
    """An unidentified parameter sitting at its prior is fine. Presenting it as fitted is not."""
    fitted = load_fitted_parameters()
    for family in fitted.theta_unidentified:
        assert fitted.theta[family] == pytest.approx(THETA_PRIOR), (
            f"{family} is reported as unidentified but does not sit at the a priori floor"
        )


def test_identifiability_is_reported_for_every_probed_family() -> None:
    from venomgap.model.calibrate import FITTED_THETA_FAMILIES

    fitted = load_fitted_parameters()
    assert set(fitted.theta_identifiability) == set(FITTED_THETA_FAMILIES)
    assert all(v >= 0.0 for v in fitted.theta_identifiability.values())


def test_a_priori_kappa_shape_reproduces_published_recognition_ordering() -> None:
    """An unfitted consistency check: the molar-equivalence kappa shape should, on its own, rank
    the families a published antivenomics study found poorly recognised below the others."""
    fitted = load_fitted_parameters()
    assert fitted.ordinal_checks > 0
    assert fitted.ordinal_violations == (), (
        f"the a priori kappa shape violates the published ordering: {fitted.ordinal_violations}"
    )


def test_fit_converged_and_is_recorded() -> None:
    fitted = load_fitted_parameters()
    assert fitted.converged
    assert fitted.ell_km is not None and fitted.ell_km > 0
    assert fitted.calibration_rmse >= 0.0
    assert fitted.excluded_study is None
