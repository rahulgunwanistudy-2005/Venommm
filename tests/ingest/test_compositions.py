"""M2 gate: the curated composition corpus honours its contract.

These tests assert on the committed data, not on synthetic fixtures, because the corpus is a
deliverable in its own right. If someone edits a row badly, this is what catches it.
"""

from __future__ import annotations

import csv

import pytest

from venomgap.config import COMPOSITIONS_CSV, FAMILIES
from venomgap.errors import DataValidationError, ProvenanceError
from venomgap.ingest.compositions import (
    MAX_RESIDUAL,
    REQUIRED_COLUMNS,
    SOURCE_FAMILY_MAP,
    corpus_summary,
    load_compositions,
)

MIN_POPULATIONS = 25
MIN_SPECIES = 6


@pytest.fixture(scope="module")
def populations() -> list:
    return load_compositions()


def test_corpus_meets_the_m2_gate(populations: list) -> None:
    indian = [p for p in populations if p.country == "India"]
    assert len(populations) >= MIN_POPULATIONS
    assert len(indian) >= MIN_POPULATIONS - 4
    assert len({p.species for p in populations}) >= MIN_SPECIES


def test_every_composition_lies_on_the_simplex(populations: list) -> None:
    for p in populations:
        total = sum(p.composition.values())
        assert total == pytest.approx(1.0, abs=1e-9), f"{p.pop_id} sums to {total}"
        assert set(p.composition) <= set(FAMILIES), f"{p.pop_id} has an off-vocabulary family"
        assert all(v >= 0.0 for v in p.composition.values())


def test_every_row_carries_provenance(populations: list) -> None:
    """A reviewer will check this. Every number must be traceable to a table in a paper."""
    for p in populations:
        assert p.provenance.doi
        assert p.provenance.table
        assert p.provenance.accessed
        assert p.provenance.method in {
            "LC-MS/MS",
            "RP-HPLC+MS",
            "transcriptome-informed",
            "other",
        }


def test_pop_ids_are_unique(populations: list) -> None:
    ids = [p.pop_id for p in populations]
    assert len(set(ids)) == len(ids)


def test_holdout_populations_are_exactly_the_preregistered_targets(populations: list) -> None:
    """The holdout set is fixed by the pre-registration. It must not drift."""
    holdout = {p.pop_id for p in populations if p.holdout}
    assert "Naja_sagittifera__Andaman" in holdout, "R2 target missing"
    assert "Bungarus_caeruleus__Punjab" in holdout, "R3 target missing"
    sochureki = {p for p in holdout if p.startswith("Echis_carinatus_sochureki__")}
    assert len(sochureki) >= 2, "R1 needs the north-west E. c. sochureki populations held out"
    # Nothing outside the R1-R3 target species is held out.
    for pop_id in holdout:
        assert pop_id.startswith(
            ("Naja_sagittifera__", "Bungarus_caeruleus__Punjab", "Echis_carinatus_sochureki__")
        )


def test_no_row_exceeds_the_residual_limit(populations: list) -> None:
    """A study that characterised less than 55% of whole venom cannot anchor a composition."""
    for p in populations:
        other = p.composition.get("other", 0.0)
        assert other <= MAX_RESIDUAL + 1e-9, f"{p.pop_id} assigns {other:.2%} to other"


def test_partial_tables_are_flagged(populations: list) -> None:
    """Whenever the residual is material the row must say so; silence here would be dishonest."""
    for p in populations:
        if p.composition.get("other", 0.0) > 0.10:
            assert "partial_table_renormalised" in p.flags, f"{p.pop_id} hides a large residual"


def test_the_r1_studies_genuinely_disagree(populations: list) -> None:
    """Pre-registration amendment A1.1 rests on this. If the corpus stops showing the conflict,
    the amendment's reasoning no longer applies and someone must revisit it."""
    by_id = {p.pop_id: p for p in populations}
    kumar = by_id["Echis_carinatus_sochureki__Sam"].composition.get("SVMP", 0.0)
    laxme = by_id["Echis_carinatus_sochureki__RajasthanPooled"].composition.get("SVMP", 0.0)
    assert kumar > 0.4, "the 2026 proteome should be SVMP-dominant"
    assert laxme < 0.1, "the 2019 proteome should report SVMP as a minor family"


def test_every_abundance_column_maps_into_the_vocabulary() -> None:
    with COMPOSITIONS_CSV.open() as fh:
        header = next(csv.reader(fh))
    for column in header:
        if column.startswith("pct_"):
            assert column.removeprefix("pct_").lower() in SOURCE_FAMILY_MAP
    for column in REQUIRED_COLUMNS:
        assert column in header, f"required column {column} missing from the CSV"


def test_summary_counts_are_computed_not_asserted(populations: list) -> None:
    summary = corpus_summary(populations)
    assert summary["populations"] == len(populations)
    assert summary["species"] == len({p.species for p in populations})
    assert summary["studies"] == len({p.provenance.doi for p in populations})


def test_missing_provenance_raises(tmp_path) -> None:  # type: ignore[no-untyped-def]
    """Provenance is enforced, not requested."""
    with COMPOSITIONS_CSV.open() as fh:
        rows = list(csv.DictReader(fh))
        header = list(rows[0])
    rows[0]["doi"] = ""
    broken = tmp_path / "broken.csv"
    with broken.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=header)
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(ProvenanceError, match="missing doi"):
        load_compositions(broken)


def test_composition_not_summing_is_rejected(tmp_path) -> None:  # type: ignore[no-untyped-def]
    with COMPOSITIONS_CSV.open() as fh:
        rows = list(csv.DictReader(fh))
        header = list(rows[0])
    for column in header:
        if column.startswith("pct_"):
            rows[0][column] = ""
    rows[0]["pct_pla2"] = "5.0"  # only 5% of whole venom characterised
    broken = tmp_path / "thin.csv"
    with broken.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=header)
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(DataValidationError, match="residual exceeds"):
        load_compositions(broken)
