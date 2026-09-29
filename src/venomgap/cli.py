"""VenomGap command line: every step of the pipeline, runnable one at a time or end to end.

    python -m venomgap.cli validate-data      # M2 gate: curated corpus honours its contract
    python -m venomgap.cli fetch-sequences    # needs network; caches UniProt/ToxProt
    python -m venomgap.cli fetch-occurrences  # needs network; caches GBIF
    python -m venomgap.cli build-crossreact   # sequence identity table -> derived/
    python -m venomgap.cli calibrate          # fit theta_f, kappa_f on the calibration split
    python -m venomgap.cli retrodict          # evaluate R1-R4, write results/retrodiction.json
    python -m venomgap.cli spatial            # fit ell by LOO-CV, write district results
    python -m venomgap.cli optimise           # coverage-vs-k curve, siting solutions
    python -m venomgap.cli sensitivity        # leave-one-study-out, sweeps, cosine ablation
    python -m venomgap.cli figures            # regenerate figures 1-9
    python -m venomgap.cli pipeline           # everything that does not need the network
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from collections.abc import Sequence

from venomgap.config import DERIVED_DIR, RESULTS_DIR

logger = logging.getLogger("venomgap")


def _configure_logging(verbose: bool) -> None:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(levelname)-7s %(name)s: %(message)s",
        stream=sys.stderr,
    )


def cmd_validate_data(_: argparse.Namespace) -> int:
    from venomgap.ingest.compositions import corpus_summary, load_compositions

    populations = load_compositions()
    summary = corpus_summary(populations)
    DERIVED_DIR.mkdir(parents=True, exist_ok=True)
    (DERIVED_DIR / "corpus_summary.json").write_text(json.dumps(summary, indent=1))

    print(f"populations      : {summary['populations']}")
    print(f"species          : {summary['species']}")
    print(f"states/regions   : {summary['states']}")
    print(f"distinct studies : {summary['studies']}")
    print(f"holdout          : {len(summary['holdout_populations'])}")
    print(f"flagged          : {len(summary['flagged_populations'])}")
    for pop_id in summary["flagged_populations"]:
        population = next(p for p in populations if p.pop_id == pop_id)
        print(f"    {pop_id}: {', '.join(population.flags)}")
    return 0


def cmd_fetch_sequences(args: argparse.Namespace) -> int:
    from venomgap.ingest.compositions import load_compositions
    from venomgap.ingest.uniprot import cache_summary, fetch_species_sequences

    populations = load_compositions()
    species = sorted({p.species for p in populations})
    # Nominate-species fallback for subspecies: UniProt indexes by species, so
    # "Echis carinatus sochureki" also pulls "Echis carinatus".
    species += sorted({" ".join(s.split()[:2]) for s in species if len(s.split()) > 2})
    for name in sorted(set(species)):
        sequences = fetch_species_sequences(name, refresh=args.refresh)
        families = sorted({s.family for s in sequences})
        print(f"{name:<34} {len(sequences):>4} sequences  {families}")
    summary = cache_summary()
    (DERIVED_DIR / "sequence_cache_summary.json").write_text(json.dumps(summary, indent=1))
    return 0


def cmd_fetch_occurrences(args: argparse.Namespace) -> int:
    from venomgap.ingest.compositions import load_compositions
    from venomgap.ingest.gbif import fetch_species_occurrences

    populations = load_compositions()
    for name in sorted({p.species for p in populations}):
        records = fetch_species_occurrences(name, refresh=args.refresh)
        print(f"{name:<34} {len(records):>5} georeferenced Indian occurrences")
    return 0


def cmd_build_crossreact(_: argparse.Namespace) -> int:
    from venomgap.ingest.compositions import load_compositions
    from venomgap.ingest.uniprot import fetch_species_sequences
    from venomgap.model.crossreact import build_identity_table

    populations = load_compositions()
    names = sorted({p.species for p in populations})
    names += sorted({" ".join(s.split()[:2]) for s in names if len(s.split()) > 2})
    sequences = {}
    for name in sorted(set(names)):
        try:
            sequences[name] = fetch_species_sequences(name)
        except Exception as exc:
            logger.warning("no sequences for %s: %s", name, exc)
    table = build_identity_table(sequences)
    table.to_json()
    print(f"species with sequences : {len(table.species)}")
    print(f"within-species values  : {len(table.within_species)}")
    print(f"between-species values : {len(table.values)}")
    return 0


def cmd_calibrate(_: argparse.Namespace) -> int:
    from venomgap.model.calibrate import run_calibration

    fitted = run_calibration()
    print(f"calibration rows : {fitted.calibration_rows}")
    print(f"RMSE             : {fitted.calibration_rmse:.4f}")
    print(f"kappa_scale      : {fitted.kappa_scale:.4f}")
    print("theta_f          :")
    for family, value in fitted.theta.items():
        print(f"    {family:<6} {value:.3f}")
    return 0


def cmd_retrodict(_: argparse.Namespace) -> int:
    from venomgap.validate.retrodiction import run_retrodiction

    outcome = run_retrodiction()
    passed, total = outcome.targets_passed, outcome.targets_total
    print(f"\n=== RETRODICTION: {passed}/{total} targets passed ===")
    for criterion in outcome.criteria:
        mark = "PASS" if criterion.passed else "FAIL"
        computed = "n/a" if criterion.computed is None else f"{criterion.computed:.4f}"
        print(f"  [{mark}] {criterion.criterion_id:<5} {computed:>9}  {criterion.detail}")
    return 0


def cmd_spatial(_: argparse.Namespace) -> int:
    from venomgap.pipeline import run_spatial

    result = run_spatial()
    print(f"fitted length scale : {result['ell_km']:.0f} km")
    print(f"districts estimated : {result['estimated']}")
    print(f"districts unknown   : {result['unknown']}")
    return 0


def cmd_optimise(_: argparse.Namespace) -> int:
    from venomgap.pipeline import run_optimisation

    result = run_optimisation()
    print(f"turnover at k = {result['turnover_k']}")
    print(f"optimality gap  = {result['optimality_gap']}")
    return 0


def cmd_sensitivity(_: argparse.Namespace) -> int:
    from venomgap.pipeline import run_sensitivity

    result = run_sensitivity()
    print(f"leave-one-study-out mean Spearman rho : {result['mean_spearman']:.3f}")
    return 0


def cmd_figures(_: argparse.Namespace) -> int:
    from venomgap.figures import render_all

    paths = render_all()
    for path in paths:
        print(path)
    return 0


def cmd_pipeline(_: argparse.Namespace) -> int:
    """Everything that can run from the committed derived tables, no network."""
    steps = [
        ("validate-data", cmd_validate_data),
        ("calibrate", cmd_calibrate),
        ("spatial", cmd_spatial),
        ("optimise", cmd_optimise),
        ("retrodict", cmd_retrodict),
        ("sensitivity", cmd_sensitivity),
        ("figures", cmd_figures),
    ]
    namespace = argparse.Namespace(refresh=False)
    for name, handler in steps:
        print(f"\n{'=' * 70}\n== {name}\n{'=' * 70}")
        code = handler(namespace)
        if code != 0:
            logger.error("step %s failed with code %d", name, code)
            return code
    print(f"\nresults written to {RESULTS_DIR}")
    return 0


COMMANDS = {
    "validate-data": cmd_validate_data,
    "fetch-sequences": cmd_fetch_sequences,
    "fetch-occurrences": cmd_fetch_occurrences,
    "build-crossreact": cmd_build_crossreact,
    "calibrate": cmd_calibrate,
    "retrodict": cmd_retrodict,
    "spatial": cmd_spatial,
    "optimise": cmd_optimise,
    "sensitivity": cmd_sensitivity,
    "figures": cmd_figures,
    "pipeline": cmd_pipeline,
}


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="venomgap", description=__doc__)
    parser.add_argument("command", choices=sorted(COMMANDS))
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument(
        "--refresh",
        action="store_true",
        help="re-download cached network resources instead of using the cache",
    )
    args = parser.parse_args(argv)
    _configure_logging(args.verbose)
    return COMMANDS[args.command](args)


if __name__ == "__main__":
    raise SystemExit(main())
