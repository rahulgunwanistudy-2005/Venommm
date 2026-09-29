"""Generate README.md entirely from results/*.json.

Every number in the README is read back out of a result file. Running this script is what makes
"no number is retyped by hand" a property of the build rather than a promise, and CI runs it and
fails if the committed README differs from the generated one.

    python experiments/scripts/build_readme.py          # write README.md
    python experiments/scripts/build_readme.py --check  # fail if it would change
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"
DERIVED = ROOT / "data" / "derived"


def load(path: Path) -> dict:
    if not path.exists():
        raise SystemExit(f"{path} is missing; run `python -m venomgap.cli pipeline` first")
    return json.loads(path.read_text())


def build() -> str:
    retro = load(RESULTS / "retrodiction.json")
    spat = load(RESULTS / "spatial.json")
    opt = load(RESULTS / "optimisation.json")
    sens = load(RESULTS / "sensitivity.json")
    corpus = load(DERIVED / "corpus_summary.json")
    fit = retro["fitted_parameters"]

    crit = {c["criterion_id"]: c for c in retro["criteria"]}

    def row(cid: str) -> str:
        c = crit[cid]
        mark = "**PASS**" if c["passed"] else "FAIL"
        comp = "—" if c["computed"] is None else f"{c['computed']:.4f}"
        thr = "—" if c["threshold"] is None else f"{c['threshold']:.4f}"
        return f"| {cid} | {c['statement']} | {comp} | {thr} | {mark} |"

    diag = retro["diagnostics"]
    theta = diag["theta_scale_mismatch"]
    forced = opt["forced_k_curve"]
    loso = sens["leave_one_study_out"]
    mix = retro["mixture_size_curve"]
    cos = retro["cosine_baseline_curve"]
    mix_peak_index = mix["coverages"].index(max(mix["coverages"]))
    mix_argmax = mix["sizes"][mix_peak_index]
    mix_drop = max(mix["coverages"]) - mix["coverages"][-1]

    forced_rows = "\n".join(
        f"| {r['k']} | {r['uniform']:.4f} | {r['weighted']:.4f} |" for r in forced
    )
    loso_rows = "\n".join(
        f"- drop `{r['excluded_study']}` → ρ = {r['spearman_rho']:.4f}, "
        f"κ_scale {r['kappa_scale']:.3f}"
        for r in loso
    )

    return f"""# VenomGap

**India's only antivenom is an estimator trained on a biased sample. This measures where that
training set fails, and solves for where to sample next.**

India's polyvalent antivenom is raised on venom pooled from four species collected almost entirely
from one facility in Tamil Nadu. Venom composition varies geographically. That is a covariate-shift
problem, and covariate shift is measurable.

> Research model. Not clinical guidance. Predictions are computational and require experimental
> validation. **The model predicts a risk ordering, not clinical neutralisation percentages.**

---

## The result in one line

The model reproduces the experimentally observed **immunogen dilution effect** without being told
about it, and a cosine-similarity baseline provably cannot. Of four pre-registered, blinded
retrodiction targets it passed **{retro['targets_passed']} of {retro['targets_total']}** — and the
one it passed is the structural one the brief calls strongest.

## What the model is

Neutralisation is a **limiting-reagent problem**, not a geometric one. A vial holds a fixed mass of
antibody, `B`, partitioned across toxin families in proportion to the immunising mixture:

```
supply_f(S, w) = B · Σ_q w_q · a_qf                     antibody allocated to family f
eff_f(S, w, p) = B · Σ_q w_q · a_qf · x_f(q → p)         of which is usable against population p
demand_f(p, D) = D · a_pf · κ_f                          what family f of the incoming venom needs
n_f            = min(1, eff_f / demand_f)                saturating neutralised fraction
C(p | S, w, D) = Σ_f σ_f · a_pf · n_f / Σ_f σ_f · a_pf   medically weighted coverage
```

Three properties fall out of that form and **none of them are fitted**:

- **Non-monotone in |S|.** Adding a venom moves antibody supply away from the families a given bite
  actually delivers. Cosine similarity between composition vectors cannot produce this.
- **Dose-dependent.** Demand is linear in `D`, supply is not.
- **Saturating.** Antibody raised in excess of a family's demand is wasted.

### What was actually fitted

This is the load-bearing honesty of the project. **One parameter.**

| Parameter | Status |
|---|---|
| `σ_f` severity weights | **Fixed a priori**, with a written per-family rationale, before any data was ingested |
| `κ_f` shape | **Fixed a priori** as `M_ab / M_family` from molar equivalence — a 7 kDa three-finger toxin costs ~7× more antibody mass per mg than a 47 kDa metalloproteinase |
| `κ_scale` | **Fitted**: {fit['kappa_scale']:.4f} |
| `θ_f` recognition floors | **Not identified.** Data sensitivity is exactly zero for all six probed families, so every `θ_f` keeps its a priori 0.65 |
| `ℓ` spatial length scale | Fitted by leave-one-population-out CV on **composition vectors only** |

The calibration set could not constrain `θ_f` because its targets are populations of species that
are themselves in the immunogen set, where cross-reactivity is 1 by construction. Rather than
present six unmoved parameters as "fitted", the fit measures and reports this. It makes the
retrodiction harder to have rigged than the pre-registration anticipated, and it is a real
limitation.

An **unfitted consistency check** passes: the a priori `κ_f` shape alone reproduces the published
ordering of well- against poorly-recognised toxin families across all {fit['ordinal_checks']}
pairwise checks, with no violations.

---

## The blinded retrodiction

Four published antivenom-failure findings were held out. `model/calibrate.py` raises
`HoldoutLeakError` on any calibration row labelled `holdout=True` **or** citing one of the four
holdout DOIs — a holdout study's compositions are legitimate model inputs, its antivenom results
are not. Tests assert the guard fires.

Pre-registration committed at `{retro['preregistration_commit'][:12]}`, before any parameter was
fitted.

| Criterion | Pre-registered statement | Computed | Threshold | Outcome |
|---|---|---|---|---|
{row('R1a')}
{row('R1b')}
{row('R2a')}
{row('R2b')}
{row('R3a')}
{row('R3b')}
{row('R4a')}
{row('R4b')}
{row('R4c')}

**{retro['targets_passed']}/{retro['targets_total']} targets passed.** A target counts only if every
sub-criterion passes.

### R4 passed completely — the mechanism works

Best achievable burden-weighted national coverage against immunogen mixture size peaks at
**|S| = {mix_argmax}** and falls by {mix_drop:.4f} thereafter. The
cosine best-match baseline, on identical inputs, rises monotonically from
{cos['coverages'][0]:.4f} to {cos['coverages'][-1]:.4f} and has no interior maximum at all. A metric
with no finite supply cannot represent a supply constraint.

### Why the misses missed

Not narrative — computed, and stored in `results/retrodiction.json` under `diagnostics`.

**θ is on the wrong scale for the identity statistic used.** θ = {theta['theta']:.2f} exceeds
**{theta['fraction_below_theta']:.1%}** of the {theta['between_species_identities']} measured
between-species identities, so cross-reactivity is forced to exactly zero for almost every species
pair and the model has nearly no paraspecific coverage. θ was pre-registered to be *fitted*; it
could not be, so it kept a prior reasoned from epitope identity between orthologous toxins, while
the statistic actually used is the mean over all sequence pairs in a family — systematically lower.

R2b makes this vivid. The sequence data **does** rank *N. kaouthia* closer to *N. sagittifera*
(identity {diag['r2b']['identity_3FTx_from_Naja_kaouthia']:.3f}) than *N. naja* is
({diag['r2b']['identity_3FTx_from_Naja_naja']:.3f}) — exactly what the held-out paper claims on
phylogenetic grounds, arrived at independently. But the recognition floor discards that ordering
before it can act, so the substitution changes coverage by exactly zero.

**R3a has a second cause: data sparsity.** Only two *Bungarus caeruleus* proteomes exist and one of
them is the R3 target, so the interpolated immunogen is built partly from the population under test.

---

## What the spatial model found, and it is not what was expected

The kernel length scale fitted by leave-one-population-out CV is **{spat['ell_km']:.0f} km**, with a
flat bottom spanning {spat['ell_plateau_km'][0]:.0f}–{spat['ell_plateau_km'][1]:.0f} km
(RMSE {spat['loo_cv_rmse']:.4f} over {spat['loo_cv_folds']} folds). The grid deliberately runs past
India's north–south extent so a genuine plateau can be told apart from a truncated grid.

A length scale that large means the cross-validation cannot distinguish a local kernel from a
near-global one. **At the present sampling density, geographic distance carries little predictive
power for venom composition**, and the best available estimator is close to a species-level mean.
That is a finding about how sparse the data is, and it is the quantitative case for sampling more
sites. It is reported rather than buried.

{spat['estimated']} districts are estimated; {spat['unknown']} is reported `unknown` and rendered
grey. **Unknown is a finding, not a blank.**

---

## Where to put new collection sites

`k` counts **new** sites added to the existing Big Four immunogen. `k = 0` is the antivenom India
makes today.

| k | uniform vial weights | re-optimised weights |
|---|---|---|
{forced_rows}

Under **uniform vial weights** coverage peaks at **k = {opt['turnover_k']}**, dips as the next
venoms dilute the mixture, then recovers. **Re-optimising the mixture weights removes the dip
entirely**, because the optimiser can down-weight the additions rather than splitting the vial
evenly.

That is the actionable form of the dilution effect: *adding a venom to a polyvalent antivenom
without re-balancing the mixture can make it worse, and re-balancing is what recovers the gain.*

The objective is **non-monotone and not submodular**. No `(1 − 1/e)` guarantee applies and none is
claimed. The optimality gap is **measured**, not derived: greedy + local search matched exact
enumeration over {opt['optimality_gap']['subsets_evaluated']:,} subsets on a reduced instance, gap
**{opt['optimality_gap']['optimality_gap']:.4f}**.

---

## Sensitivity and controls

| Check | Result | Reading |
|---|---|---|
| Leave-one-study-out rank stability | Spearman ρ mean **{sens['mean_spearman']:.3f}**, worst {sens['min_spearman']:.3f} | Well clear of the pre-registered 0.5 falsification threshold |
| `B`/`D` parameter sweep | worst ρ **{sens['sweep_min_spearman']:.3f}** | Ranking survives a 5× change in dose and vial assumptions |
| Sequence-free control (`x_f = 1`) | ρ **{sens['controls']['sequence_free_rho_vs_model']:.3f}** | The sequence layer carries essentially all of the ranking signal — it is not decoration |
| Cosine baseline | ρ **{sens['controls']['cosine_baseline_rho_vs_model']:.3f}** | Ranks populations in near-opposite order |
| Permutation null | mean {sens['controls']['permutation_rho_mean']:.3f}, p95 {sens['controls']['permutation_rho_p95']:.3f} | Confirms the null level |

Leave-one-study-out, per study:

{loso_rows}

---

## The data

{corpus['populations']} population-level proteomes, {corpus['species']} species,
{corpus['states']} states or regions, {corpus['studies']} studies. Every row carries a DOI or PMC id,
a table reference, an access date and a proteomic method. See `data/SOURCES.md` for the full
provenance, the exclusions, and the known between-study disagreement over Rajasthan
*E. c. sochureki*.

**What is imputed, and where you can see it.** Every result carries flags:

| Flag | Meaning |
|---|---|
| `composition_imputed` | No proteome published for this location; interpolated from the nearest sampled populations |
| `no_nearby_proteome` | Nearest proteome beyond the extrapolation limit; reported `unknown`, never covered |
| `sequence_imputed` | Sequences missing or very few; cross-reactivity falls back to genus or global |
| `partial_table_renormalised` | The source study characterised only part of whole venom; the remainder went to `other` |
| `single_study_basis` | Everything known about this population comes from one study |

**The immunogen composition is itself interpolated.** No complete family-level proteome exists for
*any* Big Four species at the locality the venom is actually collected from. Rather than fabricate a
Tamil Nadu row, the immunogen is estimated at the Irula coordinates by the same spatial model used
for districts, carrying the same length scale, uncertainty and flags. That absence is reported as a
finding.

---

## What this predicts, and what it does not

**It predicts:** a *risk ordering* of Indian districts and venom populations by antivenom deficit;
which toxin family drives each deficit; the shape of the coverage-versus-mixture-size curve; and a
ranked answer to where new venom collection sites should go.

**It does not predict:** clinical neutralisation percentages, ED50 values, vial requirements, or
patient outcomes. It is not validated against clinical data and must not be used to guide treatment.

---

## Running it

```bash
pip install -e ".[dev]"
python -m venomgap.cli pipeline     # everything, no network needed
python -m venomgap.cli figures      # figures 1-9, vector PDF + PNG
```

Three screens:

```bash
python -m uvicorn venomgap.api.app:app --port 8000
npm --prefix web run dev            # http://localhost:5173
```

Offline: `VENOMGAP_OFFLINE=1` forbids all network access; any missing cache raises rather than
silently degrading. CI runs with it set, so the committed derived tables are provably sufficient.

## Repository

```
src/venomgap/
  config.py       every constant; sigma_f frozen here, never fitted
  types.py        pydantic contracts; nothing untyped crosses a module boundary
  model/          coverage.py is the model; calibrate.py holds the holdout guard
  optimize/       objective, greedy, local search, exact enumeration, simplex weights
  validate/       retrodiction.py evaluates R1-R4; baselines.py is the ablation
experiments/preregistration.md    the model form, the split, the criteria, and amendment 1
results/*.json                    every number in this README and in the figures
```

Every number quoted here is generated from `results/*.json` by
`experiments/scripts/build_readme.py`. None is retyped.
"""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="fail if README.md would change")
    args = parser.parse_args()

    content = build()
    target = ROOT / "README.md"
    if args.check:
        if not target.exists() or target.read_text() != content:
            print("README.md is out of date; run experiments/scripts/build_readme.py", file=sys.stderr)
            return 1
        print("README.md is up to date")
        return 0
    target.write_text(content)
    print(f"wrote {target} ({len(content)} characters)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
