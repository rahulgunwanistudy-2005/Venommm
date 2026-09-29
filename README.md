# VenomGap

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
retrodiction targets it passed **1 of 4** — and the
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
| `κ_scale` | **Fitted**: 0.7839 |
| `θ_f` recognition floors | **Not identified.** Data sensitivity is exactly zero for all six probed families, so every `θ_f` keeps its a priori 0.65 |
| `ℓ` spatial length scale | Fitted by leave-one-population-out CV on **composition vectors only** |

The calibration set could not constrain `θ_f` because its targets are populations of species that
are themselves in the immunogen set, where cross-reactivity is 1 by construction. Rather than
present six unmoved parameters as "fitted", the fit measures and reports this. It makes the
retrodiction harder to have rigged than the pre-registration anticipated, and it is a real
limitation.

An **unfitted consistency check** passes: the a priori `κ_f` shape alone reproduces the published
ordering of well- against poorly-recognised toxin families across all 8
pairwise checks, with no violations.

---

## The blinded retrodiction

Four published antivenom-failure findings were held out. `model/calibrate.py` raises
`HoldoutLeakError` on any calibration row labelled `holdout=True` **or** citing one of the four
holdout DOIs — a holdout study's compositions are legitimate model inputs, its antivenom results
are not. Tests assert the guard fires.

Pre-registration committed at `a3cb96c963bc`, before any parameter was
fitted.

| Criterion | Pre-registered statement | Computed | Threshold | Outcome |
|---|---|---|---|---|
| R1a | Each north-west E. c. sochureki population ranks in the top decile of predicted deficit among curated Indian populations; passes on a strict majority. | 0.0000 | 2.0000 | FAIL |
| R1b | SVMP is the single worst limiting family for each such population; passes on a strict majority. | 3.0000 | 2.0000 | **PASS** |
| R2a | Naja sagittifera ranks in the top decile of predicted deficit. | 0.1364 | 0.1000 | FAIL |
| R2b | Substituting Naja kaouthia for Naja naja in the immunogen set, holding everything else fixed, increases predicted coverage of Naja sagittifera by at least 0.05. | 0.0000 | 0.0500 | FAIL |
| R3a | North Indian Bungarus caeruleus, a Big Four species, has a predicted deficit in the top tertile of curated Indian populations. | 0.7273 | 0.3333 | FAIL |
| R3b | Its predicted deficit exceeds that of the southern reference B. caeruleus population nearest the immunogen source by at least 0.10 absolute. | 0.1456 | 0.1000 | **PASS** |
| R4a | Mean burden-weighted national coverage as a function of immunogen mixture size, each size at its own best subset, is non-monotone with an interior maximum. | 4.0000 | — | **PASS** |
| R4b | The fall from the maximum to the largest mixture is at least 0.02 absolute coverage, so the turnover is not numerical noise. | 0.0268 | 0.0200 | **PASS** |
| R4c | The cosine-similarity baseline, on identical inputs, is monotone non-decreasing in mixture size and therefore cannot reproduce R4a. | 6.0000 | — | **PASS** |

**1/4 targets passed.** A target counts only if every
sub-criterion passes.

### R4 passed completely — the mechanism works

Best achievable burden-weighted national coverage against immunogen mixture size peaks at
**|S| = 4** and falls by 0.0268 thereafter. The
cosine best-match baseline, on identical inputs, rises monotonically from
0.7288 to 0.9489 and has no interior maximum at all. A metric
with no finite supply cannot represent a supply constraint.

### Why the misses missed

Not narrative — computed, and stored in `results/retrodiction.json` under `diagnostics`.

**θ is on the wrong scale for the identity statistic used.** θ = 0.65 exceeds
**92.7%** of the 82 measured
between-species identities, so cross-reactivity is forced to exactly zero for almost every species
pair and the model has nearly no paraspecific coverage. θ was pre-registered to be *fitted*; it
could not be, so it kept a prior reasoned from epitope identity between orthologous toxins, while
the statistic actually used is the mean over all sequence pairs in a family — systematically lower.

R2b makes this vivid. The sequence data **does** rank *N. kaouthia* closer to *N. sagittifera*
(identity 0.388) than *N. naja* is
(0.278) — exactly what the held-out paper claims on
phylogenetic grounds, arrived at independently. But the recognition floor discards that ordering
before it can act, so the substitution changes coverage by exactly zero.

**R3a has a second cause: data sparsity.** Only two *Bungarus caeruleus* proteomes exist and one of
them is the R3 target, so the interpolated immunogen is built partly from the population under test.

---

## What the spatial model found, and it is not what was expected

The kernel length scale fitted by leave-one-population-out CV is **1600 km**, with a
flat bottom spanning 1000–10000 km
(RMSE 0.0726 over 29 folds). The grid deliberately runs past
India's north–south extent so a genuine plateau can be told apart from a truncated grid.

A length scale that large means the cross-validation cannot distinguish a local kernel from a
near-global one. **At the present sampling density, geographic distance carries little predictive
power for venom composition**, and the best available estimator is close to a species-level mean.
That is a finding about how sparse the data is, and it is the quantitative case for sampling more
sites. It is reported rather than buried.

593 districts are estimated; 1 is reported `unknown` and rendered
grey. **Unknown is a finding, not a blank.**

---

## Where to put new collection sites

`k` counts **new** sites added to the existing Big Four immunogen. `k = 0` is the antivenom India
makes today.

| k | uniform vial weights | re-optimised weights |
|---|---|---|
| 0 | 0.6615 | 0.6615 |
| 1 | 0.7064 | 0.7135 |
| 2 | 0.6954 | 0.7438 |
| 3 | 0.6900 | 0.7438 |
| 4 | 0.7075 | 0.7513 |
| 5 | 0.7239 | 0.7535 |
| 6 | 0.7348 | 0.7614 |

Under **uniform vial weights** coverage peaks at **k = 1**, dips as the next
venoms dilute the mixture, then recovers. **Re-optimising the mixture weights removes the dip
entirely**, because the optimiser can down-weight the additions rather than splitting the vial
evenly.

That is the actionable form of the dilution effect: *adding a venom to a polyvalent antivenom
without re-balancing the mixture can make it worse, and re-balancing is what recovers the gain.*

The objective is **non-monotone and not submodular**. No `(1 − 1/e)` guarantee applies and none is
claimed. The optimality gap is **measured**, not derived: greedy + local search matched exact
enumeration over 6,195 subsets on a reduced instance, gap
**0.0000**.

---

## Sensitivity and controls

| Check | Result | Reading |
|---|---|---|
| Leave-one-study-out rank stability | Spearman ρ mean **0.994**, worst 0.989 | Well clear of the pre-registered 0.5 falsification threshold |
| `B`/`D` parameter sweep | worst ρ **0.649** | Ranking survives a 5× change in dose and vial assumptions |
| Sequence-free control (`x_f = 1`) | ρ **0.052** | The sequence layer carries essentially all of the ranking signal — it is not decoration |
| Cosine baseline | ρ **-0.249** | Ranks populations in near-opposite order |
| Permutation null | mean -0.018, p95 0.413 | Confirms the null level |

Leave-one-study-out, per study:

- drop `Kalita2017Echis` → ρ = 1.0000, κ_scale 0.786
- drop `SenjiLaxme2021Daboia` → ρ = 0.9887, κ_scale 0.938
- drop `SenjiLaxme2021Naja` → ρ = 0.9944, κ_scale 0.830

---

## The data

29 population-level proteomes, 9 species,
14 states or regions, 6 studies. Every row carries a DOI or PMC id,
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
