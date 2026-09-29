# VenomGap — build plan

## Goal
A mechanistic, stoichiometric model of antivenom coverage across India that (a) reproduces the
immunogen-dilution effect without being told about it, (b) retrodicts four pre-registered,
held-out published failure findings, and (c) solves for where to site *k* new venom collection
centres to maximise burden-weighted national coverage.

## Build order is strict
M1 (coverage model on synthetic inputs) before any data work. If the model does not reproduce
the dilution effect analytically, nothing downstream matters.

---

### M0 — Scaffold, types, config, CI
- [x] Directory skeleton per §5 of the brief
- [x] `pyproject.toml`, pinned `requirements.lock`
- [x] `config.py` — every constant, fixed `sigma_f`, `kappa_f` priors, `B`, `D`
- [x] `types.py` — pydantic contracts from §6, nothing untyped crosses a module boundary
- [x] `errors.py` — typed error hierarchy
- [x] ruff + mypy config, GitHub Actions CI
- [x] Gate: `pytest` green, `mypy src` clean, `ruff check` clean

### M1 — The coverage model (`model/coverage.py`)
- [x] `supply_f`, `eff_f`, `demand_f`, `n_f`, `C(p | S, w, D)` exactly as §1.3
- [x] Analytic unit tests, not smoke tests:
  - [x] identical `p` and singleton `S={p}` with sufficient budget -> coverage 1.0
  - [x] disjoint family support -> coverage ~ 0
  - [x] **adding a dilutant venom lowers coverage** (the R4 mechanism)
  - [x] dose-response monotone decreasing in `D`
  - [x] saturation: doubling `B` past sufficiency does not change coverage
  - [x] invariance: coverage independent of units of `B` and `D` at fixed `B/D`
- [x] Gate: all six analytic properties pass

### M2 — Curated composition corpus
- [x] `data/curated/compositions.csv` — every row carries DOI + table + accessed date + method
- [x] >= 25 Indian populations across >= 6 species
- [x] Schema validation: keys subset of FAMILIES, each row sums to 1 +/- 1e-6
- [x] `data/SOURCES.md` with licence + access date + exact query/DOI per source
- [x] Gate: `venomgap validate-data` passes; row count and species count asserted in tests

### M3 — UniProt/ToxProt ingest + cross-reactivity
- [x] `ingest/uniprot.py` — REST fetch, on-disk cache, family mapping from `protein_families`
- [x] `model/crossreact.py` — pairwise identity -> `x_f(q->p)` with fitted floor `theta_f`
- [x] Gate: within-species cross-reactivity > between-genus; cache makes run reproducible offline

### M4 — Calibration
- [x] `data/curated/antivenomics_calibration.csv` (non-holdout) and `holdout_targets.csv`
- [x] `model/calibrate.py` fits `theta_f`, `kappa_f`; **raises on any `holdout=True` row**
- [x] Gate: a test asserts `HoldoutLeakError` is raised

### M5 — Retrodiction run
- [x] `experiments/preregistration.md` committed BEFORE fitting, hash recorded
- [x] `validate/retrodiction.py` evaluates R1-R4 automatically
- [x] Gate: `results/retrodiction.json` written; outcome reported whatever it is

### M6 — Spatial model
- [x] `model/spatial.py` — Gaussian-kernel interpolation, `ell` by leave-one-population-out CV
- [x] Uncertainty grows with distance; districts past cut-off -> `status="unknown"` (grey, never green)
- [x] Gate: fitted `ell` reported in km; unknown districts identified and counted

### M7 — Optimiser
- [x] `optimize/objective.py`, `greedy.py`, `local_search.py`, `exact.py`, `weights.py`
- [x] Coverage-vs-`k` curve with turnover point
- [x] Optimality gap measured against exact enumeration on a reduced instance
- [x] Gate: gap reported; non-submodularity stated, no (1-1/e) claim made

### M8 — API + atlas web app
- [x] FastAPI serving precomputed district results, siting solutions, mixture evaluation
- [x] Vite + React + TS: Atlas, District detail, Mixture explorer
- [x] Gate: `k` slider 0..6 redraws; district click shows deficit, limiting families, flags

### M9 — Sensitivity + ablation
- [x] Leave-one-study-out rank stability (Spearman rho)
- [x] Cosine-similarity baseline: must fail R4 while the stoichiometric model passes
- [x] Parameter sweeps over `B`, `D`
- [x] Gate: `results/sensitivity.json`, `results/ablation.json`

### M10 — Figures, paper skeleton, demo
- [x] Figures 1-9 from `venomgap figures`, vector PDF + PNG, units on every axis
- [x] `demo/demo_script.md`, offline fallback bundle
- [x] Gate: full pipeline reproducible offline from `data/curated/`

---

## Risks (the §10 list — how this loses)
- **Cosine similarity as the coverage metric.** Cannot produce the dilution effect. Mitigation:
  cosine exists only in `validate/baselines.py` as an ablation that is *shown to fail*.
- **Tuning `sigma_f`, `theta_f`, `kappa_f` after seeing R1-R4.** Mitigation: `sigma_f` frozen in
  `config.py` with a literature rationale per family; `theta_f`/`kappa_f` fitted only through
  `calibrate.py`, which raises on holdout rows; pre-registration committed before fitting.
- **Colouring data-free districts green.** Mitigation: `status` is a required field on
  `DistrictResult`; the web legend has an explicit grey `unknown` class; a test asserts unknown
  districts are never assigned a deficit colour.
- **Claiming clinical neutralisation percentages.** Mitigation: README, figures and API all say
  the output is a *risk ordering*; retrodiction criteria are rank-based, not value-based.
- **Treating GBIF occurrence counts as abundance.** Mitigation: occurrences are used only to
  *exclude* districts from a species range; a comment and a test document this.
- **Shipping a map with no mechanism.** Mitigation: every district payload carries
  `limiting_families`, so the map always answers "why", not only "where".
- **Silently dropping a study that does not fit.** Mitigation: ingest counts rows in and rows out
  and raises if they differ; every exclusion is explicit in `SOURCES.md`.
- **"AI-powered" language.** Mitigation: none present; this is a mechanistic model.

## Additional risks specific to the data situation
- **Proteome numbers live in PDFs.** Family-level percentages are transcribed by hand from the
  published tables. Every row carries DOI + table. Rows whose abundances were re-normalised from
  a partial table (families reported summing to < 100%) carry the residual in `other` and are
  flagged. Rows reconstructed from a figure rather than a table carry
  `method="other"` and are listed explicitly in `SOURCES.md`.
- **Krait (`Bungarus caeruleus`) has few reviewed UniProt entries.** Cross-reactivity for KSPI and
  3FTx in krait falls back to genus consensus with `sequence_imputed=True`, which is surfaced in
  every affected `CoverageResult.flags`.
- **District geometry.** If an open district GeoJSON cannot be retrieved, the atlas degrades to
  state level and says so; the model itself is district-resolved regardless.

## Done criteria
- [x] M0-M10 gates passed and evidenced here
- [x] `results/retrodiction.json` with pre-registered criteria and the honest outcome
- [x] Figures 1-9 regenerate from one command
- [x] Demo runs offline
- [x] README states what the model predicts, what it does not, and which inputs are imputed

---

## Review

### What was built and verified

All eleven gates M0–M10 passed, each demonstrated rather than claimed. Final state:
ruff clean, mypy strict clean on 30 source files, **95 tests green**, `tsc --noEmit` clean,
`vite build` clean, and the full pipeline runs end to end with `VENOMGAP_OFFLINE=1` and exits 0.

**The headline result.** The coverage model reproduces the immunogen dilution effect from a finite
antibody budget with nothing fitted to produce it, and the cosine best-match baseline — run in its
strongest form, not as a straw man — is monotone non-decreasing and provably cannot. R4 passed on
all three sub-criteria.

**The retrodiction: 1 of 4.** R4 passed completely. R1b and R3b passed as sub-criteria. R1a, R2a,
R2b and R3a failed, and the failures are diagnosed with computed numbers stored in
`results/retrodiction.json`, not with narrative. The dominant cause is a single identified scale
mismatch: theta = 0.65 exceeds 92.7% of the 82 measured between-species sequence identities, so
cross-reactivity is forced to zero almost everywhere. theta was pre-registered to be fitted and
turned out to be unidentifiable from the available calibration data.

**What was actually fitted: one parameter.** kappa_scale = 0.784. sigma_f, the kappa_f shape and
every theta_f retain their a priori values. An unfitted consistency check passes: the a priori
kappa shape alone reproduces the published ordering of well- against poorly-recognised toxin
families across all 8 pairwise checks.

**Notable negative result, reported rather than buried.** The spatial length scale fitted by
leave-one-population-out CV is 1600 km with a plateau running past the length of the country. At
the present sampling density geographic distance carries almost no predictive power for venom
composition. That is the quantitative case for sampling more sites.

### What was deferred or is honestly weak

- **theta_f is unvalidated.** The calibration set cannot constrain it. The fix is a cross-reactivity
  statistic on the same scale as the prior — best-match identity across a family rather than the
  mean over all pairs — but changing it after seeing the retrodiction would be post-hoc, so it is
  recorded as the primary diagnosis and left for a pre-registered follow-up.
- **Only two *Bungarus caeruleus* proteomes exist**, one of them the R3 target, so the interpolated
  immunogen is partly built from the population under test. No amount of modelling fixes this; it
  needs a South Indian krait proteome.
- **District burden allocation uses district area as a rural-population proxy**, because 2011 census
  district rural populations were not retrievable without an API key. Declared on every row and
  swept in the sensitivity analysis.
- **The atlas map is inline SVG from district centroids, not a MapLibre tile choropleth.** A tile
  server is a network dependency and the demo must run offline. Same information, different
  rendering; stated in the commit and in the demo script.
- **No live polygon geometry is committed** — the 34 MB GeoJSON is fetched and reduced to centroids
  and areas rather than stored.
