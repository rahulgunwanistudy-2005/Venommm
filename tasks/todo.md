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
- [ ] Directory skeleton per §5 of the brief
- [ ] `pyproject.toml`, pinned `requirements.lock`
- [ ] `config.py` — every constant, fixed `sigma_f`, `kappa_f` priors, `B`, `D`
- [ ] `types.py` — pydantic contracts from §6, nothing untyped crosses a module boundary
- [ ] `errors.py` — typed error hierarchy
- [ ] ruff + mypy config, GitHub Actions CI
- [ ] Gate: `pytest` green, `mypy src` clean, `ruff check` clean

### M1 — The coverage model (`model/coverage.py`)
- [ ] `supply_f`, `eff_f`, `demand_f`, `n_f`, `C(p | S, w, D)` exactly as §1.3
- [ ] Analytic unit tests, not smoke tests:
  - [ ] identical `p` and singleton `S={p}` with sufficient budget -> coverage 1.0
  - [ ] disjoint family support -> coverage ~ 0
  - [ ] **adding a dilutant venom lowers coverage** (the R4 mechanism)
  - [ ] dose-response monotone decreasing in `D`
  - [ ] saturation: doubling `B` past sufficiency does not change coverage
  - [ ] invariance: coverage independent of units of `B` and `D` at fixed `B/D`
- [ ] Gate: all six analytic properties pass

### M2 — Curated composition corpus
- [ ] `data/curated/compositions.csv` — every row carries DOI + table + accessed date + method
- [ ] >= 25 Indian populations across >= 6 species
- [ ] Schema validation: keys subset of FAMILIES, each row sums to 1 +/- 1e-6
- [ ] `data/SOURCES.md` with licence + access date + exact query/DOI per source
- [ ] Gate: `venomgap validate-data` passes; row count and species count asserted in tests

### M3 — UniProt/ToxProt ingest + cross-reactivity
- [ ] `ingest/uniprot.py` — REST fetch, on-disk cache, family mapping from `protein_families`
- [ ] `model/crossreact.py` — pairwise identity -> `x_f(q->p)` with fitted floor `theta_f`
- [ ] Gate: within-species cross-reactivity > between-genus; cache makes run reproducible offline

### M4 — Calibration
- [ ] `data/curated/antivenomics_calibration.csv` (non-holdout) and `holdout_targets.csv`
- [ ] `model/calibrate.py` fits `theta_f`, `kappa_f`; **raises on any `holdout=True` row**
- [ ] Gate: a test asserts `HoldoutLeakError` is raised

### M5 — Retrodiction run
- [ ] `experiments/preregistration.md` committed BEFORE fitting, hash recorded
- [ ] `validate/retrodiction.py` evaluates R1-R4 automatically
- [ ] Gate: `results/retrodiction.json` written; outcome reported whatever it is

### M6 — Spatial model
- [ ] `model/spatial.py` — Gaussian-kernel interpolation, `ell` by leave-one-population-out CV
- [ ] Uncertainty grows with distance; districts past cut-off -> `status="unknown"` (grey, never green)
- [ ] Gate: fitted `ell` reported in km; unknown districts identified and counted

### M7 — Optimiser
- [ ] `optimize/objective.py`, `greedy.py`, `local_search.py`, `exact.py`, `weights.py`
- [ ] Coverage-vs-`k` curve with turnover point
- [ ] Optimality gap measured against exact enumeration on a reduced instance
- [ ] Gate: gap reported; non-submodularity stated, no (1-1/e) claim made

### M8 — API + atlas web app
- [ ] FastAPI serving precomputed district results, siting solutions, mixture evaluation
- [ ] Vite + React + TS: Atlas, District detail, Mixture explorer
- [ ] Gate: `k` slider 0..6 redraws; district click shows deficit, limiting families, flags

### M9 — Sensitivity + ablation
- [ ] Leave-one-study-out rank stability (Spearman rho)
- [ ] Cosine-similarity baseline: must fail R4 while the stoichiometric model passes
- [ ] Parameter sweeps over `B`, `D`
- [ ] Gate: `results/sensitivity.json`, `results/ablation.json`

### M10 — Figures, paper skeleton, demo
- [ ] Figures 1-9 from `venomgap figures`, vector PDF + PNG, units on every axis
- [ ] `demo/demo_script.md`, offline fallback bundle
- [ ] Gate: full pipeline reproducible offline from `data/curated/`

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
- [ ] M0-M10 gates passed and evidenced here
- [ ] `results/retrodiction.json` with pre-registered criteria and the honest outcome
- [ ] Figures 1-9 regenerate from one command
- [ ] Demo runs offline
- [ ] README states what the model predicts, what it does not, and which inputs are imputed

---

## Review
Filled in at the end of the build — see bottom of this file.
