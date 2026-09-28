# VenomGap — pre-registration

**Status: committed before any parameter fitting and before any holdout target was inspected
numerically.** The commit that introduces this file contains no fitted parameters, no calibration
data, and no results. The hash of that commit is recorded in `results/retrodiction.json` under
`preregistration_commit`.

Author: Rahul Gunwani
Date filed: 2026-09-29

---

## 1. Hypothesis

India's polyvalent antivenom is an estimator trained on a biased sample. Its immunogen set is
venom pooled from four species collected from essentially one locality in Tamil Nadu. Venom
composition varies geographically. Therefore antivenom efficacy across India is a covariate-shift
problem, and the resulting *deficit* is predictable from toxin-family composition, sequence-derived
cross-reactivity, and a finite antibody budget — without any clinical outcome data as input.

We claim two things:

- **H1 (structural).** A coverage model built on a finite antibody budget is non-monotone in the
  number of immunogens. Mean coverage rises then falls as venoms are added to the immunising
  mixture. This follows from the supply equation alone and is not fitted.
- **H2 (predictive).** The *rank ordering* of predicted deficit across Indian regions recovers
  independently published antivenom-failure findings that the model never saw during fitting.

We do **not** claim the model predicts clinical neutralisation percentages, ED50 values, or vial
requirements. It predicts a risk ordering.

---

## 2. Model form (frozen)

Family vocabulary, fixed at 16 families:

```
3FTx, PLA2, SVMP, SVSP, LAAO, CRISP, CTL, KSPI, DIS, NP, VEGF, NUC, PDE, AChE, CVF, other
```

A venom population `p` is a composition vector `a_p` on the simplex over these families.

An antivenom raised on immunogen set `S` with mixture weights `w` (`sum(w) = 1`) has a **finite
total antibody budget** `B`, allocated across families in proportion to the immunising mixture:

```
supply_f(S, w)        = B * sum_{q in S} w_q * a_qf
eff_f(S, w, p)        = B * sum_{q in S} w_q * a_qf * x_f(q -> p)
demand_f(p, D)        = D * a_pf * kappa_f
n_f                   = min(1, eff_f / demand_f)            (n_f := 1 where demand_f == 0)
C(p | S, w, D)        = sum_f sigma_f * a_pf * n_f  /  sum_f sigma_f * a_pf
deficit               = 1 - C
```

Cross-reactivity from sequence identity, with a family-specific recognition floor `theta_f`:

```
x_f(q -> p) = clip( (I_f(q,p) - theta_f) / (1 - theta_f), 0, 1 )
```

where `I_f(q,p)` is the mean pairwise identity between the family-`f` mature toxin sequences of the
two populations, computed by global pairwise alignment (BLOSUM62, gap open -11, gap extend -1,
identity over aligned columns). `x_f(q -> q) = 1` by construction.

**No other functional form will be substituted after unblinding.** If the model fails, it fails.

---

## 3. Parameters and how each is set

| Symbol | Meaning | Status |
|---|---|---|
| `a_p` | family composition per population | **Data.** Transcribed from published proteomes, DOI + table per row. |
| `sigma_f` | clinical severity weight | **Fixed a priori** in `config.py` from the toxin-pathology literature. Never fitted, never adjusted. |
| `kappa_f` | antibody mass per unit venom mass | **Prior fixed a priori** as `M_ab / M_f` (molar-equivalence argument: a 100 kDa F(ab')2 binding a 7 kDa 3FTx costs more antibody mass per unit toxin mass than binding a 50 kDa SVMP). A global scale `kappa_scale` and per-family log-multipliers for the eight families with calibration support are **fitted on the calibration split only**, L2-regularised toward the prior. |
| `theta_f` | recognition floor | **Fitted on the calibration split only**, bounded to [0.40, 0.85], initialised at 0.65 (conformational epitopes of 15-20 residues lose recognition below roughly 60-70% identity). |
| `B`, `D` | antibody budget, delivered venom dose | **Fixed** to published vial protein content and delivered-dose estimates; swept in sensitivity analysis. |
| `ell` | spatial kernel length scale | **Fitted by leave-one-population-out CV on composition vectors only.** No coverage or outcome information enters this fit. |

`sigma_f` being frozen before any fitting is the load-bearing commitment of this pre-registration.
It is the parameter most open to an accusation of tuning.

---

## 4. Split

**Calibration split.** Published antivenomics / immunorecognition datasets for Indian polyvalent
antivenom, excluding every study underlying R1-R4 and excluding every population labelled
`holdout=True`.

**Holdout split.** The four findings below, plus every population whose deficit ranking those
findings concern.

**Enforcement is in code, not convention.** `venomgap.model.calibrate.fit_parameters` inspects the
`holdout` column of its input and raises `HoldoutLeakError` if any row is `True`. A test asserts
that it raises. Leave-one-study-out sensitivity re-runs the whole fit per excluded study.

---

## 5. Held-out findings and pass criteria

Criteria are stated as rank or sign tests, decided now, and evaluated by
`venomgap.validate.retrodiction` with no manual intervention.

### R1 — *Echis carinatus sochureki*, Rajasthan / north-west India
Source: PMC12845857. Indian polyvalent antivenom poorly recognises *E. c. sochureki* venom,
especially SVMP classes II/III, with poor neutralisation of procoagulant activity.

**Passes if both hold:**
- R1a: the *E. c. sochureki* north-west population ranks in the **top decile of predicted deficit**
  among all curated Indian populations.
- R1b: **SVMP is the single worst limiting family** for that population, i.e. SVMP has the lowest
  `n_f` among families with `sigma_f * a_pf >= 0.02`.

### R2 — *Naja sagittifera*, Andaman and Nicobar Islands
Source: PMC8573199. Polyvalent antivenom fails against *N. sagittifera*; the species is
phylogenetically closer to *N. kaouthia* than to mainland *N. naja*.

**Passes if both hold:**
- R2a: the *N. sagittifera* population ranks in the **top decile of predicted deficit**.
- R2b: substituting *N. kaouthia* for *N. naja* in the immunogen set, holding everything else
  fixed, **increases** predicted coverage of *N. sagittifera* by at least 0.05 absolute.

### R3 — North Indian *Bungarus caeruleus*
Source: PMC6894822. Antivenoms fail against several neglected species and against the North Indian
population of *B. caeruleus* — a Big Four species.

**Passes if both hold:**
- R3a: the North Indian *B. caeruleus* population's predicted deficit is in the **top tertile** of
  all curated Indian populations, despite the species being in the immunogen set.
- R3b: its predicted deficit **exceeds** that of the South Indian (immunogen-source) *B. caeruleus*
  population by at least 0.10 absolute.

### R4 — Immunogen dilution is counterproductive past an optimum
Source: PLOS NTD 2021, doi:10.1371/journal.pntd.0009659.

**Passes if all three hold:**
- R4a: mean burden-weighted national coverage as a function of mixture size `|S|`, each `|S|`
  evaluated at its own **best** subset under a fixed solver, is **non-monotone** — it has an
  interior maximum at some `1 < |S|* < |S|_max`.
- R4b: the fall from the maximum to `|S|_max` is at least 0.02 absolute coverage, i.e. the
  turnover is not numerical noise.
- R4c: the cosine-similarity baseline, run on identical inputs, is **monotone non-decreasing** in
  `|S|` and therefore cannot produce R4a. This is the ablation that shows the mechanism is the
  finite budget, not the data.

### Aggregate criterion
The study is reported as **`k`/4 passed** with every sub-criterion listed individually. A miss is
reported with its diagnosis. **Three of four with an explained miss is reported as three of four.**
No criterion above will be softened, re-scoped or removed after unblinding. If a criterion is found
to be ill-defined on inspection, it is marked `ill_posed` and counted as a **failure**, and the
reason is recorded.

---

## 6. Pre-specified negative controls

- **Permutation control.** Shuffling the composition vectors across populations must destroy the
  R1-R3 rank results (expected pass rate at chance).
- **Cosine ablation.** The cosine baseline must fail R4 (this is criterion R4c) and is expected to
  perform worse than the stoichiometric model on R1-R3 rank recovery.
- **Sequence-free control.** Setting `x_f = 1` for all pairs must degrade R1-R3, showing
  cross-reactivity carries signal rather than being decoration.

---

## 7. What would falsify the framework

- Coverage is monotone in `|S|` under every reasonable `B` (falsifies H1 and the mechanism).
- R1-R3 targets land in the middle of the deficit distribution, or the permutation control performs
  as well as the real model (falsifies H2).
- Fitted `theta_f` hits a bound for most families, indicating the identity signal carries no
  information and the fit is absorbing it into an edge case.
- Leave-one-study-out Spearman rho on district deficit ranks falls below 0.5, indicating the ranking
  is an artefact of a single study.

---

## 8. Reporting commitment

`results/retrodiction.json` will contain, for every criterion: the criterion as stated here, the
computed value, the threshold, and the boolean outcome — plus the preregistration commit hash, the
fitted parameters, and the data version. Every number quoted in the README, the figures and any
slide is read from `results/*.json`. Nothing is retyped by hand.
