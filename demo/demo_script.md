# VenomGap — 3 minute demo

**Runs entirely offline.** Both servers read committed data; nothing reaches the network.
Rehearse with the Wi-Fi off — that is the point of the fallback, not a precaution.

## Before you start

```bash
export VENOMGAP_OFFLINE=1
python -m uvicorn venomgap.api.app:app --port 8000 &
npm --prefix web run dev
```

Check `http://localhost:8000/health` returns `{"ok": true, "missing_results": []}`.
Open `http://localhost:5173` on the **Atlas** tab with the `k` slider at **0**.
Have `figures/fig03_dilution_curve.pdf` open in a second window as the backup for beat 2.

---

## 0:00 – 0:25 · The problem

> "India's only antivenom is raised on venom from four species, collected almost entirely from one
> facility in Tamil Nadu. Venom composition varies across the country. So the antivenom is an
> **estimator trained on a biased sample**, and that is a measurable problem, not a rhetorical one."

**Show:** the Atlas as loaded. Point at the red ring in the south — the collection locality — and at
the dark north-east and islands.

> "Coverage is 66% nationally. It is highest around the one place the venom comes from and falls
> away from it. The grey district is one where the nearest published proteome is too far away to
> extrapolate from. **Grey means unknown. Unknown is a finding, not a blank.**"

---

## 0:25 – 1:10 · The mechanism  ← the important beat

**Do:** switch to **Mixture explorer**. It opens with four venoms, one per Big Four species.

> "This is the immunising mixture. Let me add venoms to it."

**Do:** tick *Bungarus fasciatus* — West Bengal. Wait for the number to settle.
**Do:** tick *Bungarus sindanus* — Rajasthan.

**Show:** coverage falls from **60.0% to 51.2%**, the bars turn orange, and the callout appears.

> "Coverage went **down**. The model was never told that adding venoms can hurt. It falls out of one
> equation: a vial holds a fixed mass of antibody, so antibody raised against one more venom is
> antibody *not* raised against the toxins a given bite actually delivers.
>
> That is exactly what an experimental study reported in 2021 — and that study was **held out**. The
> model never saw it."

> "A cosine-similarity metric — which is what the existing literature uses — cannot do this. It has
> no budget, so it rises monotonically no matter how many venoms you add. We ran it on identical
> inputs: it goes from 0.73 to 0.95 and never turns over."

*(Backup: `figures/fig03_dilution_curve.pdf` shows both curves on one axis.)*

---

## 1:10 – 2:00 · The blinded retrodiction

**Do:** switch to the **Retrodiction** tab.

> "We hid four published antivenom-failure findings before fitting anything. The calibration code
> doesn't just avoid them — it **raises an exception** if a held-out row reaches a parameter fit,
> and there are tests that assert it raises."

> "We passed **one of four**. I want to be straight about that, and about why it's more interesting
> than four of four would have been."

**Point at R4:** all three sub-criteria pass.

> "R4 is the structural one and it passed completely."

**Point at the diagnostics block at the bottom.**

> "The other three missed for one identifiable reason, and it's in the results file as a computed
> number, not a story. The recognition threshold sits above **92.7%** of every between-species
> sequence identity we measured, so cross-reactivity is forced to zero almost everywhere and the
> model has almost no paraspecific coverage. That threshold was pre-registered to be *fitted* — and
> the calibration data turned out to be unable to constrain it at all, so it kept its prior."

**Point at R2b.**

> "Here's the sharpest version. The sequence data **does** rank the monocled cobra closer to the
> Andaman cobra than the mainland cobra is — 0.388 against 0.278 — which is precisely what the
> held-out paper claims on phylogenetic grounds, and we got there independently. But the
> recognition floor throws that ordering away before it can act, so the substitution moves coverage
> by exactly zero. The signal is in the data. The model's threshold discards it. **That is a
> specific, fixable diagnosis, not a shrug.**"

> "Only one parameter was ever fitted: a single global scale. The severity weights, the
> stoichiometry, every recognition threshold — all fixed before we looked."

---

## 2:00 – 2:40 · The answer

**Do:** back to **Atlas**. Drag `k` from 0 to 3, slowly.

**Show:** pins drop, coverage climbs to **74.4%**, "+8.2 points", "23 districts out of high deficit".

> "Given a budget for *k* new collection sites, this solves for where they go."

**Point at the coverage-vs-k chart, the blue line.**

> "But look at the blue line — the vial split evenly across its venoms. It peaks at **one** added
> site, then **dips**. That's the dilution effect again, now in the siting problem."

**Point at the green line.**

> "The green line re-optimises the mixture weights, and the dip disappears. So the operational
> finding is sharper than 'collect more venom': **adding a venom to a polyvalent antivenom without
> re-balancing the mixture can make it worse. Re-balancing is what recovers the gain.**"

> "And we measured the solver rather than trusting it. The objective is non-monotone and not
> submodular, so there is no `1 − 1/e` guarantee to lean on — we enumerated 6,195 subsets exactly on
> a reduced instance. The heuristic matched it. Gap zero."

---

## 2:40 – 3:00 · Close

> "There is one more result, and it is the one I'd act on first."

> "We fitted how far a venom proteome generalises, by cross-validation, from composition alone. The
> answer is **1600 km with a plateau running past the length of the country**. That means distance
> carries almost no predictive power at the sampling density that exists today. You cannot infer
> regional venom composition from what has been published. You have to go and measure it."

> "India has eight antivenom manufacturers and essentially one venom source. This says which regions
> to sample first, how many sites is enough, and — measured, not asserted — how little we currently
> know."

---

## If something breaks

| Symptom | Do this |
|---|---|
| API 500s on load | The model builds lazily on first request. Wait 10 s and reload. |
| Coverage doesn't move in the explorer | Requests are superseded by ticket; click one venom and pause. |
| Anything else | Open `figures/*.pdf`. Every figure is generated from `results/*.json` and tells the same story without a server. |

## Rehearsal checklist

- [ ] `VENOMGAP_OFFLINE=1` exported in both shells
- [ ] Wi-Fi **off**, full run-through completes
- [ ] `python -m venomgap.cli pipeline` green from a clean checkout
- [ ] `figures/` regenerated and open in a second window
- [ ] Screen capture recorded as the final fallback
