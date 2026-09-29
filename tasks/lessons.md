# Lessons

Format: `[YYYY-MM-DD] — [what went wrong] → [rule going forward]`

- [2026-09-29] — The first cosine ablation used the weighted-mean form, which declines with mixture
  size for a purely arithmetic reason (it averages fixed per-venom scores). It looked like the
  baseline reproduced the dilution effect, when in fact it was failing the criterion for the wrong
  reason. → **Always implement an ablation in its strongest form.** A straw man that fails proves
  nothing; the comparison is only worth making against the version a critic would defend.

- [2026-09-29] — `limiting_families` tie-broke alphabetically. Where cross-reactivity is zero many
  families sit at exactly `n_f = 0`, so the "worst limiting family" was decided by the alphabet and
  R1b reported CTL instead of SVMP, which contributed 76% of the deficit. → **When sorting on a
  quantity that ties often, decide the tie-break deliberately.** A tie-break is part of the
  definition, not an implementation detail.

- [2026-09-29] — The siting optimiser selected `k` populations as the *entire* immunogen, so `k=1`
  scored below `k=0` and the coverage-versus-k curve was meaningless. The bug was invisible until
  the curve was plotted in the web app. → **Plot the headline curve early.** A number in a JSON file
  can be wrong in a way a chart makes obvious in one second.

- [2026-09-29] — `theta_f` was pre-registered as "fitted on the calibration split" but turned out to
  have exactly zero data sensitivity, because the calibration targets share species with the
  immunogen set where cross-reactivity is 1 by construction. Left unexamined it would have looked
  like a tuned parameter sitting at a convenient value. → **Measure identifiability, don't assume
  it.** Probe the data term alone, not the objective including its regulariser — the first probe
  returned the L2 penalty and reported every parameter as identified.

- [2026-09-29] — The leave-one-population-out CV for the spatial length scale picked the largest
  value on the grid, which looked like a real optimum until the grid was extended. → **Extend the
  grid past the plausible range before reporting an optimum.** An edge solution and a plateau are
  indistinguishable from inside a short grid, and a flat-bottomed curve must be reported as a
  plateau, not as a point estimate.

- [2026-09-29] — An automated line-wrapper applied to long lines broke a multi-line f-string and a
  function signature. → **Never bulk-reformat source with a regex over line lengths.** Wrap by hand
  or use a real formatter that parses the code.
