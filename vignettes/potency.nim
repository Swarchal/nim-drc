## # Several curves and relative potency
##
## Fit several dose-response curves at once, test whether they are parallel,
## and estimate how potent each test sample is relative to a reference
## standard, with delta-method, profile-likelihood and bootstrap intervals.
##
## ## The data
##
## A simulated potency assay: a reference standard and two test samples,
## each at nine half-log concentrations (0.01 to 100 ng/ml) in duplicate.
## The samples were simulated as 4× and 1/5× as potent as the reference,
## with the same plateaus and slope.
##
## A joint fit takes all the observations in one set of vectors, with
## `curve[i]` saying which curve observation i belongs to (0, 1, 2, …).

import std/[math, sequtils, strformat]
import drc
import drcplot

let response = [
  0.037, 0.07, 0.059, 0.157, 0.207, 0.161, 0.498, 0.508, 1.048, 1.013, 1.632, 1.587,
  1.832, 1.768, 2.002, 1.956, 1.987, 2.026, 0.138, 0.129, 0.269, 0.263, 0.574, 0.492,
  1.175, 1.148, 1.645, 1.593, 1.868, 1.904, 2.021, 1.963, 2.006, 1.988, 1.942, 1.981,
  0.036, 0.05, 0.101, 0.088, 0.069, 0.066, 0.167, 0.162, 0.306, 0.305, 0.799, 0.815,
  1.375, 1.415, 1.789, 1.749, 1.944, 1.885,
]
var doses: seq[float]
var curve: seq[int] # which curve each observation belongs to
for c in 0 .. 2: # 0: reference, 1: sample A, 2: sample B
  for k in 0 .. 8:
    let dose = pow(10.0, -2.0 + 0.5 * float(k))
    doses.add [dose, dose]
    curve.add [c, c]

let samples = ["reference", "sample A", "sample B"]

## ## Separate curves
##
## With no `shared` parameters, `fitLL4Multi` fits each curve its own four
## parameters. That is the same as three separate `fitLL4` fits, except
## that σ is pooled over all three curves, so every curve's standard errors
## rest on 42 residual degrees of freedom rather than 14.
##
## `fit.curves[c]` is curve c as an ordinary `DrcFit`, so everything from
## the [getting started vignette](getting_started.md) works on it:

let separate = fitLL4Multi(doses, response, curve)
for c, sample in samples:
  let f = separate.curves[c]
  echo &"{sample:<10} Bottom {f.params[Bottom]:.3f}  Top {f.params[Top]:.3f}  " &
    &"Hill {f.params[Hill]:.3f}  EC50 {f.ec50:.3f} ng/ml"

## The plateaus and slopes look much alike, and only the EC50s differ by
## much. Plotted, with a helper that draws any joint fit:

let axis = doseAxis("Concentration (ng/ml)", 0.006, 160.0)
let color = seriesColors(samples, [blue, orange, aqua], legendAt = "top-left")
let sampleOfObservation = curve.mapIt(samples[it])

proc curvesFigure(fit: MultiFit, title, name, alt: string) =
  var curves: seq[LabelledFit]
  for c, sample in samples:
    curves.add labelled(sample, fit.curves[c])
  layer(
    axis.curveLines(curves, "Response", color),
    axis.seriesDots(
      axis.observationFrame(doses, response, labels = sampleOfObservation), color
    ),
  )
    .finish(title)
    .figure(name, alt)

separate.curvesFigure(
  "Separate 4PL curves", "separate", "three 4PL curves fitted separately"
)

## ## Parallel curves
##
## A relative potency is only well defined if the curves are *parallel*:
## the same shape, shifted along the log-dose axis. Then one sample behaves
## exactly like a dilution of the other, and the shift is the potency. In
## a joint fit that means sharing Bottom, Top and Hill, so only LogEC50
## differs between curves:

let parallel = fitLL4Multi(doses, response, curve, shared = {Bottom, Top, Hill})
echo parallel.summary

## Shared parameters appear once, the others once per curve (`LogEC50[c]`).
## Is the simpler model good enough? `fTest` compares two nested joint fits
## to the same data, as R's `anova`:

let parallelTest = fTest(parallel, separate)
echo &"parallel vs separate: F = {parallelTest.f:.3f} on {parallelTest.df1} " &
  &"and {parallelTest.df2} df, p = {parallelTest.p:.3f}"

## The six extra parameters of the separate fit (two plateaus and a slope
## for each sample) do not improve it by more than chance would, so there
## is no evidence against parallelism. The test can also go the other way:
## a single curve for all the data (every curve id 0) against separate
## curves asks whether the samples differ at all:

let pooled = fitLL4Multi(doses, response, newSeq[int](doses.len))
let pooledTest = fTest(pooled, separate)
echo &"one curve vs separate: F = {pooledTest.f:.1f} on {pooledTest.df1} " &
  &"and {pooledTest.df2} df, p = {pooledTest.p:.3g}"

## They clearly do. Note that an F test that does not reject parallelism is
## not proof of it: with few data it may simply lack power. `lackOfFit`
## checks a different thing, whether the curves describe the dose means at
## all, using the scatter between duplicates:

let lof = parallel.lackOfFit
echo &"lack of fit: F = {lof.f:.2f} on {lof.df1} and {lof.df2} df, p = {lof.p:.3f}"

## Nor is there evidence of lack of fit.
##
## The parallel curves, drawn the same way:

parallel.curvesFigure(
  "Parallel 4PL curves (shared Bottom, Top and Hill)", "parallel",
  "three parallel 4PL curves fitted jointly",
)

## ## Relative potency
##
## `relativePotency(a, b)` is ECx of curve a divided by ECx of curve b. A
## sample that needs a quarter of the reference's dose for the same effect
## is 4× as potent, so the potency of sample c relative to the reference
## is `relativePotency(0, c)`. Its interval comes from the delta method on
## log10 of the ratio, using the joint covariance, so it is asymmetric.
##
## Two other intervals are available: `relativePotencyProfileCI`
## profiles the likelihood, and `bootstrap` refits resampled data (2000
## times by default) so `bootCI` can take percentiles of any statistic:

let boot = parallel.bootstrap()
echo &"bootstrap: {boot.replicates.len} refits, {boot.failed} failed"

for c in 1 .. 2:
  let rp = parallel.relativePotency(0, c)
  let profile = parallel.relativePotencyProfileCI(0, c)
  let percentile = boot.bootCI(
    proc(f: MultiFit): float =
      f.relativePotency(0, c).ratio
  )
  echo &"{samples[c]}: potency {rp.ratio:.3f}"
  for (methodName, ci) in [
    ("delta", rp.ci), ("profile", profile), ("bootstrap", percentile)
  ]:
    echo &"  {methodName:<10} [{ci.lower:.3f}, {ci.upper:.3f}]"

## On data this clean the three agree closely, and all of them contain the
## simulated potencies of 4 and 0.2. They part company on sparse or partial
## curves: see the [confidence intervals vignette](intervals.md).
##
## `bootValues` gives the statistic for every refit, which shows what the
## bootstrap interval is made of. The dashed lines are its ends, the 2.5th
## and 97.5th percentiles of the refits:

proc potencyOfA(f: MultiFit): float =
  f.relativePotency(0, 1).ratio

let ends = boot.bootCI(potencyOfA)
layer(
  chart(dataFrame(col("potency", boot.bootValues(potencyOfA))))
    .markBar(color = blue, opacity = 0.7)
    .encode(
      x = "potency:Q".bin(maxbins = 40).title(
          "Potency of sample A relative to the reference"
        ),
      y = "count():Q".title("Refits"),
    ),
  chart(dataFrame(col("potency", @[ends.lower, ends.upper])))
    .markRule(color = ink, strokeDash = @[4.0, 4.0])
    .encode(x = "potency:Q"),
)
  .finish("Bootstrap distribution of the potency of sample A", h = 220)
  .figure(
    "bootstrap", "histogram of the potency of sample A over 2000 bootstrap refits"
  )

## ## Potency at other response levels
##
## `relativePotency` takes `pct`, to compare EC20s or EC80s instead of
## EC50s. For parallel curves the answer does not depend on it: shifting a
## curve along the log-dose axis moves every ECx by the same amount, so
## every ratio is the same. For the separate curves, whose slopes and
## plateaus differ slightly, the ratio drifts with the response level:

for pct in [20.0, 50.0, 80.0]:
  echo &"EC{pct:g}: parallel {parallel.relativePotency(0, 1, pct).ratio:.3f}, " &
    &"separate {separate.relativePotency(0, 1, pct).ratio:.3f}"

## That drift is what non-parallel curves look like. Here it is no more than
## noise (the F test above found no evidence against parallelism), so the
## parallel model's single potency stands. When the drift is large and the
## F test rejects, a single potency figure does not describe the samples.
##
## ## The joint estimates
##
## `param(c, p)` is curve c's estimate of p, and `cov(a, pa, b, pb)` a
## covariance from the joint fit, including between curves. Shared
## parameters have one value for all curves:

echo &"Hill: {parallel.param(0, Hill):.4f} (reference), " &
  &"{parallel.param(2, Hill):.4f} (sample B)"
let corr =
  parallel.cov(0, LogEC50, 1, LogEC50) /
  sqrt(parallel.cov(0, LogEC50, 0, LogEC50) * parallel.cov(1, LogEC50, 1, LogEC50))
echo &"correlation of LogEC50[0] and LogEC50[1]: {corr:.3f}"

## The LogEC50s of different curves are correlated because they lean on
## the same shared plateaus and slope; the potency interval accounts for
## that. Each `parallel.curves[c]` carries curve c's block of the joint
## covariance and the joint df and σ, so single-curve tools such as
## `ecxCI` and `confidenceBand` work on it unchanged:

for c, sample in samples:
  let ci = parallel.curves[c].ec50CI
  echo &"{sample:<10} EC50 {parallel.curves[c].ec50:.3f} ng/ml  " &
    &"[{ci.lower:.3f}, {ci.upper:.3f}]"

## ## Next
##
## - [Confidence intervals](intervals.md): Wald, profile-likelihood and
##   bootstrap intervals for a single curve, and when they disagree.
## - [Choosing a model](model_choice.md): 3PL, 4PL or 5PL. `fitLL3Multi`
##   and `fitLL5Multi` fit those models jointly.
## - [Standard curves](standard_curves.md): reading unknown doses off a
##   fitted curve.
