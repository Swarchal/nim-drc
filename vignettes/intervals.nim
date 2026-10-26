## # Confidence intervals
##
## drc gives three kinds of interval for a parameter or an effective dose:
## Wald (delta-method) intervals, profile-likelihood intervals and bootstrap
## intervals. On a well-designed experiment they agree. On a partial curve
## they do not, and the differences say how much the data really pin down.
##
## ## Three intervals for one EC50
##
## Start with the ryegrass data from [Getting started](getting_started.md),
## which cover both plateaus with replicates at every dose.

import std/[math, sequtils, strformat]
import drc
import drcplot

## Its doses and responses, `ryeDoses` and `rootLength`, are as there and
## not repeated here.

#| hide
let ryeDoses = [
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.94, 0.94, 0.94, 1.88, 1.88, 1.88, 3.75, 3.75, 3.75,
  7.5, 7.5, 7.5, 15.0, 15.0, 15.0, 30.0, 30.0, 30.0,
]
let rootLength = [
  7.58, 8.0, 8.3285714, 7.25, 7.375, 7.9625, 8.3555556, 6.9142857, 7.75, 6.8714286,
  6.45, 5.9222222, 1.925, 2.8857143, 4.2333333, 1.1875, 0.8571429, 1.0571429, 0.6875,
  0.525, 0.825, 0.25, 0.22, 0.44,
]

## - **Wald** (`confint`, `ecxCI`, `ec50CI`): estimate ± t·SE, from the
##   parameter covariance. For an ECx it is computed on log10 dose and
##   back-transformed. Instant, but it assumes the sum of squares is
##   quadratic around the estimate.
## - **Profile likelihood** (`profileCI`, `ecxProfileCI`, `ec50ProfileCI`):
##   holds the parameter at each candidate value, refits the others, and keeps
##   the values the data do not reject, as R's `confint` on an `nls` fit.
## - **Bootstrap** (`bootstrap`, then `bootCI`): refits many resampled copies
##   of the data and takes percentiles of any statistic over the refits.

let rye = fitLL4(ryeDoses, rootLength)
let ryeBoot = rye.bootstrap() # 2000 refits, seed 1
echo &"{ryeBoot.replicates.len} refits, {ryeBoot.failed} failed"

proc show(name: string, ci: ConfInt) =
  echo &"{name:<12} [{ci.lower:.3f}, {ci.upper:.3f}]"

echo &"EC50 = {rye.ec50:.3f} mM"
show "Wald", rye.ec50CI
show "profile", rye.ec50ProfileCI
show "bootstrap", ryeBoot.bootCI(ec50)

## All three agree to within a few hundredths of a millimolar. When they do,
## the cheap Wald interval is all you need.
##
## ## A partial curve
##
## Now a compound whose response has only fallen to about two thirds by the
## top dose: the data never reach the lower plateau.

let doses = [1e-9, 1e-9, 1e-8, 1e-8, 1e-7, 1e-7, 1e-6, 1e-6, 1e-5, 1e-5] # molar
let response = [100.4, 99.1, 98.7, 97.2, 93.0, 91.2, 80.1, 82.4, 68.0, 66.3]
let fit = fitLL4(doses, response)
let boot = fit.bootstrap()
echo &"{boot.replicates.len} refits, {boot.failed} failed"

proc row(name: string, est: float, w, p, b: ConfInt, scale = 1.0) =
  proc fmt(ci: ConfInt): string =
    &"[{ci.lower * scale:7.2f}, {ci.upper * scale:7.2f}]"

  echo &"{name:<10}{est * scale:>8.2f}   {fmt(w)}   {fmt(p)}   {fmt(b)}"

echo "          estimate   Wald 95%             profile 95%          bootstrap 95%"
for param in [Bottom, Top, Hill]:
  let paramOf = proc(f: DrcFit): float =
    f.params[param]
  row(
    $param,
    fit.params[param],
    fit.confint(param),
    fit.profileCI(param),
    boot.bootCI(paramOf),
  )
row("EC50 (µM)", fit.ec50, fit.ec50CI, fit.ec50ProfileCI, boot.bootCI(ec50), 1e6)

## The Wald intervals look reassuring: Bottom somewhere between 18 and 80,
## an EC50 between 0.2 and 34 µM. The profile intervals disagree. Nothing in
## the data rules out a Bottom far below the observations, and with it an
## EC50 far beyond the top dose, so those ends are infinite. An infinite end
## means the data do not bound the parameter within 10·t standard errors of
## the estimate (R's `confint` reports `NA` there). Top and Hill, which the
## data do constrain, get much the same interval from every method.
##
## The bootstrap shows the same thing another way. Here are 100 of the
## refits, each a curve that fits a resampled copy of the data about as
## well as the original, drawn out to doses far beyond the experiment:

let axis = doseAxis("Dose (µM)", 5e-4, 1e4, fittedToDisplay = 1e6)
let replicateCurves = toSeq(0 ..< 100).mapIt(labelled($it, boot.replicates[it]))
layer(
  chart(axis.curveFrame(replicateCurves))
    .markLine(color = blue, strokeWidth = 1, opacity = 0.2)
    .encode(
      x = axis.x,
      y = "response:Q".title("Response (%)").scale(domain = vals(-40.0, 110.0)),
      detail = "series:N",
    ),
  axis.curveLine(fit, "Response (%)", color = ink),
  axis.dots(axis.observationFrame(doses, response)),
)
  .finish("100 bootstrap refits of a partial curve (fitted curve in grey)")
  .figure(
    "partial-replicates",
    "Bootstrap refits of a partial curve fan out beyond the top dose",
  )

## Where there are data the refits agree. Beyond the top dose they fan out,
## because the data say little about where the curve levels off. Some refits
## failed to converge altogether (`boot.failed`, above); they are left out
## of every bootstrap interval, so on data like these the bootstrap
## understates the uncertainty, if anything.
##
## `bootValues` gives a statistic for every refit, so its distribution can be
## drawn. Here is the log10 EC50 of the refits, with the three 95% intervals
## below it on the same scale. Refits beyond 10⁵ µM are piled into the last
## bin, and an arrow marks an interval that runs off the axis:

const (lo, hi) = (-1.0, 5.0) # the axis, in log10 µM

let logEc50s = boot.bootValues(
  proc(f: DrcFit): float =
    log10(f.ec50 * 1e6)
)
echo "refits with EC50 above 100 µM: ",
  logEc50s.countIt(it > 2.0), " of ", logEc50s.len
let logAxis = "logEC50:Q".scale(domain = vals(lo, hi), nice = false)
let histogram = chart(dataFrame(col("logEC50", logEc50s.mapIt(min(it, hi - 1e-6)))))
  .markBar(color = blue, opacity = 0.6)
  .encode(
    x = logAxis.bin(step = 0.25, extent = @[lo, hi]).axis(title = ""),
    y = "count():Q".title("Refits"),
  )

proc intervalBars(cis: openArray[(string, ConfInt)]): Chart =
  ## One horizontal bar per EC50 interval, on log10 µM, clipped to the axis,
  ## with an arrow where it is clipped.
  var names, clipped: seq[string]
  var lower, upper, arrowAt: seq[float]
  for (name, ci) in cis:
    let (l, u) = (log10(ci.lower * 1e6), log10(ci.upper * 1e6))
    names.add name
    lower.add max(l, lo)
    upper.add min(u, hi)
    if u > hi:
      clipped.add name
      arrowAt.add hi
  let y = "method:N".sort(names.mapIt(v(it))).title("")
  let bars = chart(
      dataFrame(col("method", names), col("logEC50", lower), col("upper", upper))
    )
    .markRule(color = ink, strokeWidth = 3)
    .encode(x = logAxis.axis(title = "log10 EC50 (µM)"), x2 = "upper:Q", y = y)
  let arrows = chart(dataFrame(col("method", clipped), col("logEC50", arrowAt)))
    .markPoint(
      filled = true, color = ink, size = 80, raw = %*{"shape": "triangle-right"}
    )
    .encode(x = logAxis, y = y)
  layer(bars, arrows)

vconcat(
  histogram.properties(width = width, height = 200),
  intervalBars(
    [
      ("Wald", fit.ec50CI),
      ("profile", fit.ec50ProfileCI),
      ("bootstrap", boot.bootCI(ec50)),
    ]
  )
    .properties(width = width, height = 90),
)
  .properties(title = "Bootstrap EC50s of the partial curve, and three 95% intervals")
  .configureBackground("#ffffff")
  .configureView(noStroke = true)
  .configureAxis(gridColor = "#eeeeee")
  .figure(
    "partial-ec50",
    "Bootstrap EC50 distribution with Wald, profile and bootstrap intervals",
  )

## The distribution has a long right tail: about 6% of the refits put the
## EC50 above 100 µM, ten times the top dose. The symmetric Wald
## interval cuts that tail off and reaches too far to the left; the profile
## and bootstrap intervals both follow it.
##
## **What to do about it.** No interval method can recover a plateau the
## experiment never reached. Either extend the dose range, or fix or bound
## the plateau from prior knowledge (see [Choosing a model](model_choice.md)),
## and say so when reporting the EC50.
##
## ## Residuals or cases?
##
## `bootstrap` resamples in one of two ways. `Residuals` (the default) adds
## residuals drawn with replacement to the fitted values: it assumes the
## model and the weights are right. `Cases` draws whole observations with
## replacement, which assumes neither, but a resample can miss a dose
## altogether, so it needs more data:

for (name, f) in [("ryegrass", rye), ("partial", fit)]:
  for resample in [Residuals, Cases]:
    let b = f.bootstrap(resample = resample)
    let ci = b.bootCI(ec50)
    let unit = if name == "ryegrass": 1.0 else: 1e6
    echo &"{name:<9} {$resample:<10} EC50 95% CI [{ci.lower * unit:6.2f}, {ci.upper * unit:8.2f}]   {b.failed:>3} of 2000 failed"

## On the ryegrass data the case interval is somewhat wider, as it often is
## with only three replicates per dose. On the partial curve, with ten
## observations and two per dose, the schemes disagree at both ends: with
## so little data neither interval should be taken at face value.
##
## ## Other confidence levels
##
## Every interval takes a `level` (0.95 by default):

show "Wald 90%", rye.ec50CI(level = 0.90)
show "profile 90%", rye.ec50ProfileCI(level = 0.90)
show "boot 90%", ryeBoot.bootCI(ec50, level = 0.90)

## `bootCI` can reuse one bootstrap for any level and any statistic, since
## the refits are kept in `boot.replicates`.
##
## ## Relative and absolute IC50
##
## `ec50` is the *relative* IC50: the dose half-way between the two fitted
## plateaus. For responses normalised to % of control, the *absolute* IC50,
## where the curve crosses 50%, is often wanted instead. `doseAt` gives it,
## and `doseAtCI` its interval. The two differ when the fitted plateaus are
## not at 0 and 100:

let viabilityDoses = [1e-9, 1e-9, 1e-8, 1e-8, 1e-7, 1e-7, 1e-6, 1e-6, 1e-5, 1e-5]
let viability = [101.3, 97.8, 96.1, 99.0, 78.2, 74.9, 31.5, 36.0, 4.1, 1.2]
let v = fitLL4(viabilityDoses, viability)
let (rel, ab) = (v.ec50CI, v.doseAtCI(50.0))
echo &"plateaus: {v.params[Bottom]:.1f} and {v.params[Top]:.1f}"
echo &"relative IC50 = {v.ec50 * 1e9:.1f} nM (95% CI {rel.lower * 1e9:.1f}–{rel.upper * 1e9:.1f})"
echo &"absolute IC50 = {v.doseAt(50.0) * 1e9:.1f} nM (95% CI {ab.lower * 1e9:.1f}–{ab.upper * 1e9:.1f})"

let vAxis = doseAxis("Dose (nM)", 0.5, 2e4, fittedToDisplay = 1e9)
let midpoint = (v.params[Bottom] + v.params[Top]) / 2

proc responseLabel(response: float, label: string, dy: float): Chart =
  ## `label` at the left edge of the plot, just above or below `response`.
  chart(
    dataFrame(
      col("dose", @[vAxis.lo]), col("response", @[response]), col("label", @[label])
    )
  )
    .markText(align = "left", dx = 4, dy = dy, color = ink, fontSize = 12)
    .encode(x = vAxis.x, y = "response:Q", text = "label:N")

layer(
  responseRule(50.0),
  responseLabel(50.0, "50%", -8),
  responseRule(midpoint),
  responseLabel(midpoint, &"half-way, {midpoint:.1f}%", 12),
  vAxis.doseRule(v.ec50, &"relative {round(v.ec50 * 1e9).int} nM", dy = -110),
  vAxis.doseRule(
    v.doseAt(50.0), &"absolute {round(v.doseAt(50.0) * 1e9).int} nM", dy = 130
  ),
  vAxis.curveLine(v, "Viability (% of control)"),
  vAxis.dots(vAxis.observationFrame(viabilityDoses, viability)),
)
  .finish("Relative IC50 (half-way between plateaus) and absolute IC50 (50%)")
  .figure("relative-absolute", "Relative and absolute IC50 marked on a fitted curve")

## The fitted Bottom is below 0, so the half-way response is a little below
## 50%, and the relative IC50 is at a higher dose than the absolute one. The
## relative IC50's interval does not depend on how well the plateaus are
## known (it is always half-way between them), but the absolute IC50's does.
## Which one to report is a scientific choice, not a statistical one; say
## which it is.
##
## ## Summary
##
## | | Wald | Profile | Bootstrap |
## |---|---|---|---|
## | Functions | `confint`, `ecxCI`, `doseAtCI` | `profileCI`, `ecxProfileCI` | `bootstrap` + `bootCI` |
## | Cost | none | about one fit per interval | 2000 refits, reusable |
## | Shape | symmetric (on log dose for ECx) | follows the data | follows the data |
## | Covers | parameters, ECx, `doseAt` | parameters, ECx | any function of a fit |
##
## Use Wald intervals by default, and check a profile or bootstrap interval
## whenever the curve is incomplete, the data sparse, or the Wald interval
## reaches outside the tested doses.
##
## ## Next
##
## - [Choosing a model](model_choice.md): fixed plateaus and bounds, for
##   when the data cannot place a plateau on their own.
## - [Relative potency](potency.md): the same three kinds of interval for a
##   potency ratio between curves.
## - [Standard curves](standard_curves.md): intervals for doses read off a
##   calibration curve.
