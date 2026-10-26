## # Fitting many curves: a compound screen
##
## Fit one curve per compound across a whole screen, decide which compounds
## are real hits from the fits' warnings and intervals, and rank them. The
## data are simulated, so the estimates can be checked against the truth.
##
## ## A simulated screen
##
## 32 compounds, each tested at eight concentrations in a 3-fold dilution
## series from 30 µM down to 14 nM, in duplicate. Responses are % of the
## untreated control, simulated from a 4PL with Bottom 0 and Top 100,
##
## ```
## response = 100 / (1 + (dose / IC50)^h) + noise,   noise ~ Normal(0, 6²)
## ```
##
## in three kinds of compound: 16 actives with an IC50 between 30 nM and
## 3 µM, 8 weak ones with an IC50 between 30 and 300 µM (beyond the top
## dose, so only the start of the curve is seen), and 8 inactives whose
## response is flat at 100%. Within the screen their order is shuffled.

import std/[algorithm, math, random, sequtils, strformat]
import drc
import drcplot

type Compound = object
  name: string
  kind: string ## "active", "weak" or "inactive": the truth, for checking
  trueIc50: float ## molar; NaN for an inactive
  doses, responses: seq[float]

var rng = initRand(2024)
var kinds = newSeqWith(16, "active") & newSeqWith(8, "weak") & newSeqWith(8, "inactive")
rng.shuffle(kinds)

var screen: seq[Compound]
for i, kind in kinds:
  var c = Compound(name: &"C{i + 1:02}", kind: kind, trueIc50: NaN)
  let hill = rng.rand(0.8 .. 1.6)
  case kind
  of "active":
    c.trueIc50 = pow(10.0, rng.rand(-7.5 .. -5.5))
  of "weak":
    c.trueIc50 = pow(10.0, rng.rand(-4.5 .. -3.5))
  for k in 0 .. 7:
    let dose = 30e-6 / pow(3.0, float(k))
    for rep in 0 .. 1:
      let mean =
        if kind == "inactive":
          100.0
        else:
          100 / (1 + pow(dose / c.trueIc50, hill))
      c.doses.add dose
      c.responses.add mean + rng.gauss(0.0, 6.0)
  screen.add c

echo &"{screen.len} compounds, {screen[0].doses.len} wells each"

## ## Fitting every compound
##
## The responses are normalised to the control, so the upper plateau is
## known: Top is fixed at 100. The lower plateau is not (a partial
## inhibitor stops above 0), so Bottom is fitted, but bounded to 0–100: a
## noisy or flat curve could otherwise put it anywhere. With Top fixed this
## is a 4PL with three free parameters.
##
## A fit never raises because the data are bad (only for invalid input,
## such as a negative dose), so there is no need for `try` in the loop.
## Problems are reported in `fit.warnings` instead, and the triage below
## uses them.

var fits: seq[DrcFit]
for c in screen:
  fits.add fitLL4(
    c.doses, c.responses, fixed = [(Top, 100.0)], bounds = [(Bottom, 0.0 .. 100.0)]
  )

for warning in FitWarning:
  let count = fits.countIt(warning in it.warnings)
  if count > 0:
    echo &"{count:>2} fits: {warning}"
echo &"{fits.countIt(not it.converged)} fits did not converge"

## Every fit converged, but that alone says little. The singular fits are
## flat curves, where nothing locates an EC50. The EC50 warning means the
## curve's midpoint was extrapolated beyond the tested doses, as for half of
## the weak compounds. The bound warnings are curves whose Bottom ended at
## 0, as expected for full inhibitors, since the noise can pull a free
## Bottom below 0.
##
## Warnings are not enough on their own. A flat curve with a slight dip
## can fit a small step (Bottom near 100) with an EC50 inside the tested
## range and no warning at all: its EC50 is meaningless, and its interval is
## (0, ∞).
##
## For a screen of thousands of compounds, `FitOptions(multiStart: false)`
## fits from a single start, which is about 13 times faster. The
## multi-start guards against local minima on sparse or noisy curves, so
## keep it for the hits you go on to report.
##
## ## Triage
##
## So a compound is called a hit only when its fit says where its IC50 is
## (the EC50 lies within the tested doses and its 95% interval spans less
## than a factor of 10) *and* it has a real effect (Bottom, the maximal
## response, at most 50%). Everything else is set aside, as "weak" if the
## fitted curve has fallen by at least 20% at the top dose, and as
## "inactive" if not.

type Call = enum
  Hit = "hit"
  Weak = "weak"
  Inactive = "inactive"

proc call(fit: DrcFit, topDose: float): Call =
  let ci = fit.ec50CI
  if fwEc50OutsideRange notin fit.warnings and fit.converged and ci.upper / ci.lower < 10 and
      fit.params[Bottom] <= 50:
    Hit
  elif fit.predict(topDose) <= 80:
    Weak
  else:
    Inactive

var calls: seq[Call]
for i, fit in fits:
  calls.add fit.call(screen[i].doses.max)

echo "truth      called hit   weak   inactive"
for kind in ["active", "weak", "inactive"]:
  var counts: array[Call, int]
  for i, c in screen:
    if c.kind == kind:
      inc counts[calls[i]]
  echo &"{kind:<10} {counts[Hit]:>10} {counts[Weak]:>6} {counts[Inactive]:>10}"

## Every active was found, with no false hits. Two of the weak compounds
## fell by less than 20% within the tested doses, so they look inactive:
## a screen can only see what its dose range reaches.
##
## ## The curves
##
## Eight of the compounds (four actives, two weak, two inactive), one
## panel each, titled with the call the triage made. The data for small
## multiples go in one long frame, with a `compound`
## column to split the panels on and a `kind` column that says which rows
## are the fitted curve and which the observations.

let axis = doseAxis("Dose (µM)", 0.01, 40.0, fittedToDisplay = 1e6)
var shown: seq[int] # the first four actives, two weak and two inactive
for (kind, count) in [("active", 4), ("weak", 2), ("inactive", 2)]:
  shown.add toSeq(0 ..< screen.len).filterIt(screen[it].kind == kind)[0 ..< count]
shown.sort

var panel, rowKind: seq[string]
var dose, response: seq[float]
for i in shown:
  let label = &"{screen[i].name}: {calls[i]}"
  for d in axis.doseGrid(100):
    panel.add label
    rowKind.add "fit"
    dose.add axis.toDisplay(d)
    response.add fits[i].predict(d)
  for j, d in screen[i].doses:
    panel.add label
    rowKind.add "observed"
    dose.add axis.toDisplay(d)
    response.add screen[i].responses[j]
let frame = dataFrame(
  col("compound", panel),
  col("kind", rowKind),
  col("dose", dose),
  col("response", response),
)
let yAxis = "response:Q".title("Response (% of control)").scale(
    domain = vals(-20.0, 130.0), nice = false
  )
# A facet splits the data of the chart inside it, so the frame goes on the
# layer, and the two marks (`Chart()`, with no data) take their rows from it.
layer(
  Chart()
    .transformFilter("datum.kind == 'fit'")
    .markLine(color = blue, strokeWidth = 2)
    .encode(x = axis.x, y = yAxis),
  Chart()
    .transformFilter("datum.kind == 'observed'")
    .markPoint(filled = true, color = ink, size = 20)
    .encode(x = axis.x, y = yAxis),
)
  .withData(inlineData(frame))
  .properties(width = 180, height = 130)
  .facet(facet = "compound:N".title(""), columns = 4)
  .configureBackground("#ffffff")
  .configureAxis(gridColor = "#eeeeee")
  .figure("curves", "Fitted curves for eight compounds of the screen, one per panel")

## The hits' curves are well determined. The weak compounds have only begun
## to fall at the top doses, and the inactives are flat, so the EC50s
## fitted to them mean nothing.
##
## ## Ranking the hits
##
## Hits ranked by potency, with the true IC50 alongside (which a real
## screen would not have). Hill is negative because the response falls
## with dose:

proc nm(dose: float): string =
  ## A molar dose in whole nM.
  $int(round(dose * 1e9))

var hits = toSeq(0 ..< screen.len).filterIt(calls[it] == Hit)
hits.sort(
  proc(a, b: int): int =
    cmp(fits[a].ec50, fits[b].ec50)
)
echo "compound   IC50 (nM)   95% CI (nM)      Hill   Bottom   true IC50 (nM)"
for i in hits:
  let fit = fits[i]
  let ci = fit.ec50CI
  let interval = nm(ci.lower) & "–" & nm(ci.upper)
  echo &"{screen[i].name:<10} {nm(fit.ec50):>9}   {interval:<15}  " &
    &"{fit.params[Hill]:>5.2f}  {fit.params[Bottom]:>6.1f}   {nm(screen[i].trueIc50):>14}"

## The estimated IC50 against the true one, on log scales, with each hit's
## 95% interval. Points on the dashed diagonal are exactly right:

var trueNm, estNm, lowerNm, upperNm: seq[float]
for i in hits:
  let ci = fits[i].ec50CI
  trueNm.add screen[i].trueIc50 * 1e9
  estNm.add fits[i].ec50 * 1e9
  lowerNm.add ci.lower * 1e9
  upperNm.add ci.upper * 1e9
let logAxis = proc(field, title: string): Encoding =
  field.scale(scaleType = scLog, domain = vals(50.0, 5000.0), nice = false).axis(
    title = title, values = vals(100.0, 1000.0), format = "~g"
  )
let estimates = dataFrame(
  col("true", trueNm),
  col("estimate", estNm),
  col("lower", lowerNm),
  col("upper", upperNm),
)
layer(
  chart(dataFrame(col("true", @[50.0, 5000.0]), col("estimate", @[50.0, 5000.0])))
    .markLine(color = muted, strokeDash = @[4.0, 4.0])
    .encode(x = logAxis("true:Q", "True IC50 (nM)"), y = "estimate:Q"),
  chart(estimates).markRule(color = blue, strokeWidth = 1.5).encode(
    x = "true:Q", y = logAxis("lower:Q", "Estimated IC50 (nM)"), y2 = "upper:Q"
  ),
  chart(estimates).markPoint(filled = true, color = blue, size = dotSize).encode(
    x = "true:Q", y = "estimate:Q"
  ),
)
  .finish("Estimated vs true IC50 of the hits, with 95% intervals")
  .figure("ic50", "Estimated against true IC50 for each hit, with 95% intervals")

var covered = 0
for k in 0 ..< hits.len:
  if lowerNm[k] <= trueNm[k] and trueNm[k] <= upperNm[k]:
    inc covered
echo &"{covered} of {hits.len} intervals contain the true IC50"

## That is fewer than the 15 of 16 that 95% intervals would contain on
## average, though not far outside chance for 16 intervals. These are
## Wald intervals, from the curvature at the estimate; the [confidence
## intervals vignette](intervals.md) compares them with profile-likelihood
## and bootstrap intervals.
##
## ## Next
##
## - [Getting started](getting_started.md): one curve in detail, with
##   residual checks.
## - [Choosing a model](model_choice.md): fixed plateaus and bounds, and
##   when a 3PL or 5PL suits better.
## - [Confidence intervals](intervals.md): the wide or open-ended intervals
##   of partial curves like the weak compounds here.
## - [Outliers](outliers.md): robust fits for screens with bad wells.
