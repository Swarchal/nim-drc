## Render the README figures to `docs/images/`.
##
## This uses the `plotter` package (github.com/Swarchal/nim-plotter), which is
## *not* a dependency of drc: it is only needed to rebuild these images, and
## `nimble figures` finds it through `--path`. PNG output also needs
## `vl-convert` on PATH.
##
## Each figure proc does three things, in order:
##
## 1. **Data and fit.** These are copied from the matching README example, so
##    the figure shows exactly the fit that the README prints.
## 2. **Dose axis.** The dose range to show, and the unit to show it in.
## 3. **Layers.** Built one per `let` and stacked bottom to top by `layer`.
##
## The helpers turn fits and observations into data frames and chart layers.
## Every data frame has the same three columns: `dose`, `response` and `series`.

import std/[math, os, sequtils, strformat]
import drc
import plotter

# ----------------------------------------------------------------- style ----

const
  width = 520
  height = 320
  pngScale = 2.0 # render at twice the size, for high-DPI screens
  # The first three categorical slots, which stay distinguishable under
  # colour-vision deficiency even when every pair is on screen at once.
  blue = "#2a78d6"
  orange = "#eb6834"
  aqua = "#1baf7a"
  ink = "#52514e" # observations and annotations
  muted = "#9a9893" # reference lines
  dotSize = 40 # area in px², as Vega-Lite measures point size

# ------------------------------------------------------------- dose axis ----

type DoseAxis = object
  ## A log10 dose axis. A fit uses whatever dose unit it was given (e.g. molar);
  ## the axis shows each dose multiplied by `fittedToDisplay` (1e9 shows molar
  ## as nM). `lo` and `hi` are the ends of the axis, in display units.
  title: string
  lo, hi: float
  fittedToDisplay: float

proc doseAxis(title: string, lo, hi: float, fittedToDisplay = 1.0): DoseAxis =
  DoseAxis(title: title, lo: lo, hi: hi, fittedToDisplay: fittedToDisplay)

proc toDisplay(a: DoseAxis, fittedDose: float): float =
  fittedDose * a.fittedToDisplay

proc x(a: DoseAxis): Encoding =
  ## The x encoding: a log scale over exactly [lo, hi], with ticks only at
  ## whole decades. The 1e-9 keeps an end that is already a decade, such as
  ## 1e3, from being rounded past.
  let firstDecade = int(ceil(log10(a.lo) - 1e-9))
  let lastDecade = int(floor(log10(a.hi) + 1e-9))
  var ticks: seq[Value]
  for e in firstDecade .. lastDecade:
    ticks.add v(pow(10.0, float(e)))
  "dose:Q".scale(scaleType = scLog, domain = @[v(a.lo), v(a.hi)], nice = false).axis(
    title = a.title, grid = false, values = ticks, format = "~g"
  )

proc doseGrid(a: DoseAxis, n = 200): seq[float] =
  ## n doses across the whole axis, evenly spaced in log10. They are in
  ## *fitted* units, ready to pass to `predict` or a band.
  let lo = log10(a.lo / a.fittedToDisplay)
  let hi = log10(a.hi / a.fittedToDisplay)
  for i in 0 ..< n:
    result.add pow(10.0, lo + (hi - lo) * float(i) / float(n - 1))

# ----------------------------------------------------------- data frames ----

type Series = object ## A fitted curve and its legend label.
  label: string
  fit: DrcFit

proc series(label: string, fit: DrcFit): Series =
  Series(label: label, fit: fit)

proc curveFrame(a: DoseAxis, curves: openArray[Series]): DataFrame =
  ## Each curve's predicted response across the axis, one row per point, with
  ## the curves stacked one after another.
  var dose, response: seq[float]
  var label: seq[string]
  for curve in curves:
    for d in a.doseGrid:
      dose.add a.toDisplay(d)
      response.add curve.fit.predict(d)
      label.add curve.label
  dataFrame(col("dose", dose), col("response", response), col("series", label))

proc observationFrame(
    a: DoseAxis,
    doses, responses: openArray[float],
    labels: openArray[string] = [],
    keep: openArray[bool] = [],
): DataFrame =
  ## The observations to plot, one row each. `labels` gives each observation's
  ## series (default: none), and `keep` selects which observations to include
  ## (default: all). Zero doses cannot sit on a log axis, so they are always
  ## left out, but the fits still use them.
  var dose, response: seq[float]
  var label: seq[string]
  for i, d in doses:
    let wanted = keep.len == 0 or keep[i]
    if d > 0 and wanted:
      dose.add a.toDisplay(d)
      response.add responses[i]
      label.add (if labels.len > 0: labels[i] else: "")
  dataFrame(col("dose", dose), col("response", response), col("series", label))

# ---------------------------------------------------------------- layers ----

proc seriesColors(labels, colors: openArray[string], legendAt: string): Encoding =
  ## Colour by series, giving each label a fixed colour, so that a series keeps
  ## its colour wherever it appears. `legendAt` is a Vega-Lite legend `orient`.
  "series:N".scale(domain = labels.mapIt(v(it)), range = colors.mapIt(v(it))).legend(
    title = "", orient = legendAt
  )

proc curveLines(
    a: DoseAxis, curves: openArray[Series], yTitle: string, color: Encoding
): Chart =
  ## One line per fitted curve. This layer also sets the y-axis title.
  chart(a.curveFrame(curves)).markLine(strokeWidth = 2).encode(
    x = a.x, y = "response:Q".title(yTitle), color = color
  )

proc dots(a: DoseAxis, observations: DataFrame): Chart =
  ## Observations as solid grey dots.
  chart(observations).markPoint(filled = true, color = ink, size = dotSize).encode(
    x = a.x, y = "response:Q"
  )

proc seriesDots(a: DoseAxis, observations: DataFrame, color: Encoding): Chart =
  ## Observations as solid dots, coloured by series.
  chart(observations).markPoint(filled = true, size = dotSize, opacity = 0.9).encode(
    x = a.x, y = "response:Q", color = color
  )

proc crosses(a: DoseAxis, observations: DataFrame): Chart =
  ## Observations as open grey ×s, for points that a fit has rejected.
  chart(observations)
    .markPoint(
      filled = false,
      color = ink,
      size = 90,
      strokeWidth = 1.5,
      raw = %*{"shape": "cross", "angle": 45},
    )
    .encode(x = a.x, y = "response:Q")

proc finish(c: Chart, title: string): Chart =
  ## The size, title and styling that every figure shares.
  c
    .properties(width = width, height = height, title = title)
    .configureBackground("#ffffff")
    .configureView(noStroke = true)
    .configureAxis(gridColor = "#eeeeee")

# --------------------------------------------------------------- figures ----

proc ryegrassFigure(): Chart =
  ## README, "Fit a curve and print a summary".
  # Ryegrass root length vs ferulic acid concentration (mM), from the R drc package
  let doses = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.94, 0.94, 0.94, 1.88, 1.88, 1.88, 3.75, 3.75, 3.75,
    7.5, 7.5, 7.5, 15.0, 15.0, 15.0, 30.0, 30.0, 30.0,
  ]
  let rootLength = [
    7.58, 8.0, 8.3285714, 7.25, 7.375, 7.9625, 8.3555556, 6.9142857, 7.75, 6.8714286,
    6.45, 5.9222222, 1.925, 2.8857143, 4.2333333, 1.1875, 0.8571429, 1.0571429, 0.6875,
    0.525, 0.825, 0.25, 0.22, 0.44,
  ]
  let fit = fitLL4(doses, rootLength)

  let axis = doseAxis("Ferulic acid (mM)", 0.5, 50.0)

  # The fitted line and both bands come from `confidenceBand` and
  # `predictionBand`, so they share one frame with a column for each edge.
  let grid = axis.doseGrid
  let confidence = fit.confidenceBand(grid)
  let prediction = fit.predictionBand(grid)
  let bands = dataFrame(
    col("dose", grid.mapIt(axis.toDisplay(it))),
    col("fitted", confidence.mapIt(it.fit)),
    col("confidenceLower", confidence.mapIt(it.lower)),
    col("confidenceUpper", confidence.mapIt(it.upper)),
    col("predictionLower", prediction.mapIt(it.lower)),
    col("predictionUpper", prediction.mapIt(it.upper)),
  )
  let ec50 = dataFrame(
    col("dose", @[axis.toDisplay(fit.ec50)]),
    col("label", @[&"EC50 = {fit.ec50:.2f} mM"]),
  )

  let predictionArea = chart(bands).markArea(color = blue, opacity = 0.10).encode(
      x = axis.x,
      y = "predictionLower:Q".title("Root length (cm)"),
      y2 = "predictionUpper:Q",
    )
  let confidenceArea = chart(bands).markArea(color = blue, opacity = 0.25).encode(
      x = axis.x, y = "confidenceLower:Q", y2 = "confidenceUpper:Q"
    )
  let fittedLine = chart(bands).markLine(color = blue, strokeWidth = 2).encode(
      x = axis.x, y = "fitted:Q"
    )
  let ec50Rule =
    chart(ec50).markRule(color = muted, strokeDash = @[4.0, 4.0]).encode(x = axis.x)
  # With no y encoding the label sits at mid-height; dy lifts it near the top.
  let ec50Label = chart(ec50)
    .markText(align = "left", dx = 6, dy = -100, color = ink, fontSize = 12)
    .encode(x = axis.x, text = "label:N")
  let observations = axis.dots(axis.observationFrame(doses, rootLength))

  layer(predictionArea, confidenceArea, fittedLine, ec50Rule, ec50Label, observations)
    .finish("4PL fit with 95% confidence and prediction bands")

proc asymmetricFigure(): Chart =
  ## README, "Asymmetric curves (5PL)".
  var doses = @[0.0, 0.0] # two untreated controls
  for k in 0 .. 12: # 1 nM to 1 mM in half-log steps, duplicates
    let dose = pow(10.0, -9.0 + 0.5 * float(k))
    doses.add [dose, dose]
  let responses = [
    5.0, 5.6, 8.34, 7.1, 10.65, 9.57, 16.18, 18.74, 22.64, 22.38, 37.01, 36.75, 54.87,
    52.8, 75.55, 77.0, 86.17, 87.94, 89.67, 90.89, 90.97, 94.18, 92.39, 95.46, 95.3,
    94.61, 89.96, 93.92,
  ]
  let four = fitLL4(doses, responses)
  let five = fitLL5(doses, responses)

  let axis = doseAxis("Dose (µM)", 7e-4, 1.5e3, fittedToDisplay = 1e6)

  let curves = [
    series(&"4PL (AIC {four.aic:.1f})", four), series(&"5PL (AIC {five.aic:.1f})", five)
  ]
  let color =
    seriesColors(curves.mapIt(it.label), [orange, blue], legendAt = "top-left")

  let lines = axis.curveLines(curves, "Response (%)", color)
  let observations = axis.dots(axis.observationFrame(doses, responses))

  layer(lines, observations).finish("An asymmetric curve: 4PL vs 5PL")

proc robustFigure(): Chart =
  ## README, "Robust fitting (outliers)".
  # % activity of an inhibitor, triplicates, with three wells gone wrong
  let doses = [
    1e-9, 1e-9, 1e-9, 3.16e-9, 3.16e-9, 3.16e-9, 1e-8, 1e-8, 1e-8, 3.16e-8, 3.16e-8,
    3.16e-8, 1e-7, 1e-7, 1e-7, 3.16e-7, 3.16e-7, 3.16e-7, 1e-6, 1e-6, 1e-6, 3.16e-6,
    3.16e-6, 3.16e-6, 1e-5, 1e-5, 1e-5, 3.16e-5, 3.16e-5, 3.16e-5, 1e-4, 1e-4, 1e-4,
  ]
  let activity = [
    104.03, 90.24, 99.16, 95.92, 136.26, 96.97, 90.44, 95.81, 93.91, 102.27, 92.98,
    91.24, 77.88, 76.72, 75.56, 48.83, 16.45, 49.28, 24.15, 20.67, 21.35, 12.34, 9.33,
    6.18, 2.95, 35.12, 9.3, 1.57, 1.65, 5.39, -0.56, 1.22, 4.74,
  ]
  var curves: seq[Series]
  for loss in [Squared, Huber, Tukey]:
    let fit = fitLL4(doses, activity, opts = FitOptions(loss: loss))
    curves.add series(&"{loss} (IC50 {round(fit.ec50 * 1e9).int} nM)", fit)

  # Tukey's bisquare gives weight exactly 0 to a residual beyond its cut-off,
  # so these observations play no part in the Tukey fit.
  let tukey = curves[^1].fit
  let rejected = tukey.robustWeights.mapIt(it == 0)
  let used = rejected.mapIt(not it)

  let axis = doseAxis("Inhibitor (nM)", 0.7, 1.5e5, fittedToDisplay = 1e9)
  let color =
    seriesColors(curves.mapIt(it.label), [orange, aqua, blue], legendAt = "top-right")

  let lines = axis.curveLines(curves, "Activity (% of control)", color)
  let usedObservations = axis.dots(axis.observationFrame(doses, activity, keep = used))
  let rejectedObservations =
    axis.crosses(axis.observationFrame(doses, activity, keep = rejected))

  layer(lines, usedObservations, rejectedObservations).finish(
    "Least squares vs robust fits (× = given weight 0 by Tukey)"
  )

proc multiFigure(): Chart =
  ## README, "Several curves: shared parameters and relative potency".
  # R's DNase assay, runs 1, 4 and 8: duplicate wells at eight concentrations
  let conc = [0.04882812, 0.1953125, 0.390625, 0.78125, 1.5625, 3.125, 6.25, 12.5]
  let density = [
    0.017, 0.018, 0.121, 0.124, 0.206, 0.215, 0.377, 0.374, 0.614, 0.609, 1.019, 1.001,
    1.334, 1.364, 1.73, 1.71, 0.011, 0.016, 0.118, 0.108, 0.2, 0.206, 0.364, 0.36, 0.62,
    0.64, 0.979, 0.973, 1.424, 1.399, 1.74, 1.732, 0.054, 0.054, 0.152, 0.148, 0.226,
    0.222, 0.392, 0.383, 0.658, 0.644, 1.043, 1.002, 1.466, 1.381, 1.743, 1.724,
  ]
  var doses: seq[float]
  var curve: seq[int] # which curve each observation belongs to
  for c in 0 .. 2:
    for d in conc:
      doses.add [d, d]
      curve.add [c, c]
  let parallel = fitLL4Multi(doses, density, curve, shared = {Bottom, Top, Hill})

  let axis = doseAxis("DNase concentration (ng/ml)", 0.03, 20.0)

  # `parallel.curves[c]` is curve c's own fit, so it plots like any other.
  let runs = ["run 1", "run 4", "run 8"]
  var curves: seq[Series]
  for c, run in runs:
    curves.add series(run, parallel.curves[c])
  let runOfObservation = curve.mapIt(runs[it])
  let color = seriesColors(runs, [blue, orange, aqua], legendAt = "top-left")

  let lines = axis.curveLines(curves, "Optical density", color)
  let observations = axis.seriesDots(
    axis.observationFrame(doses, density, labels = runOfObservation), color
  )

  layer(lines, observations).finish("Parallel 4PL curves (shared Bottom, Top and Hill)")

# ------------------------------------------------------------------ main ----

when isMainModule:
  if findExe(vlConvertExe).len == 0:
    quit(&"`{vlConvertExe}` is not on PATH; it is needed to render PNGs.", 1)
  let outDir = currentSourcePath().parentDir.parentDir / "docs" / "images"
  createDir(outDir)
  # In README order. Each name is the PNG that the README links to.
  for (name, figure) in [
    ("ryegrass", ryegrassFigure()),
    ("asymmetric", asymmetricFigure()),
    ("robust", robustFigure()),
    ("multi", multiFigure()),
  ]:
    figure.savePng(outDir / &"{name}.png", scale = pngScale)
    echo &"  {name}.png"
