## Render the README figures to `docs/images/`.
##
## This uses the `plotter` package (github.com/Swarchal/nim-plotter), which is
## *not* a dependency of drc: it is only needed to rebuild these images, and
## `nimble figures` finds it through `--path`. PNG output also needs
## `vl-convert` on PATH.
##
## The plotting helpers (`doseAxis`, `curveLines`, `dots`, `finish`, …) are in
## `vignettes/drcplot.nim`, which the vignettes share.
##
## Each figure proc does three things, in order:
##
## 1. **Data and fit.** These are copied from the matching README example, so
##    the figure shows exactly the fit that the README prints.
## 2. **Dose axis.** The dose range to show, and the unit to show it in.
## 3. **Layers.** Built one per `let` and stacked bottom to top by `layer`.

import std/[math, os, sequtils, strformat]
import drc
import drcplot

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
    labelled(&"4PL (AIC {four.aic:.1f})", four),
    labelled(&"5PL (AIC {five.aic:.1f})", five),
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
  var curves: seq[LabelledFit]
  for loss in [Squared, Huber, Tukey]:
    let fit = fitLL4(doses, activity, opts = FitOptions(loss: loss))
    curves.add labelled(&"{loss} (IC50 {round(fit.ec50 * 1e9).int} nM)", fit)

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
  # Simulated potency assay: a reference standard and two test samples, each at
  # nine half-log concentrations (0.01 to 100 ng/ml) in duplicate.
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
  let parallel = fitLL4Multi(doses, response, curve, shared = {Bottom, Top, Hill})

  let axis = doseAxis("Concentration (ng/ml)", 0.006, 160.0)

  # `parallel.curves[c]` is curve c's own fit, so it plots like any other.
  let samples = ["reference", "sample A", "sample B"]
  var curves: seq[LabelledFit]
  for c, sample in samples:
    curves.add labelled(sample, parallel.curves[c])
  let sampleOfObservation = curve.mapIt(samples[it])
  let color = seriesColors(samples, [blue, orange, aqua], legendAt = "top-left")

  let lines = axis.curveLines(curves, "Response", color)
  let observations = axis.seriesDots(
    axis.observationFrame(doses, response, labels = sampleOfObservation), color
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
