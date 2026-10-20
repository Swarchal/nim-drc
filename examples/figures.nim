## Render the README figures to `docs/images/`.
##
## This uses the `plotter` package (github.com/Swarchal/nim-plotter), which is
## *not* a dependency of drc: it is only needed to rebuild these images, and
## `nimble figures` finds it through `--path`. PNG output also needs
## `vl-convert` on PATH. The data are the README's examples, so the figures
## show the same fits that the README prints.

import std/[math, os, strformat]
import drc
import plotter

const
  width = 520
  height = 320
  pngScale = 2.0
  # The first three categorical slots, which stay distinguishable under
  # colour-vision deficiency even when every pair is on screen at once.
  blue = "#2a78d6"
  orange = "#eb6834"
  aqua = "#1baf7a"
  ink = "#52514e"
  muted = "#9a9893"

proc logGrid(lo, hi: float; n = 200): seq[float] =
  ## n doses evenly spaced in log10 between lo and hi.
  let (a, b) = (log10(lo), log10(hi))
  for i in 0 ..< n:
    result.add pow(10.0, a + (b - a) * float(i) / float(n - 1))

proc nonZero(doses, resp: openArray[float]): (seq[float], seq[float]) =
  ## Zero doses cannot sit on a log axis, so the plots leave them out
  ## (the fits still use them).
  for i, d in doses:
    if d > 0:
      result[0].add d
      result[1].add resp[i]

proc logX(title: string; lo, hi: float): Encoding =
  ## A log10 dose axis over [lo, hi], labelled at whole decades only.
  var decades: seq[Value]
  for e in int(ceil(log10(lo) - 1e-9)) .. int(floor(log10(hi) + 1e-9)):
    decades.add v(pow(10.0, float(e)))
  "dose:Q".scale(scaleType = scLog, domain = @[v(lo), v(hi)], nice = false)
          .axis(title = title, grid = false, values = decades, format = "~g")

proc finish(c: Chart): Chart =
  c.properties(width = width, height = height)
   .configureBackground("#ffffff")
   .configureView(noStroke = true)
   .configureAxis(gridColor = "#eeeeee")

proc ryegrassFigure(): Chart =
  let doses = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.94, 0.94, 0.94, 1.88, 1.88, 1.88,
               3.75, 3.75, 3.75, 7.5, 7.5, 7.5, 15.0, 15.0, 15.0, 30.0, 30.0, 30.0]
  let rootLength = [7.58, 8.0, 8.3285714, 7.25, 7.375, 7.9625,
                    8.3555556, 6.9142857, 7.75, 6.8714286, 6.45, 5.9222222,
                    1.925, 2.8857143, 4.2333333, 1.1875, 0.8571429, 1.0571429,
                    0.6875, 0.525, 0.825, 0.25, 0.22, 0.44]
  let fit = fitLL4(doses, rootLength)
  let grid = logGrid(0.5, 50.0)
  let conf = fit.confidenceBand(grid)
  let pred = fit.predictionBand(grid)
  var fitted, confLo, confHi, predLo, predHi: seq[float]
  for i in 0 ..< grid.len:
    fitted.add conf[i].fit
    confLo.add conf[i].lower
    confHi.add conf[i].upper
    predLo.add pred[i].lower
    predHi.add pred[i].upper
  let curve = dataFrame(col("dose", grid), col("fit", fitted),
                        col("confLo", confLo), col("confHi", confHi),
                        col("predLo", predLo), col("predHi", predHi))
  let (pd, pr) = nonZero(doses, rootLength)
  let points = dataFrame(col("dose", pd), col("response", pr))
  let x = logX("Ferulic acid (mM)", 0.5, 50.0)

  let prediction = chart(curve).markArea(color = blue, opacity = 0.10).encode(x = x, y = "predLo:Q".title("Root length (cm)"), y2 = "predHi:Q")
  let confidence = chart(curve).markArea(color = blue, opacity = 0.25).encode(x = x, y = "confLo:Q", y2 = "confHi:Q")
  let line = chart(curve).markLine(color = blue, strokeWidth = 2).encode(x = x, y = "fit:Q")
  let obs = chart(points).markPoint(filled = true, color = ink, size = 40).encode(x = x, y = "response:Q")
  let ec50 = dataFrame(col("dose", @[fit.ec50]), col("label", @[&"EC50 = {fit.ec50:.2f} mM"]))
  let rule = chart(ec50).markRule(color = muted, strokeDash = @[4.0, 4.0]).encode(x = x)
  let label = chart(ec50).markText(align = "left", dx = 6, dy = -100, color = ink, fontSize = 12).encode(x = x, text = "label:N")
  layer(prediction, confidence, line, rule, label, obs).finish().properties(title = "4PL fit with 95% confidence and prediction bands")

proc robustFigure(): Chart =
  let doses = [1e-9, 1e-9, 1e-9, 3.16e-9, 3.16e-9, 3.16e-9, 1e-8, 1e-8, 1e-8,
               3.16e-8, 3.16e-8, 3.16e-8, 1e-7, 1e-7, 1e-7, 3.16e-7, 3.16e-7, 3.16e-7,
               1e-6, 1e-6, 1e-6, 3.16e-6, 3.16e-6, 3.16e-6, 1e-5, 1e-5, 1e-5,
               3.16e-5, 3.16e-5, 3.16e-5, 1e-4, 1e-4, 1e-4]
  let activity = [104.03, 90.24, 99.16, 95.92, 136.26, 96.97, 90.44, 95.81, 93.91,
                  102.27, 92.98, 91.24, 77.88, 76.72, 75.56, 48.83, 16.45, 49.28,
                  24.15, 20.67, 21.35, 12.34, 9.33, 6.18, 2.95, 35.12, 9.3,
                  1.57, 1.65, 5.39, -0.56, 1.22, 4.74]
  let grid = logGrid(0.7e-9, 1.5e-4)
  var gDose, gResp: seq[float]
  var gLoss: seq[string]
  var tukey: DrcFit
  for loss in [Squared, Huber, Tukey]:
    let fit = fitLL4(doses, activity, opts = FitOptions(loss: loss))
    if loss == Tukey: tukey = fit
    let name = &"{loss} (IC50 {round(fit.ec50 * 1e9).int} nM)"
    for d in grid:
      gDose.add d * 1e9
      gResp.add fit.predict(d)
      gLoss.add name
  let curves = dataFrame(col("dose", gDose), col("activity", gResp), col("fit", gLoss))
  var kept, rejected: seq[(float, float)]
  for i, w in tukey.robustWeights:
    if w == 0: rejected.add (doses[i] * 1e9, activity[i])
    else: kept.add (doses[i] * 1e9, activity[i])
  proc frame(pts: seq[(float, float)]): DataFrame =
    var d, a: seq[float]
    for (x, y) in pts:
      d.add x
      a.add y
    dataFrame(col("dose", d), col("activity", a))
  let x = logX("Inhibitor (nM)", 0.7, 1.5e5)
  let names = @[v(gLoss[0]), v(gLoss[grid.len]), v(gLoss[2 * grid.len])]

  let lines = chart(curves).markLine(strokeWidth = 2)
    .encode(x = x, y = "activity:Q".title("Activity (% of control)"),
            color = "fit:N".scale(domain = names, range = @[v(orange), v(aqua), v(blue)])
                           .legend(title = "", orient = "top-right"))
  let ok = chart(frame(kept)).markPoint(filled = true, color = ink, size = 40)
    .encode(x = x, y = "activity:Q")
  let bad = chart(frame(rejected)).markPoint(filled = false, color = ink, size = 90,
                                             strokeWidth = 1.5, raw = %*{"shape": "cross", "angle": 45})
    .encode(x = x, y = "activity:Q")
  layer(lines, ok, bad).finish()
    .properties(title = "Least squares vs robust fits (× = given weight 0 by Tukey)")

proc multiFigure(): Chart =
  let conc = [0.04882812, 0.1953125, 0.390625, 0.78125, 1.5625, 3.125, 6.25, 12.5]
  let density = [0.017, 0.018, 0.121, 0.124, 0.206, 0.215, 0.377, 0.374,
                 0.614, 0.609, 1.019, 1.001, 1.334, 1.364, 1.73, 1.71,
                 0.011, 0.016, 0.118, 0.108, 0.2, 0.206, 0.364, 0.36,
                 0.62, 0.64, 0.979, 0.973, 1.424, 1.399, 1.74, 1.732,
                 0.054, 0.054, 0.152, 0.148, 0.226, 0.222, 0.392, 0.383,
                 0.658, 0.644, 1.043, 1.002, 1.466, 1.381, 1.743, 1.724]
  var doses: seq[float]
  var curve: seq[int]
  for c in 0 .. 2:
    for d in conc:
      doses.add [d, d]
      curve.add [c, c]
  let parallel = fitLL4Multi(doses, density, curve, shared = {Bottom, Top, Hill})
  let runs = ["run 1", "run 4", "run 8"]
  let grid = logGrid(0.03, 20.0)
  var gDose, gResp: seq[float]
  var gRun: seq[string]
  for c in 0 .. 2:
    for d in grid:
      gDose.add d
      gResp.add parallel.curves[c].predict(d)
      gRun.add runs[c]
  var pRun: seq[string]
  for c in curve: pRun.add runs[c]
  let x = logX("DNase concentration (ng/ml)", 0.03, 20.0)
  let color = "run:N".scale(domain = @[v(runs[0]), v(runs[1]), v(runs[2])],
                            range = @[v(blue), v(orange), v(aqua)])
                     .legend(title = "", orient = "top-left")
  let lines = chart(dataFrame(col("dose", gDose), col("density", gResp), col("run", gRun)))
    .markLine(strokeWidth = 2)
    .encode(x = x, y = "density:Q".title("Optical density"), color = color)
  let pts = chart(dataFrame(col("dose", doses), col("density", @density), col("run", pRun)))
    .markPoint(filled = true, size = 40, opacity = 0.9)
    .encode(x = x, y = "density:Q", color = color)
  layer(lines, pts).finish()
    .properties(title = "Parallel 4PL curves (shared Bottom, Top and Hill)")

proc asymmetricFigure(): Chart =
  var doses = @[0.0, 0.0]
  for k in 0 .. 12:
    let dose = pow(10.0, -9.0 + 0.5 * float(k))
    doses.add [dose, dose]
  let responses = [5.0, 5.6, 8.34, 7.1, 10.65, 9.57, 16.18, 18.74, 22.64, 22.38, 37.01, 36.75,
                   54.87, 52.8, 75.55, 77.0, 86.17, 87.94, 89.67, 90.89, 90.97, 94.18, 92.39,
                   95.46, 95.3, 94.61, 89.96, 93.92]
  let four = fitLL4(doses, responses)
  let five = fitLL5(doses, responses)
  let grid = logGrid(7e-10, 1.5e-3)
  var gDose, gResp: seq[float]
  var gModel: seq[string]
  let names = [&"4PL (AIC {four.aic:.1f})", &"5PL (AIC {five.aic:.1f})"]
  for (name, fit) in [(names[0], four), (names[1], five)]:
    for d in grid:
      gDose.add d * 1e6
      gResp.add fit.predict(d)
      gModel.add name
  let (pd, pr) = nonZero(doses, responses)
  var pdUm: seq[float]
  for d in pd: pdUm.add d * 1e6
  let x = logX("Dose (µM)", 7e-4, 1.5e3)
  let lines = chart(dataFrame(col("dose", gDose), col("response", gResp), col("model", gModel)))
    .markLine(strokeWidth = 2)
    .encode(x = x, y = "response:Q".title("Response (%)"),
            color = "model:N".scale(domain = @[v(names[0]), v(names[1])],
                                    range = @[v(orange), v(blue)])
                             .legend(title = "", orient = "top-left"))
  let pts = chart(dataFrame(col("dose", pdUm), col("response", pr)))
    .markPoint(filled = true, color = ink, size = 40)
    .encode(x = x, y = "response:Q")
  layer(lines, pts).finish()
    .properties(title = "An asymmetric curve: 4PL vs 5PL")

when isMainModule:
  if findExe(vlConvertExe).len == 0:
    quit(&"`{vlConvertExe}` is not on PATH; it is needed to render PNGs.", 1)
  let outDir = currentSourcePath().parentDir.parentDir / "docs" / "images"
  createDir(outDir)
  for (name, c) in [("ryegrass", ryegrassFigure()), ("robust", robustFigure()),
                    ("multi", multiFigure()), ("asymmetric", asymmetricFigure())]:
    c.savePng(outDir / &"{name}.png", scale = pngScale)
    echo &"  {name}.png"
