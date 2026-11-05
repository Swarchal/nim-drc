## Plotting helpers for drc fits, shared by the vignettes and the README
## figures (`examples/figures.nim`).
##
## This uses the `plotter` package (github.com/Swarchal/nim-plotter), which is
## *not* a dependency of drc: it is only needed to draw figures, and the
## `nimble vignettes`/`nimble figures` tasks find it through `--path`. PNG
## output also needs `vl-convert` on PATH.
##
## Every data frame built here has the same three columns: `dose`, `response`
## and `series`. A figure is a `layer` of charts, finished with `finish`.

import std/[math, os, sequtils, strformat]
import drc
import plotter

export plotter

# ----------------------------------------------------------------- style ----

const
  width* = 520
  height* = 320
  pngScale* = 2.0 # render at twice the size, for high-DPI screens
  # The first three categorical slots, which stay distinguishable under
  # colour-vision deficiency even when every pair is on screen at once.
  blue* = "#2a78d6"
  orange* = "#eb6834"
  aqua* = "#1baf7a"
  ink* = "#52514e" # observations and annotations
  muted* = "#9a9893" # reference lines
  dotSize* = 40 # area in px², as Vega-Lite measures point size

# ------------------------------------------------------------- dose axis ----

type DoseAxis* = object
  ## A log10 dose axis. A fit uses whatever dose unit it was given (e.g. molar);
  ## the axis shows each dose multiplied by `fittedToDisplay` (1e9 shows molar
  ## as nM). `lo` and `hi` are the ends of the axis, in display units.
  title*: string
  lo*, hi*: float
  fittedToDisplay*: float

proc doseAxis*(title: string, lo, hi: float, fittedToDisplay = 1.0): DoseAxis =
  DoseAxis(title: title, lo: lo, hi: hi, fittedToDisplay: fittedToDisplay)

proc toDisplay*(a: DoseAxis, fittedDose: float): float =
  fittedDose * a.fittedToDisplay

proc x*(a: DoseAxis): Encoding =
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

proc doseGrid*(a: DoseAxis, n = 200): seq[float] =
  ## n doses across the whole axis, evenly spaced in log10. They are in
  ## *fitted* units, ready to pass to `predict` or a band.
  let lo = log10(a.lo / a.fittedToDisplay)
  let hi = log10(a.hi / a.fittedToDisplay)
  for i in 0 ..< n:
    result.add pow(10.0, lo + (hi - lo) * float(i) / float(n - 1))

# ----------------------------------------------------------- data frames ----

type LabelledFit* = object ## A fitted curve and its legend label.
  label*: string
  fit*: DrcFit

proc labelled*(label: string, fit: DrcFit): LabelledFit =
  LabelledFit(label: label, fit: fit)

proc curveFrame*(a: DoseAxis, curves: openArray[LabelledFit]): DataFrame =
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

proc observationFrame*(
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

proc seriesColors*(labels, colors: openArray[string], legendAt: string): Encoding =
  ## Colour by series, giving each label a fixed colour, so that a series keeps
  ## its colour wherever it appears. `legendAt` is a Vega-Lite legend `orient`.
  "series:N".scale(domain = labels.mapIt(v(it)), range = colors.mapIt(v(it))).legend(
    title = "", orient = legendAt
  )

proc curveLines*(
    a: DoseAxis, curves: openArray[LabelledFit], yTitle: string, color: Encoding
): Chart =
  ## One line per fitted curve. This layer also sets the y-axis title.
  chart(a.curveFrame(curves)).markLine(strokeWidth = 2).encode(
    x = a.x, y = "response:Q".title(yTitle), color = color
  )

proc curveLine*(a: DoseAxis, fit: DrcFit, yTitle: string, color = blue): Chart =
  ## A single fitted curve in one colour. This layer also sets the y-axis title.
  chart(a.curveFrame([labelled("", fit)]))
    .markLine(color = color, strokeWidth = 2)
    .encode(x = a.x, y = "response:Q".title(yTitle))

proc bands*(a: DoseAxis, fit: DrcFit, color = blue, prediction = true): Chart =
  ## The fit's 95% confidence band, and (with `prediction`) its wider 95%
  ## prediction band behind it, as shaded areas.
  let grid = a.doseGrid
  let dose = grid.mapIt(a.toDisplay(it))
  let confidence = fit.confidenceBand(grid)
  let confidenceArea = chart(
      dataFrame(
        col("dose", dose),
        col("lower", confidence.mapIt(it.lower)),
        col("upper", confidence.mapIt(it.upper)),
      )
    )
    .markArea(color = color, opacity = 0.25)
    .encode(x = a.x, y = "lower:Q", y2 = "upper:Q")
  if not prediction:
    return confidenceArea
  let pred = fit.predictionBand(grid)
  let predictionArea = chart(
      dataFrame(
        col("dose", dose),
        col("predictionLower", pred.mapIt(it.lower)),
        col("predictionUpper", pred.mapIt(it.upper)),
      )
    )
    .markArea(color = color, opacity = 0.10)
    .encode(x = a.x, y = "predictionLower:Q", y2 = "predictionUpper:Q")
  layer(predictionArea, confidenceArea)

proc dots*(a: DoseAxis, observations: DataFrame): Chart =
  ## Observations as solid grey dots.
  chart(observations).markPoint(filled = true, color = ink, size = dotSize).encode(
    x = a.x, y = "response:Q"
  )

proc seriesDots*(a: DoseAxis, observations: DataFrame, color: Encoding): Chart =
  ## Observations as solid dots, coloured by series.
  chart(observations).markPoint(filled = true, size = dotSize, opacity = 0.9).encode(
    x = a.x, y = "response:Q", color = color
  )

proc crosses*(a: DoseAxis, observations: DataFrame): Chart =
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

proc doseRule*(a: DoseAxis, dose: float, label = "", dy = -100.0): Chart =
  ## A dashed vertical line at `dose` (in fitted units), labelled at its right.
  ## With no y encoding the label sits at mid-height; `dy` moves it up or down.
  let at = @[a.toDisplay(dose)]
  let rule = chart(dataFrame(col("dose", at)))
    .markRule(color = muted, strokeDash = @[4.0, 4.0])
    .encode(x = a.x)
  if label.len == 0:
    return rule
  let text = chart(dataFrame(col("dose", at), col("label", @[label])))
    .markText(align = "left", dx = 6, dy = dy, color = ink, fontSize = 12)
    .encode(x = a.x, text = "label:N")
  layer(rule, text)

proc responseRule*(response: float): Chart =
  ## A dashed horizontal line at `response`, across the whole plot.
  chart(dataFrame(col("response", @[response])))
    .markRule(color = muted, strokeDash = @[4.0, 4.0])
    .encode(y = "response:Q")

proc finish*(c: Chart, title: string, w = width, h = height): Chart =
  ## The size, title and styling that every figure shares.
  c
    .properties(width = w, height = h, title = title)
    .configureBackground("#ffffff")
    .configureView(noStroke = true)
    .configureAxis(gridColor = "#eeeeee")

# --------------------------------------------------------------- figures ----

proc figure*(c: Chart, name: string, alt = "") =
  ## Render `c` to `<name>.png`, and print a marker line that the vignette
  ## renderer (`vignettes/knit.nim`) turns into an image link. The directory
  ## is `$DRC_FIGURE_DIR`, or the current directory if that is not set.
  let dir = getEnv("DRC_FIGURE_DIR", ".")
  createDir(dir)
  c.savePng(dir / name & ".png", scale = pngScale)
  echo &"\x1Dfigure {name}.png {alt}"
