## Dose-response models, written in log10-dose space (x = log10(dose)).
## A zero dose maps to x = -Inf and evaluates to the low-dose plateau, so
## untreated controls can be fitted directly.

import std/[algorithm, math, stats]
import linalg

type DrcParam* = enum
  Bottom
  Top
  LogEC50
  Hill
  Asymmetry

const drcParamCount* = ord(high(DrcParam)) + 1

const ln10 = ln(10.0)

func log1p(x: float): float {.importc, header: "<math.h>".}
func expm1(x: float): float {.importc, header: "<math.h>".}

func logistic(u: float): tuple[s, sc, lnS: float] =
  ## s = 1/(1+eᵘ), its complement sc = 1 − s, and ln s, computed without
  ## overflow or cancellation.
  if u > 0.0:
    let e = exp(-u)
    (e / (1.0 + e), 1.0 / (1.0 + e), -u - log1p(e))
  else:
    let e = exp(u)
    (1.0 / (1.0 + e), e / (1.0 + e), -log1p(e))

proc ll5*(
    x: float, p: array[drcParamCount, float], grad: var array[drcParamCount, float]
): float {.nimcall.} =
  ## Five-parameter logistic:
  ##   f(x) = Bottom + (Top − Bottom) / (1 + 10^((LogEC50 − x) · Hill))^Asymmetry
  ## With Asymmetry = 1 this is the four-parameter logistic (Hill equation),
  ## and LogEC50 is the dose of the half-way response; otherwise LogEC50 is a
  ## location parameter and the half-way dose differs from it.
  ##
  ## Hill > 0 gives a rising curve, Hill < 0 a falling one. For the 4PL,
  ## swapping Bottom and Top and negating Hill gives the same curve;
  ## `drcGuess` starts fits from the Bottom < Top form, so fitted parameters
  ## normally come out that way. Asymmetry must be positive (NaN otherwise,
  ## which optimisers treat as a failed step).
  let
    b = p[ord(Bottom)]
    t = p[ord(Top)]
    l = p[ord(LogEC50)]
    h = p[ord(Hill)]
    a = p[ord(Asymmetry)]
  if not (a > 0.0):
    grad = default(array[drcParamCount, float])
    return NaN
  if x == NegInf and h != 0.0:
    # Zero dose: the limit is Bottom for a rising curve and Top for a falling one.
    grad = default(array[drcParamCount, float])
    if h > 0.0:
      grad[ord(Bottom)] = 1.0
      return b
    grad[ord(Top)] = 1.0
    return t

  # With Hill = 0 the curve is flat at z = 0, zero dose included (where
  # (LogEC50 − x) · Hill would be ∞ · 0).
  let
    zeroDose = x == NegInf
    (_, sc, lnS) = logistic(
      if zeroDose:
        0.0
      else:
        ln10 * (l - x) * h
    )
    sa = exp(a * lnS) # sᴬ, the fraction of the way to Top
    dfdz = -(t - b) * a * sa * sc * ln10 # ∂f/∂z with z = (LogEC50 − x) · Hill
  grad[ord(Bottom)] = -expm1(a * lnS) # 1 − sᴬ
  grad[ord(Top)] = sa
  grad[ord(LogEC50)] = dfdz * h
  grad[ord(Hill)] =
    if zeroDose:
      0.0
    else:
      dfdz * (l - x)
  grad[ord(Asymmetry)] = (t - b) * sa * lnS
  b + (t - b) * sa

proc drcGuess*(
    x, y: openArray[float],
    fixedValues: array[drcParamCount, float] = default(array[drcParamCount, float]),
    fixed: set[DrcParam] = {},
): array[drcParamCount, float] =
  ## Starting values for a fit, using parameters in `fixed` as given. The
  ## curve is assumed symmetric (Asymmetry = 1) unless Asymmetry is fixed.
  ##
  ## * Bottom/Top: the smallest and largest mean response over replicate groups.
  ## * Hill direction: sign of the linear regression slope of y on x.
  ## * LogEC50: where the group means cross (Bottom + Top) / 2, by linear interpolation.
  ## * Hill magnitude: slope of the linearised curve log10((Top − y)/(y − Bottom)) vs x.
  assert x.len == y.len and x.len > 0

  # Mean response for each distinct dose, ordered by dose.
  var points = newSeq[(float, float)](x.len)
  for i in 0 ..< x.len:
    points[i] = (x[i], y[i])
  points.sort()
  var gx, gy: seq[float]
  var i = 0
  while i < points.len:
    var j = i
    var total = 0.0
    while j < points.len and points[j][0] == points[i][0]:
      total += points[j][1]
      inc j
    gx.add points[i][0]
    gy.add total / float(j - i)
    i = j

  var b = min(gy)
  var t = max(gy)
  if Bottom in fixed:
    b = fixedValues[ord(Bottom)]
  if Top in fixed:
    t = fixedValues[ord(Top)]

  # Only finite doses carry information about position and slope.
  var fx, fy: seq[float]
  for k in 0 ..< gx.len:
    if gx[k].isFinite:
      fx.add gx[k]
      fy.add gy[k]

  var trend: RunningRegress
  trend.push(fx, fy)
  let direction = if trend.slope < 0.0: -1.0 else: 1.0

  var l = NaN
  let mid = 0.5 * (b + t)
  for k in 0 ..< fx.len - 1:
    let (y0, y1) = (fy[k], fy[k + 1])
    if y0 != y1 and (y0 - mid) * (y1 - mid) <= 0.0:
      l = fx[k] + (mid - y0) * (fx[k + 1] - fx[k]) / (y1 - y0)
      break
  if l.isNaN:
    l =
      if fx.len > 0:
        0.5 * (fx[0] + fx[^1])
      else:
        0.0

  var h = direction
  if t != b:
    # (Top − y)/(y − Bottom) = 10^((LogEC50 − x)·Hill), so z is linear in x with slope −Hill.
    let margin = 0.05 * abs(t - b)
    var linearised: RunningRegress
    for k in 0 ..< fx.len:
      if fy[k] > min(b, t) + margin and fy[k] < max(b, t) - margin:
        linearised.push(fx[k], log10((t - fy[k]) / (fy[k] - b)))
    # NaN or infinite when fewer than two distinct doses were usable.
    let est = -linearised.slope
    if est.isFinite and est * direction > 0.0:
      h = direction * clamp(abs(est), 0.1, 10.0)

  result[ord(Bottom)] = b
  result[ord(Top)] = t
  result[ord(LogEC50)] = l
  result[ord(Hill)] = h
  result[ord(Asymmetry)] = 1.0
  for param in fixed:
    result[ord(param)] = fixedValues[ord(param)]
