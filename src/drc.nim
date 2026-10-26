## Dose-response curve fitting.
##
## .. code-block:: nim
##   let fit = fitLL4(doses, responses, fixed = [(Bottom, 0.0)])
##   # or, with the model chosen at run time: fitCurve(LL4, doses, responses, ...)
##   echo fit.ec50, " ", fit.ec50CI
##   echo fit.summary

import std/[algorithm, math, strformat]
from std/stats import
  RunningStat, push, mean, varianceS # drc/stats is a different module
import drc/[linalg, lm, models, neldermead, robust, stats]

export lm, models, neldermead, robust, stats

type
  Model* = enum
    LL3 = "3PL" ## 4PL with Bottom fixed (at 0, unless `fixed` gives it)
    LL4 = "4PL"
    LL5 = "5PL"

  FitWarning* = enum
    fwNotConverged = "the optimiser did not converge"
    fwSingular =
      "parameters are not identifiable (JᵀWJ is singular); standard errors unavailable"
    fwNoResidualDf = "no residual degrees of freedom; standard errors unavailable"
    fwEc50OutsideRange = "EC50 lies outside the tested dose range"
    fwRobustNotConverged = "the robust reweighting did not converge"
    fwAtBound =
      "a parameter is at a bound; its standard error and interval ignore the bound"

  FitOptions* = object
    lm*: LMOptions = LMOptions()
    nelderMead*: NMOptions = NMOptions()
    multiStart*: bool = true
      ## also start from shifted LogEC50 and scaled Hill, keeping the best fit
    fallback*: bool = true
      ## run Nelder–Mead when Levenberg–Marquardt does not converge
    loss*: Loss = Squared ## Huber or Tukey for a robust fit that downweights outliers
    tuning*: float = NaN
      ## the loss's tuning constant in units of the robust scale;
      ## NaN uses the 95%-efficient default (1.345 Huber, 4.685 Tukey)
    irlsMaxIter*: int = 100 ## reweighting rounds for a robust fit
    irlsTol*: float = 1e-10
      ## stop reweighting when the residuals change by less than this (relative)

  ConfInt* = tuple[lower, upper: float]

  Band* = tuple[fit, lower, upper: float] ## a fitted response and its interval

  FTest* = tuple[f: float, df1, df2: int, p: float]

  DrcFit* = object
    model*: Model
    params*: array[DrcParam, float]
    se*: array[DrcParam, float]
      ## standard errors: 0 for fixed parameters, NaN if unavailable
    cov*: array[DrcParam, array[DrcParam, float]]
      ## parameter covariance, σ²(JᵀWJ)⁻¹
    fixed*: set[DrcParam] ## includes Asymmetry (= 1) for 3PL and 4PL fits
    lower*, upper*: array[DrcParam, float]
      ## the bounds on each parameter: −Inf and Inf when unbounded
    atBound*: set[DrcParam] ## free parameters that ended on a bound
    n*: int ## observations with non-zero weight
    df*: int ## residual degrees of freedom: n − number of free parameters
    ssr*: float ## weighted residual sum of squares
    sigma*: float
      ## residual standard error, sqrt(ssr / df); for a robust fit, the
      ## robust (MAD) scale of the residuals
    rSquared*: float ## 1 − SSR / (weighted) total sum of squares
    logLik*: float ## Gaussian log-likelihood, as R's logLik.nls (NaN for robust fits)
    aic*: float ## AIC, counting σ as a parameter, as R and drc (NaN for robust fits)
    bic*: float ## (NaN for robust fits)
    loss*: Loss
    tuning*: float ## the loss's tuning constant (robust fits)
    robustWeights*: seq[float]
      ## final ψ(u)/u per observation (robust fits; 0 where the
      ## user's weight is 0); 1 means not downweighted
    logDoseRange*: tuple[lower, upper: float] ## smallest and largest non-zero log10 dose
    converged*: bool
      ## false if the optimiser or, for a robust fit, the reweighting
      ## did not converge (see `warnings`)
    usedFallback*: bool ## Levenberg–Marquardt failed and Nelder–Mead was run
    fallbackConverged*: bool
      ## with `usedFallback`: Nelder–Mead or its polishing run converged
    iterations*: int ## iterations of the final Levenberg–Marquardt run
    status*: LMStatus
      ## how the final Levenberg–Marquardt run stopped; with the
      ## fallback, `converged` may be true when this is not
    warnings*: set[FitWarning]
    doses*, responses*: seq[float] ## the data, as given
    weights*: seq[float] ## the weights, as given (all 1 if none were)
    jointView: bool ## a curve of a joint fit (`MultiFit.curves`)

  ParamVec = array[drcParamCount, float]

  Box = tuple[lower, upper: ParamVec] ## bounds on each parameter, ±Inf if none

func toVec(p: array[DrcParam, float]): ParamVec =
  for param in DrcParam:
    result[ord(param)] = p[param]

func toParams(v: ParamVec): array[DrcParam, float] =
  for param in DrcParam:
    result[param] = v[ord(param)]

proc checkDose(dose: float) =
  if not (dose >= 0.0 and dose < Inf):
    raise newException(ValueError, "dose must be finite and non-negative, got " & $dose)

proc predict*(fit: DrcFit, dose: float): float =
  ## The fitted mean response at `dose` (0 gives the low-dose plateau).
  checkDose(dose)
  var grad: ParamVec
  ll5(log10(dose), fit.params.toVec, grad)

proc predict*(fit: DrcFit, doses: openArray[float]): seq[float] =
  ## The fitted mean response at each dose.
  for dose in doses:
    result.add fit.predict(dose)

proc tQuantile(fit: DrcFit, level: float): float =
  if not (level > 0.0 and level < 1.0):
    raise newException(ValueError, "confidence level must be in (0, 1), got " & $level)
  studentTQuantile(0.5 + 0.5 * level, float(fit.df))

func waldInterval(fit: DrcFit, param: DrcParam, tq: float): ConfInt =
  let half = tq * fit.se[param]
  (fit.params[param] - half, fit.params[param] + half)

proc confint*(fit: DrcFit, param: DrcParam, level = 0.95): ConfInt =
  ## Wald confidence interval using the t distribution with `fit.df` degrees of freedom.
  fit.waldInterval(param, fit.tQuantile(level))

proc checkPercent(pct: float) =
  if not (pct > 0.0 and pct < 100.0):
    raise newException(ValueError, "percentage must be in (0, 100), got " & $pct)

func logDoseAtFraction(
    fit: DrcFit, r: float
): tuple[x, dxdr: float, grad: array[DrcParam, float]] =
  ## The log10 dose x at which the curve is the fraction r = sᴬ of the way
  ## from Bottom to Top, where s = 1/(1 + 10^z) and z = (LogEC50 − x)·Hill,
  ## with ∂x/∂r and the gradient of x with respect to LogEC50, Hill and
  ## Asymmetry. Solving gives x = LogEC50 − z/Hill with z = log10(r^(−1/A) − 1).
  let
    h = fit.params[Hill]
    a = fit.params[Asymmetry]
    v = pow(r, -1.0 / a)
    z = log10(v - 1.0)
  result.x = fit.params[LogEC50] - z / h
  result.dxdr = v / (a * r * (v - 1.0) * ln(10.0) * h)
  result.grad[LogEC50] = 1.0
  result.grad[Hill] = z / (h * h)
  result.grad[Asymmetry] = -v * ln(r) / (a * a * (v - 1.0) * ln(10.0) * h)

func logEcxAndGradient(fit: DrcFit, pct: float): (float, array[DrcParam, float]) =
  ## log10(ECx) and its gradient with respect to the parameters.
  ##
  ## The fraction of the change from the low-dose plateau is sᴬ for a rising
  ## curve and 1 − sᴬ for a falling one, so r = pct/100 (rising) or
  ## 1 − pct/100 (falling). It does not depend on Bottom or Top.
  let r =
    if fit.params[Hill] > 0.0:
      pct / 100.0
    else:
      1.0 - pct / 100.0
  let at = fit.logDoseAtFraction(r)
  (at.x, at.grad)

func logDoseAtAndGradient(
    fit: DrcFit, response: float
): tuple[x, dxdy: float, grad: array[DrcParam, float]] =
  ## log10 of the dose at which the fitted curve equals `response`, with
  ## ∂x/∂response and the gradient with respect to all the parameters. NaN
  ## unless `response` lies strictly between the plateaus of a curve that
  ## is not flat.
  let
    b = fit.params[Bottom]
    t = fit.params[Top]
    r = (response - b) / (t - b)
  if not (r > 0.0 and r < 1.0) or fit.params[Hill] == 0.0:
    result.x = NaN
    result.dxdy = NaN
    for param in DrcParam:
      result.grad[param] = NaN
    return
  let at = fit.logDoseAtFraction(r)
  # r = (response − Bottom)/(Top − Bottom), so ∂r/∂Bottom = −(1 − r)/(Top − Bottom)
  # and ∂r/∂Top = −r/(Top − Bottom).
  result.x = at.x
  result.dxdy = at.dxdr / (t - b)
  result.grad = at.grad
  result.grad[Bottom] = -at.dxdr * (1.0 - r) / (t - b)
  result.grad[Top] = -at.dxdr * r / (t - b)

proc ecx*(fit: DrcFit, pct: float): float =
  ## Dose giving `pct`% of the change from the low-dose plateau to the
  ## high-dose plateau (drc's relative ED).
  checkPercent(pct)
  pow(10.0, fit.logEcxAndGradient(pct)[0])

proc ec50*(fit: DrcFit): float =
  ## The half-way dose. For 3PL and 4PL fits this is 10^LogEC50.
  fit.ecx(50.0)

func deltaVariance(
    fit: DrcFit, grad: array[DrcParam, float], params: Slice[DrcParam]
): float =
  ## gᵀ·cov·g over `params`, the delta-method variance of a function with gradient g.
  for i in params:
    for j in params:
      result += grad[i] * fit.cov[i][j] * grad[j]

func ecxInterval(fit: DrcFit, pct, tq: float): ConfInt =
  # Delta method on log10(ECx), which does not depend on Bottom or Top.
  let (center, grad) = fit.logEcxAndGradient(pct)
  let half = tq * sqrt(fit.deltaVariance(grad, LogEC50 .. Asymmetry))
  (pow(10.0, center - half), pow(10.0, center + half))

proc ecxCI*(fit: DrcFit, pct: float, level = 0.95): ConfInt =
  ## Confidence interval for `ecx`, by the delta method on log10 dose and
  ## then back-transformed, so it is asymmetric on the dose scale.
  checkPercent(pct)
  fit.ecxInterval(pct, fit.tQuantile(level))

proc ec50CI*(fit: DrcFit, level = 0.95): ConfInt =
  fit.ecxCI(50.0, level)

proc residualVariance(fit: DrcFit, weight: float): float =
  if not (weight > 0.0 and weight < Inf):
    raise newException(ValueError, "weight must be finite and positive, got " & $weight)
  fit.sigma ^ 2 / weight

proc checkResponse(response: float) =
  if not response.isFinite:
    raise newException(ValueError, "response must be finite, got " & $response)

proc doseAt*(fit: DrcFit, response: float): float =
  ## The dose at which the fitted mean response equals `response` (inverse
  ## prediction). For data normalised to % of control, `doseAt(50.0)` is the
  ## absolute IC50, where `ec50` is the relative one (half-way between the
  ## fitted plateaus). NaN unless `response` lies strictly between the
  ## plateaus. Doses outside `logDoseRange` are extrapolated.
  checkResponse(response)
  pow(10.0, fit.logDoseAtAndGradient(response).x)

proc doseAtInterval(fit: DrcFit, response, level, observedVariance: float): ConfInt =
  ## Delta method on log10 dose, back-transformed. `observedVariance` is the
  ## variance of an observed response, 0 for a mean response.
  checkResponse(response)
  let tq = fit.tQuantile(level)
  let (center, dxdy, grad) = fit.logDoseAtAndGradient(response)
  let half =
    tq *
    sqrt(
      fit.deltaVariance(grad, DrcParam.low .. DrcParam.high) +
        dxdy * dxdy * observedVariance
    )
  (pow(10.0, center - half), pow(10.0, center + half))

proc doseAtCI*(fit: DrcFit, response: float, level = 0.95): ConfInt =
  ## Confidence interval for `doseAt`: where the dose lies at which the true
  ## mean response is `response` (an absolute ECx). It uses the delta method
  ## on log10 dose and is back-transformed, as `ecxCI`. For the dose of a
  ## sample whose response was measured, use `calibrationCI`.
  fit.doseAtInterval(response, level, 0.0)

proc calibrationCI*(
    fit: DrcFit, response: float, level = 0.95, weight = 1.0, replicates = 1
): ConfInt =
  ## Interval for the unknown dose of a sample with measured `response` (the
  ## mean of `replicates` observations, each of relative `weight`), read off
  ## the fitted curve as a standard curve. As `doseAtCI`, plus the variance
  ## σ²/(`weight`·`replicates`) of the measured response.
  if replicates < 1:
    raise newException(ValueError, "replicates must be at least 1, got " & $replicates)
  fit.doseAtInterval(response, level, fit.residualVariance(weight * float(replicates)))

proc bandAt(fit: DrcFit, dose, tq, extraVariance: float): Band =
  ## The fitted response at `dose` ± tq·sqrt(Var(f̂) + extraVariance), with
  ## Var(f̂) = gᵀ·cov·g by the delta method.
  checkDose(dose)
  var grad: ParamVec
  let value = ll5(log10(dose), fit.params.toVec, grad)
  let half =
    tq *
    sqrt(
      fit.deltaVariance(grad.toParams, DrcParam.low .. DrcParam.high) + extraVariance
    )
  (value, value - half, value + half)

proc confidenceBand*(fit: DrcFit, doses: openArray[float], level = 0.95): seq[Band] =
  ## Pointwise confidence band for the fitted curve: where the true mean
  ## response lies at each dose.
  let tq = fit.tQuantile(level)
  for dose in doses:
    result.add fit.bandAt(dose, tq, 0.0)

proc confidenceBand*(fit: DrcFit, dose: float, level = 0.95): Band =
  ## The band at one dose. For many doses, the `openArray` overload computes
  ## the t quantile once.
  fit.bandAt(dose, fit.tQuantile(level), 0.0)

proc predictionBand*(
    fit: DrcFit, doses: openArray[float], level = 0.95, weight = 1.0
): seq[Band] =
  ## Pointwise prediction band: where a new observation at each dose is
  ## expected to fall. It adds the residual variance σ²/`weight` to the
  ## confidence band, with `weight` on the same relative scale as the fit's.
  let extra = fit.residualVariance(weight)
  let tq = fit.tQuantile(level)
  for dose in doses:
    result.add fit.bandAt(dose, tq, extra)

proc predictionBand*(fit: DrcFit, dose: float, level = 0.95, weight = 1.0): Band =
  fit.bandAt(dose, fit.tQuantile(level), fit.residualVariance(weight))

proc fitted*(fit: DrcFit): seq[float] =
  ## The fitted mean response at each observation's dose.
  fit.predict(fit.doses)

proc residuals*(fit: DrcFit): seq[float] =
  ## Response residuals yᵢ − f(xᵢ), as R's `residuals(type = "response")`.
  ## Zero-weight observations are included.
  result = fit.fitted
  for i in 0 ..< result.len:
    result[i] = fit.responses[i] - result[i]

proc pearsonResiduals*(fit: DrcFit): seq[float] =
  ## √wᵢ·(yᵢ − f(xᵢ))/σ, as R's `residuals(type = "pearson")`: roughly unit
  ## variance under the fitted weights, so values beyond ±2–3 stand out. For
  ## a robust fit, σ is the robust scale. 0 for zero-weight observations.
  result = fit.residuals
  for i in 0 ..< result.len:
    result[i] = sqrt(fit.weights[i]) * result[i] / fit.sigma

proc pureError(
    doses, responses, weights: openArray[float]
): tuple[ss: float, groups: int] =
  ## The weighted sum of squares about each dose's weighted mean response,
  ## and the number of doses, over observations with non-zero weight.
  var points: seq[(float, int)] # (dose, observation), sorted by dose
  for i in 0 ..< doses.len:
    if weights[i] > 0.0:
      points.add (doses[i], i)
  points.sort()
  var k = 0
  while k < points.len:
    var j = k
    var sumW, sumWY = 0.0
    while j < points.len and points[j][0] == points[k][0]:
      let i = points[j][1]
      sumW += weights[i]
      sumWY += weights[i] * responses[i]
      inc j
    for (_, i) in points[k ..< j]:
      result.ss += weights[i] * (responses[i] - sumWY / sumW) ^ 2
    inc result.groups
    k = j

proc lackOfFitTest(ssr: float, df, n: int, pe: tuple[ss: float, groups: int]): FTest =
  ## F test of the curve against one mean per dose: lack of fit
  ## (SSR − SSpe on groups − p df) against pure error (SSpe on n − groups df).
  let
    df1 = pe.groups - (n - df)
    df2 = n - pe.groups
  if df2 < 1:
    raise newException(
      ValueError, "a lack-of-fit test needs replicate observations at some dose"
    )
  if df1 < 1:
    raise newException(
      ValueError,
      "a lack-of-fit test needs more doses than fitted parameters, got " & $pe.groups &
        " doses for " & $(n - df) & " parameters",
    )
  let f = ((ssr - pe.ss) / float(df1)) / (pe.ss / float(df2))
  (f, df1, df2, fTail(f, float(df1), float(df2)))

proc lackOfFit*(fit: DrcFit): FTest =
  ## Lack-of-fit F test: does the curve describe the mean responses at each
  ## dose? It compares the fit with one free mean per dose, using the
  ## scatter among replicates as pure error, as R's
  ## `anova(fit, lm(y ~ factor(dose)))` and drc's `modelFit`. A small p means
  ## the curve misses the dose means by more than replicate scatter explains.
  ## Needs replicates at some dose and more doses than free parameters. Not
  ## available for robust fits or for the curves of a joint fit (use
  ## `lackOfFit` on the `MultiFit`).
  if fit.loss != Squared:
    raise newException(ValueError, "a lack-of-fit test needs a least-squares fit")
  if fit.jointView:
    raise
      newException(ValueError, "for a curve of a joint fit, test the joint fit instead")
  lackOfFitTest(
    fit.ssr, fit.df, fit.n, pureError(fit.doses, fit.responses, fit.weights)
  )

proc profileTau[F, V, B](
    fn: F,
    x, y, w: openArray[float],
    j: int,
    beta: float,
    start: V,
    pinned: B,
    lower, upper: V,
    ssrHat, sigma: float,
    lm: LMOptions,
    params: var V,
): float =
  ## The profile t statistic τ(β) = sqrt(S(β) − Ŝ)/σ, where S(β) is the SSR
  ## refitted from `start` with parameter j held at β; `params` gets the refit.
  ## Inf if the refit has no finite SSR.
  var p = start
  p[j] = beta
  let r = levenbergMarquardt(fn, x, y, w, p, pinned, lower, upper, lm)
  params = r.params
  if r.ssr.isNaN:
    Inf
  else:
    sqrt(max(r.ssr - ssrHat, 0.0)) / sigma

proc profileEnd[F, V, B](
    fn: F,
    x, y, w: openArray[float],
    pHat: V,
    j: int,
    free: B,
    lower, upper: V,
    ssrHat, sigma, tq, step, direction: float,
    lm: LMOptions,
    asymmetrySlot = ord(Asymmetry),
): float =
  ## One end of a profile-likelihood interval: the value β of parameter j,
  ## beyond pHat[j] in `direction` (±1), at which τ(β) reaches `tq`. It steps
  ## outwards from pHat[j] (by `step`, growing 1.5× each time) until τ passes
  ## `tq`, then finds the crossing by the Illinois method. Each refit starts
  ## from the nearest point inside the crossing, so the profile follows one
  ## valley. If τ stays below `tq` up to a bound, the bound is the end (0 for
  ## Asymmetry, which must be positive: `asymmetrySlot` says where it is, or
  ## −1 if parameter j is not an Asymmetry). As R's `profile.nls`, the profile is
  ## followed only to 10·`tq` steps from the estimate, and if τ has not
  ## reached `tq` by then the end is ±Inf (a finite bound beyond that is
  ## still tried): far out, the refits are poorly conditioned (on a partial
  ## curve, a plateau runs off to infinity) and can stop short of the
  ## minimum, which would look like a crossing.
  var pinned = free
  pinned[j] = false
  let open = j == asymmetrySlot and direction < 0.0 and not (lower[j] > 0.0)
  let limit =
    if direction > 0.0:
      upper[j]
    elif open:
      0.0
    else:
      lower[j]
  let reach = 10.0 * tq * step
  var
    inner = pHat[j]
    innerP = pHat
    innerTau = 0.0
    outer, outerTau: float
    found = false
    h = step
    p = pHat
  for _ in 0 ..< 60:
    var beta = inner + direction * h
    var atLimit = false
    if (beta - limit) * direction >= 0.0:
      if open:
        beta = 0.5 * (inner + limit) # approach 0 without reaching it
      else:
        (beta, atLimit) = (limit, true)
    if not atLimit and abs(beta - pHat[j]) > reach:
      if not limit.isFinite:
        return direction * Inf
      (beta, atLimit) = (limit, true)
    let t = profileTau(
      fn, x, y, w, j, beta, innerP, pinned, lower, upper, ssrHat, sigma, lm, p
    )
    if t >= tq:
      (outer, outerTau, found) = (beta, t, true)
      break
    if atLimit:
      return limit
    (inner, innerP, innerTau) = (beta, p, t)
    h *= 1.5
  if not found:
    return limit # only reachable when approaching Asymmetry's 0

  var
    (a, ga) = (inner, innerTau - tq) # ga < 0
    (b, gb) = (outer, outerTau - tq) # gb ≥ 0
    side = 0
  for _ in 0 ..< 200:
    if abs(b - a) <= 1e-12 * step:
      break
    var c = (a * gb - b * ga) / (gb - ga)
    if not (c > min(a, b) and c < max(a, b)):
      c = 0.5 * (a + b) # e.g. gb = Inf
    let gc =
      profileTau(fn, x, y, w, j, c, innerP, pinned, lower, upper, ssrHat, sigma, lm, p) -
      tq
    if gc == 0.0:
      return c
    if gc < 0.0:
      (a, ga, innerP) = (c, gc, p)
      if side == -1:
        gb *= 0.5
      side = -1
    else:
      (b, gb) = (c, gc)
      if side == 1:
        ga *= 0.5
      side = 1
  0.5 * (a + b)

proc checkProfile(fit: DrcFit) =
  if fit.loss != Squared:
    raise
      newException(ValueError, "profile-likelihood intervals need a least-squares fit")
  if fit.jointView:
    raise newException(
      ValueError,
      "profile-likelihood intervals are not available for curves of a joint fit",
    )

proc profileData(
    fit: DrcFit
): tuple[x: seq[float], free: array[drcParamCount, bool], box: Box] =
  result.x = newSeq[float](fit.doses.len)
  for i, dose in fit.doses:
    result.x[i] = log10(dose)
  for param in DrcParam:
    result.free[ord(param)] = param notin fit.fixed
  result.box = (fit.lower.toVec, fit.upper.toVec)

proc profileStep(fit: DrcFit, param: DrcParam): float =
  ## The first step out from the estimate: its standard error, if usable.
  if fit.se[param] > 0.0 and fit.se[param] < Inf:
    fit.se[param]
  else:
    0.1 * max(abs(fit.params[param]), 1.0)

proc profileCI*(
    fit: DrcFit, param: DrcParam, level = 0.95, opts = FitOptions()
): ConfInt =
  ## Profile-likelihood confidence interval for `param`, as R's
  ## `confint(nls fit)`: the values β at which
  ## τ(β) = sqrt((S(β) − Ŝ)/σ²) equals the t quantile, where S(β) is the SSR
  ## with `param` held at β and the other free parameters refitted. Unlike
  ## the Wald interval (`confint`), it follows the curvature of the SSR, so
  ## it can be asymmetric. An end is a bound if the data do not exclude the
  ## bound (0 for Asymmetry), and ±Inf if the data do not bound the
  ## parameter within 10·t standard errors of the estimate (where R gives NA). A fixed parameter's interval is its value; NaN without
  ## residual degrees of freedom. `opts.lm` controls the refits.
  ## Least-squares fits only, and not for the curves of a joint fit.
  fit.checkProfile
  let tq = fit.tQuantile(level)
  let v = fit.params[param]
  if param in fit.fixed:
    return (v, v)
  if not (fit.sigma > 0.0 and fit.sigma < Inf):
    return (NaN, NaN)
  let (x, free, box) = fit.profileData
  let step = fit.profileStep(param)
  var bounds: array[2, float]
  for k, direction in [-1.0, 1.0]:
    bounds[k] = profileEnd(
      ll5,
      x,
      fit.responses,
      fit.weights,
      fit.params.toVec,
      ord(param),
      free,
      box.lower,
      box.upper,
      fit.ssr,
      fit.sigma,
      tq,
      step,
      direction,
      opts.lm,
    )
  (bounds[0], bounds[1])

proc ecxProfileCI*(
    fit: DrcFit, pct: float, level = 0.95, opts = FitOptions()
): ConfInt =
  ## Profile-likelihood confidence interval for `ecx`: the model is
  ## reparametrised with log10(ECx) in place of LogEC50, and that is profiled
  ## as in `profileCI`, then back-transformed. An unbounded end is 0 or Inf.
  ## Not available when LogEC50 is fixed or bounded.
  checkPercent(pct)
  fit.checkProfile
  let tq = fit.tQuantile(level)
  if LogEC50 in fit.fixed:
    raise newException(ValueError, "LogEC50 is fixed, so ECx cannot be profiled")
  if fit.lower[LogEC50] > -Inf or fit.upper[LogEC50] < Inf:
    raise newException(ValueError, "LogEC50 is bounded, so ECx cannot be profiled")
  if not (fit.sigma > 0.0 and fit.sigma < Inf):
    return (NaN, NaN)
  let r =
    if fit.params[Hill] > 0.0:
      pct / 100.0
    else:
      1.0 - pct / 100.0
  # q is the parameter vector with q[LogEC50] = log10(ECx). From
  # `logDoseAtFraction`, LogEC50 = log10(ECx) + z/Hill, z = log10(r^(−1/A) − 1).
  let model = proc(x: float, q: ParamVec, grad: var ParamVec): float =
    let
      h = q[ord(Hill)]
      a = q[ord(Asymmetry)]
      v = pow(r, -1.0 / a)
      z = log10(v - 1.0)
    var p = q
    p[ord(LogEC50)] = q[ord(LogEC50)] + z / h
    result = ll5(x, p, grad)
    let dl = grad[ord(LogEC50)]
    grad[ord(Hill)] -= dl * z / (h * h)
    grad[ord(Asymmetry)] += dl * v * ln(r) / (a * a * (v - 1.0) * ln(10.0) * h)
  let (x, free, box) = fit.profileData
  var q = fit.params.toVec
  let (logEcx, grad) = fit.logEcxAndGradient(pct)
  q[ord(LogEC50)] = logEcx
  var step = sqrt(fit.deltaVariance(grad, LogEC50 .. Asymmetry))
  if not (step > 0.0 and step < Inf):
    step = 0.1
  var bounds: array[2, float]
  for k, direction in [-1.0, 1.0]:
    bounds[k] = pow(
      10.0,
      profileEnd(
        model,
        x,
        fit.responses,
        fit.weights,
        q,
        ord(LogEC50),
        free,
        box.lower,
        box.upper,
        fit.ssr,
        fit.sigma,
        tq,
        step,
        direction,
        opts.lm,
      ),
    )
  (bounds[0], bounds[1])

proc ec50ProfileCI*(fit: DrcFit, level = 0.95, opts = FitOptions()): ConfInt =
  fit.ecxProfileCI(50.0, level, opts)

func better[R](a, b: R): bool =
  ## A lower SSR wins, and any finite SSR beats NaN. Within rounding of each
  ## other, a converged fit beats an unconverged one; otherwise the earlier
  ## (default) start keeps its place.
  if b.ssr.isNaN:
    not a.ssr.isNaN
  elif a.ssr < b.ssr * (1.0 - 1e-10):
    true
  elif a.ssr <= b.ssr * (1.0 + 1e-10):
    a.status.converged and not b.status.converged
  else:
    false

proc multiStart(
    x, y, w: openArray[float],
    p0: ParamVec,
    free: array[drcParamCount, bool],
    box: Box,
    opts: FitOptions,
): LMResult[drcParamCount] =
  ## Levenberg–Marquardt from `p0` and, with `opts.multiStart`, from starts with
  ## LogEC50 shifted by ±1 and Hill scaled by 0.4 or 2.5; the best fit wins.
  ## Starts outside `box` are clamped into it.
  result = levenbergMarquardt(ll5, x, y, w, p0, free, box.lower, box.upper, opts.lm)
  if not opts.multiStart:
    return
  for shift in [0.0, -1.0, 1.0]:
    for scale in [1.0, 0.4, 2.5]:
      if (shift == 0.0 and scale == 1.0) or (shift != 0.0 and not free[ord(LogEC50)]) or
          (scale != 1.0 and not free[ord(Hill)]):
        continue
      var p = p0
      p[ord(LogEC50)] += shift
      p[ord(Hill)] *= scale
      let r = levenbergMarquardt(ll5, x, y, w, p, free, box.lower, box.upper, opts.lm)
      if r.better(result):
        result = r

proc searchLM(
    x, y, w: openArray[float],
    p0: ParamVec,
    free: array[drcParamCount, bool],
    box: Box,
    opts: FitOptions,
): LMResult[drcParamCount] =
  ## `p0` must lie within `box`, since the symmetric stage of a 5PL holds
  ## Asymmetry at its starting value.
  if not free[ord(Asymmetry)]:
    return multiStart(x, y, w, p0, free, box, opts)
  # 5PL: fit the symmetric curve first, then release Asymmetry from there.
  var symmetric = free
  symmetric[ord(Asymmetry)] = false
  let base = multiStart(x, y, w, p0, symmetric, box, opts)
  result = base
  let released =
    levenbergMarquardt(ll5, x, y, w, base.params, free, box.lower, box.upper, opts.lm)
  if released.better(result):
    result = released
  if opts.multiStart:
    for a in [0.5, 2.0]:
      var p = base.params
      p[ord(Asymmetry)] = a
      let r = levenbergMarquardt(ll5, x, y, w, p, free, box.lower, box.upper, opts.lm)
      if r.better(result):
        result = r

proc standardisedResiduals(x, y, w: openArray[float], p: ParamVec): seq[float] =
  ## √wᵢ·(yᵢ − f(xᵢ)), and 0 for observations with zero weight.
  result = newSeq[float](x.len)
  var grad: ParamVec
  for i in 0 ..< x.len:
    if w[i] > 0.0:
      result[i] = sqrt(w[i]) * (y[i] - ll5(x[i], p, grad))

type IrlsResult = object
  fit: LMResult[drcParamCount]
  scale: float
  residuals: seq[float] ## standardised, at fit.params
  converged: bool

proc irls(
    x, y, w: openArray[float],
    start: LMResult[drcParamCount],
    free: array[drcParamCount, bool],
    box: Box,
    loss: Loss,
    c: float,
    scaleFloor: float,
    tol: float,
    opts: FitOptions,
): IrlsResult =
  ## M-estimation by iteratively reweighted least squares, as R's
  ## robustbase::nlrob. Each round takes the robust scale s of the current
  ## residuals (at least `scaleFloor`), weights observation i by wᵢ·ψ(uᵢ)/uᵢ
  ## with uᵢ = rᵢ/s, and refits. It stops when a converged refit changes the
  ## residuals by less than `tol` (relative).
  result.fit = start
  result.residuals = standardisedResiduals(x, y, w, start.params)
  var
    used: seq[float]
    robustW = newSeq[float](x.len)
  for _ in 1 .. opts.irlsMaxIter:
    used.setLen 0
    for i in 0 ..< x.len:
      if w[i] > 0.0:
        used.add result.residuals[i]
    # Without the floor, a fit through more than half the points to rounding
    # error would get a scale of ~1e-12 and reject every other point.
    result.scale = max(madScale(used), scaleFloor)
    if not (result.scale > 0.0):
      # Only when the responses are all equal and fitted exactly: nothing to reweight.
      result.converged = true
      return
    for i in 0 ..< x.len:
      robustW[i] = w[i] * psiWeight(loss, result.residuals[i] / result.scale, c)
    result.fit = levenbergMarquardt(
      ll5, x, y, robustW, result.fit.params, free, box.lower, box.upper, opts.lm
    )
    if result.fit.status in {lmStalled, lmBadStart}:
      return # these weights cannot be fitted; unchanged residuals would look converged
    let residuals = standardisedResiduals(x, y, w, result.fit.params)
    var change, size = 0.0
    for i in 0 ..< x.len:
      change += (residuals[i] - result.residuals[i]) ^ 2
      size += result.residuals[i] ^ 2
    result.residuals = residuals
    if result.fit.status.converged and sqrt(change / max(size, 1e-20)) <= tol:
      result.converged = true
      return

proc checkObservation(dose, response, weight: float) =
  checkDose(dose)
  if not response.isFinite:
    raise newException(ValueError, "response must be finite, got " & $response)
  if not (weight >= 0.0 and weight < Inf):
    raise
      newException(ValueError, "weight must be finite and non-negative, got " & $weight)

proc resolveFixed(
    model: Model, fixed: openArray[(DrcParam, float)]
): tuple[values: ParamVec, params: set[DrcParam]] =
  ## The fixed parameters and their values, including Asymmetry = 1 for 3PL/4PL
  ## and, for a 3PL that does not give Bottom, Bottom = 0.
  for (param, value) in fixed:
    if param == Asymmetry and model != LL5:
      raise newException(ValueError, "Asymmetry can only be fixed in a 5PL fit")
    if not value.isFinite:
      raise
        newException(ValueError, "fixed " & $param & " must be finite, got " & $value)
    result.values[ord(param)] = value
    result.params.incl param
  if model == LL3 and Bottom notin result.params:
    result.values[ord(Bottom)] = 0.0
    result.params.incl Bottom
  if model != LL5:
    result.values[ord(Asymmetry)] = 1.0
    result.params.incl Asymmetry
  if Asymmetry in result.params and not (result.values[ord(Asymmetry)] > 0.0):
    raise newException(ValueError, "Asymmetry must be positive")

proc resolveBounds(
    model: Model, fixed: set[DrcParam], bounds: openArray[(DrcParam, Slice[float])]
): Box =
  ## The bounds on each parameter, ±Inf where none is given.
  for j in 0 ..< drcParamCount:
    (result.lower[j], result.upper[j]) = (-Inf, Inf)
  for (param, range) in bounds:
    if param == Asymmetry and model != LL5:
      raise newException(ValueError, "Asymmetry can only be bounded in a 5PL fit")
    if param in fixed:
      raise newException(ValueError, $param & " is fixed and cannot also be bounded")
    if range.a.isNaN or range.b.isNaN:
      raise newException(ValueError, "bounds on " & $param & " must not be NaN")
    if not (range.a < range.b):
      raise newException(
        ValueError,
        "the lower bound on " & $param & " must be below the upper one, got " & $range.a &
          " .. " & $range.b & " (to hold it at one value, fix it)",
      )
    (result.lower[ord(param)], result.upper[ord(param)]) = (range.a, range.b)

proc fitCurve*(
    model: Model,
    doses, responses: openArray[float],
    weights: openArray[float] = [],
    fixed: openArray[(DrcParam, float)] = [],
    bounds: openArray[(DrcParam, Slice[float])] = [],
    opts = FitOptions(),
): DrcFit =
  ## Fits `model` to `responses` at `doses` (linear scale; zero doses are
  ## treated as untreated controls). `fitLL3`, `fitLL4` and `fitLL5` are
  ## this with the model given.
  ##
  ## * `weights` are relative: the error variance of observation i is taken
  ##   to be σ²/wᵢ, with σ² estimated from the fit. A weight of 0 excludes
  ##   an observation.
  ## * Parameters in `fixed` are held at the given value. A 3PL fixes Bottom
  ##   at 0 unless `fixed` gives it; 3PL and 4PL fits fix Asymmetry at 1.
  ## * Parameters in `bounds` are kept within their range, e.g.
  ##   `bounds = [(Hill, -10.0 .. -0.1), (Top, -Inf .. 120.0)]`. A parameter
  ##   that ends on a bound is listed in `fit.atBound` and raises `fwAtBound`;
  ##   its standard error is computed as if it were unbounded, as R's
  ##   `nls(algorithm = "port")` does.
  result.model = model
  if doses.len != responses.len:
    raise newException(ValueError, "doses and responses differ in length")
  if weights.len != 0 and weights.len != doses.len:
    raise newException(ValueError, "weights and doses differ in length")

  var
    x, w = newSeq[float](doses.len)
    guessX, guessY: seq[float] # observations with non-zero weight
    sumW, sumWY, sumLogW = 0.0
    lo = Inf
    hi = NegInf
  for i in 0 ..< doses.len:
    x[i] = log10(doses[i])
    w[i] =
      if weights.len == 0:
        1.0
      else:
        weights[i]
    checkObservation(doses[i], responses[i], w[i])
    if w[i] > 0.0:
      inc result.n
      guessX.add x[i]
      guessY.add responses[i]
      sumW += w[i]
      sumWY += w[i] * responses[i]
      sumLogW += ln(w[i])
      if x[i] > NegInf:
        lo = min(lo, x[i])
        hi = max(hi, x[i])
  result.logDoseRange =
    if lo <= hi:
      (lo, hi)
    else:
      (NaN, NaN)
  result.doses = @doses
  result.responses = @responses
  result.weights = w

  let (fixedValues, fixedParams) = resolveFixed(model, fixed)
  result.fixed = fixedParams
  var free: array[drcParamCount, bool]
  for param in DrcParam:
    free[ord(param)] = param notin fixedParams
  let box = resolveBounds(model, fixedParams, bounds)
  for param in DrcParam:
    result.lower[param] = box.lower[ord(param)]
    result.upper[param] = box.upper[ord(param)]

  result.loss = opts.loss
  if opts.loss != Squared:
    result.tuning =
      if opts.tuning.isNaN:
        defaultTuning(opts.loss)
      else:
        opts.tuning
    if not (result.tuning > 0.0 and result.tuning < Inf):
      raise newException(
        ValueError, "tuning must be finite and positive, got " & $opts.tuning
      )
    if opts.irlsMaxIter < 1:
      raise newException(
        ValueError, "irlsMaxIter must be at least 1, got " & $opts.irlsMaxIter
      )
    if not (opts.irlsTol >= 0.0):
      raise
        newException(ValueError, "irlsTol must be non-negative, got " & $opts.irlsTol)

  let nFree = drcParamCount - card(result.fixed)
  if result.n < max(nFree, 1):
    raise newException(
      ValueError,
      "need at least " & $max(nFree, 1) & " observations with " &
        "non-zero weight to fit " & $nFree & " parameters, got " & $result.n,
    )

  var p0 = drcGuess(guessX, guessY, fixedValues, result.fixed)
  for j in 0 ..< drcParamCount:
    p0[j] = clamp(p0[j], box.lower[j], box.upper[j])
  var best = searchLM(x, responses, w, p0, free, box, opts)
  result.converged = best.status.converged
  if opts.fallback and not result.converged:
    # Nelder–Mead needs no derivatives; Levenberg–Marquardt then polishes its result.
    let y = @responses
    let start = if best.ssr.isFinite: best.params else: p0
    let nm = nelderMead(
      proc(p: ParamVec): float =
        weightedSsr(ll5, x, y, w, p),
      start,
      free,
      box.lower,
      box.upper,
      opts.nelderMead,
    )
    # Each stage only goes downhill from the best fit so far, so the polished
    # fit is never worse.
    best =
      levenbergMarquardt(ll5, x, y, w, nm.params, free, box.lower, box.upper, opts.lm)
    result.usedFallback = true
    result.fallbackConverged = best.status.converged or nm.converged
    result.converged = result.fallbackConverged

  let meanY = sumWY / sumW
  var sst = 0.0
  for i in 0 ..< responses.len:
    sst += w[i] * (responses[i] - meanY) ^ 2

  # A robust fit reweights from the least-squares fit.
  var
    robustFit: IrlsResult
    u: seq[float] # residuals in units of the robust scale
  if result.loss != Squared:
    # 1 ppm of the responses' (weighted RMS) spread: far below any real scatter.
    let scaleFloor = 1e-6 * sqrt(sst / float(result.n))
    var start = best
    if result.loss == Tukey:
      # Tukey's redescending ψ can stop in a local minimum, so start from Huber,
      # converged only as far as nlrob's default tolerance.
      start = irls(
        x,
        responses,
        w,
        best,
        free,
        box,
        Huber,
        defaultTuning(Huber),
        scaleFloor,
        max(opts.irlsTol, 1e-6),
        opts,
      ).fit
    robustFit = irls(
      x, responses, w, start, free, box, result.loss, result.tuning, scaleFloor,
      opts.irlsTol, opts,
    )
    best = robustFit.fit
    result.converged = best.status.converged
    if not robustFit.converged:
      result.warnings.incl fwRobustNotConverged
    u = newSeq[float](x.len)
    result.robustWeights = newSeq[float](x.len)
    for i in 0 ..< x.len:
      if w[i] > 0.0:
        if robustFit.scale > 0.0:
          u[i] = robustFit.residuals[i] / robustFit.scale
        result.robustWeights[i] = psiWeight(result.loss, u[i], result.tuning)

  result.params = best.params.toParams
  for param in DrcParam:
    if param notin result.fixed and
        result.params[param] in [result.lower[param], result.upper[param]]:
      result.atBound.incl param
  if result.atBound != {}:
    result.warnings.incl fwAtBound
  var
    jtwj, unscaled: array[drcParamCount, array[drcParamCount, float]]
    gradient: ParamVec
  # Under the user weights, so robust fits report their plain weighted SSR too.
  result.ssr = normalEquations(ll5, x, responses, w, best.params, jtwj, gradient)
  result.iterations = best.iterations
  result.status = best.status
  result.df = result.n - nFree

  if not result.converged:
    result.warnings.incl fwNotConverged
  result.converged = result.warnings * {fwNotConverged, fwRobustNotConverged} == {}
  let logEC50 = log10(result.ec50)
  if not (logEC50 >= result.logDoseRange.lower and logEC50 <= result.logDoseRange.upper):
    result.warnings.incl fwEc50OutsideRange

  # Uncertainty: cov = varianceScale · (JᵀWJ)⁻¹ over the free parameters, with
  # varianceScale = σ² for least squares. Parameters at a bound count as free,
  # as in R's nls(algorithm = "port") and SciPy's curve_fit.
  var varianceScale = NaN
  if result.loss == Squared:
    if result.df > 0:
      result.sigma = sqrt(result.ssr / float(result.df))
    else:
      result.sigma = NaN
    varianceScale = result.sigma ^ 2
  else:
    # Huber's correction, as MASS::rlm: κ²·S / mean(ψ')², where
    # S = s²·Σψ(uᵢ)²/(n − p) and κ = 1 + p·var(ψ')/(n·mean(ψ')²).
    let s = robustFit.scale
    result.sigma = s
    var
      sumPsi2 = 0.0
      derivative: RunningStat # of ψ'(uᵢ)
    for i in 0 ..< x.len:
      if w[i] > 0.0:
        sumPsi2 += psi(result.loss, u[i], result.tuning) ^ 2
        derivative.push(psiDerivative(result.loss, u[i], result.tuning))
    let
      meanD = derivative.mean
      kappa =
        1.0 + float(nFree) * derivative.varianceS / (float(result.n) * meanD * meanD)
    if result.df > 0 and meanD > 0.0:
      varianceScale = (kappa / meanD) ^ 2 * s * s * sumPsi2 / float(result.df)
  if not unscaledCovariance(jtwj, free, unscaled):
    result.warnings.incl fwSingular
  if result.df == 0:
    result.warnings.incl(fwNoResidualDf)
  for i in DrcParam:
    for j in DrcParam:
      # Fixed entries are exactly 0; a NaN σ or (JᵀWJ)⁻¹ carries through to the rest.
      if i in result.fixed or j in result.fixed:
        result.cov[i][j] = 0.0
      else:
        result.cov[i][j] = varianceScale * unscaled[ord(i)][ord(j)]
  for param in DrcParam:
    result.se[param] = sqrt(result.cov[param][param])

  # Goodness of fit.
  result.rSquared = 1.0 - result.ssr / sst
  if result.loss == Squared:
    result.logLik = gaussianLogLik(result.ssr, result.n, sumLogW)
    let k = float(nFree + 1) # free parameters plus σ
    result.aic = -2.0 * result.logLik + 2.0 * k
    result.bic = -2.0 * result.logLik + ln(float(result.n)) * k
  else:
    # A robust fit does not maximise the Gaussian likelihood.
    (result.logLik, result.aic, result.bic) = (NaN, NaN, NaN)

proc fitLL4*(
    doses, responses: openArray[float],
    weights: openArray[float] = [],
    fixed: openArray[(DrcParam, float)] = [],
    bounds: openArray[(DrcParam, Slice[float])] = [],
    opts = FitOptions(),
): DrcFit =
  ## Fits a four-parameter logistic (Hill equation); see `fitCurve`.
  fitCurve(LL4, doses, responses, weights, fixed, bounds, opts)

proc fitLL3*(
    doses, responses: openArray[float],
    weights: openArray[float] = [],
    fixed: openArray[(DrcParam, float)] = [],
    bounds: openArray[(DrcParam, Slice[float])] = [],
    opts = FitOptions(),
): DrcFit =
  ## Fits a three-parameter logistic: the 4PL with Bottom fixed, at 0 unless
  ## `fixed` gives it (e.g. `fixed = [(Bottom, 5.0)]`); see `fitCurve`.
  fitCurve(LL3, doses, responses, weights, fixed, bounds, opts)

proc fitLL5*(
    doses, responses: openArray[float],
    weights: openArray[float] = [],
    fixed: openArray[(DrcParam, float)] = [],
    bounds: openArray[(DrcParam, Slice[float])] = [],
    opts = FitOptions(),
): DrcFit =
  ## Fits a five-parameter logistic, for curves that are not symmetric about
  ## their midpoint; see `fitCurve`. See `ll5` for the parameterisation;
  ## `ec50` and `ecx` give the actual doses for a given response.
  fitCurve(LL5, doses, responses, weights, fixed, bounds, opts)

func boundNote(fit: DrcFit, param: DrcParam): string =
  if param in fit.atBound:
    if fit.params[param] == fit.lower[param]:
      "  (at lower bound)"
    else:
      "  (at upper bound)"
  else:
    ""

proc summary*(fit: DrcFit, level = 0.95): string =
  ## A printable table of estimates, standard errors and confidence intervals.
  let
    pct = &"{100.0 * level:g}%"
    tq = fit.tQuantile(level)
  result.add &"{fit.model} fit to {fit.n} observations"
  if fit.loss != Squared:
    result.add &", robust ({fit.loss}, tuning {fit.tuning:g})"
  result.add "\n\n"
  result.add &"""{"":<10}{"Estimate":>12}{"Std.Error":>12}{"Lower " & pct:>14}{"Upper " & pct:>14}""" &
    "\n"
  for param in DrcParam:
    if param == Asymmetry and fit.model != LL5:
      continue
    if param in fit.fixed:
      result.add &"{$param:<10}{fit.params[param]:>12.6g}{\"(fixed)\":>12}\n"
    else:
      let ci = fit.waldInterval(param, tq)
      result.add &"{$param:<10}{fit.params[param]:>12.6g}{fit.se[param]:>12.6g}{ci.lower:>14.6g}{ci.upper:>14.6g}"
      result.add fit.boundNote(param) & "\n"
  let ci = fit.ecxInterval(50.0, tq)
  result.add &"{\"EC50\":<10}{fit.ec50:>12.6g}{\"\":>12}{ci.lower:>14.6g}{ci.upper:>14.6g}\n\n"
  if fit.loss == Squared:
    result.add &"Residual standard error: {fit.sigma:.6g} on {fit.df} degrees of freedom\n"
    result.add &"R²: {fit.rSquared:.6g}   AIC: {fit.aic:.6g}   BIC: {fit.bic:.6g}\n"
  else:
    var strong, rejected = 0
    for rw in fit.robustWeights:
      if rw == 0.0:
        inc rejected
      elif rw < 0.5:
        inc strong
    # Observations the user excluded (weight 0) are not counted as rejected.
    rejected -= fit.robustWeights.len - fit.n
    result.add &"Robust residual scale (MAD): {fit.sigma:.6g} on {fit.df} degrees of freedom\n"
    result.add &"Robust weight below 0.5: {strong} observations; 0 (rejected): {rejected}\n"
    result.add &"R²: {fit.rSquared:.6g}\n"
  if fit.usedFallback:
    let outcome = if fit.fallbackConverged: "converged" else: "did not converge"
    result.add "Levenberg–Marquardt did not converge; the Nelder–Mead fallback " &
      outcome
    if fit.loss == Squared:
      result.add &" (final polishing run: {fit.status} after {fit.iterations} iterations)"
    else:
      result.add &". Last reweighting round: {fit.status} after {fit.iterations} iterations"
  else:
    result.add &"Stopped after {fit.iterations} iterations: {fit.status}"
  for warning in fit.warnings:
    result.add &"\nWarning: {warning}"

include drc/multi

include drc/bootstrap
