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

func logDoses(doses: openArray[float]): seq[float] =
  ## x = log10(dose) for each dose (−Inf for a zero dose).
  result = newSeq[float](doses.len)
  for i, dose in doses:
    result[i] = log10(dose)

proc predict*(fit: DrcFit, dose: float): float =
  ## The fitted mean response at `dose` (0 gives the low-dose plateau).
  checkDose(dose)
  var grad: ParamVec
  ll5(log10(dose), fit.params.toVec, grad)

proc predict*(fit: DrcFit, doses: openArray[float]): seq[float] =
  ## The fitted mean response at each dose.
  for dose in doses:
    result.add fit.predict(dose)

proc checkLevel(level: float) =
  if not (level > 0.0 and level < 1.0):
    raise newException(ValueError, "confidence level must be in (0, 1), got " & $level)

proc tQuantile(fit: DrcFit, level: float): float =
  checkLevel(level)
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

func fractionAt(hill, pct: float): float =
  ## The fraction sᴬ of the way from Bottom to Top at which a curve has made
  ## `pct`% of its change from the low-dose plateau: pct/100 for a rising
  ## curve and 1 − pct/100 for a falling one.
  if hill > 0.0:
    pct / 100.0
  else:
    1.0 - pct / 100.0

func zOffset(r, a: float): tuple[z, dzdr, dzda: float] =
  ## z = log10(r^(−1/A) − 1), with its derivatives. The curve is the fraction
  ## r = sᴬ of the way from Bottom to Top, where s = 1/(1 + 10^z), at
  ## z = (LogEC50 − x)·Hill, so at x = LogEC50 − z/Hill.
  let
    v = pow(r, -1.0 / a)
    d = (v - 1.0) * ln(10.0)
  (log10(v - 1.0), -v / (a * r * d), v * ln(r) / (a * a * d))

func logDoseAtFraction(
    fit: DrcFit, r: float
): tuple[x, dxdr: float, grad: array[DrcParam, float]] =
  ## The log10 dose x at which the curve is the fraction r = sᴬ of the way
  ## from Bottom to Top (see `zOffset`), with ∂x/∂r and the gradient of x
  ## with respect to LogEC50, Hill and Asymmetry.
  let
    h = fit.params[Hill]
    (z, dzdr, dzda) = zOffset(r, fit.params[Asymmetry])
  result.x = fit.params[LogEC50] - z / h
  result.dxdr = -dzdr / h
  result.grad[LogEC50] = 1.0
  result.grad[Hill] = z / (h * h)
  result.grad[Asymmetry] = -dzda / h

func logEcxAndGradient(fit: DrcFit, pct: float): (float, array[DrcParam, float]) =
  ## log10(ECx) and its gradient with respect to the parameters. It does not
  ## depend on Bottom or Top.
  let at = fit.logDoseAtFraction(fractionAt(fit.params[Hill], pct))
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
    result.grad.fill(NaN)
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
  ## a robust fit, σ is the robust scale. 0 for zero-weight observations;
  ## NaN for the others when σ is 0 (an exact fit) or unavailable (no
  ## residual degrees of freedom).
  result = fit.residuals
  let usable = fit.sigma > 0.0 and fit.sigma < Inf
  for i in 0 ..< result.len:
    result[i] =
      if fit.weights[i] == 0.0:
        0.0
      elif usable:
        sqrt(fit.weights[i]) * result[i] / fit.sigma
      else:
        NaN

proc pureError(
    doses, responses, weights: openArray[float]
): tuple[ss: float, groups: int] =
  ## The weighted sum of squares about each dose's weighted mean response,
  ## and the number of doses, over observations with non-zero weight. A dose
  ## whose responses are all equal adds exactly 0 (not the rounding error of
  ## its mean).
  var points: seq[(float, int)] # (dose, observation), sorted by dose
  for i in 0 ..< doses.len:
    if weights[i] > 0.0:
      points.add (doses[i], i)
  points.sort()
  var k = 0
  while k < points.len:
    var j = k
    var sumW, sumWY = 0.0
    var equal = true
    while j < points.len and points[j][0] == points[k][0]:
      let i = points[j][1]
      sumW += weights[i]
      sumWY += weights[i] * responses[i]
      equal = equal and responses[i] == responses[points[k][1]]
      inc j
    if not equal:
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
  if pe.ss == 0.0:
    raise newException(
      ValueError,
      "a lack-of-fit test needs scatter among replicates, but every dose's " &
        "replicates are identical",
    )
  let f = ((ssr - pe.ss) / float(df1)) / (pe.ss / float(df2))
  (f, df1, df2, fTail(f, float(df1), float(df2)))

proc lackOfFit*(fit: DrcFit): FTest =
  ## Lack-of-fit F test: does the curve describe the mean responses at each
  ## dose? It compares the fit with one free mean per dose, using the
  ## scatter among replicates as pure error, as R's
  ## `anova(fit, lm(y ~ factor(dose)))` and drc's `modelFit`. A small p means
  ## the curve misses the dose means by more than replicate scatter explains.
  ## Needs replicates that differ at some dose, and more doses than free
  ## parameters. Not available for robust fits or for the curves of a joint
  ## fit (use `lackOfFit` on the `MultiFit`).
  if fit.loss != Squared:
    raise newException(ValueError, "a lack-of-fit test needs a least-squares fit")
  if fit.jointView:
    raise
      newException(ValueError, "for a curve of a joint fit, test the joint fit instead")
  lackOfFitTest(
    fit.ssr, fit.df, fit.n, pureError(fit.doses, fit.responses, fit.weights)
  )

type ProfileTau[V] = proc(beta: float, start: V, params: var V): float
  ## τ(β) = sqrt(S(β) − Ŝ)/σ, where S(β) is the SSR refitted from `start`
  ## with the profiled quantity held at β; `params` gets the refit. Inf where
  ## there is no such fit.

proc profileTauOf(s, ssrHat, sigma: float): float =
  if s.isNaN:
    Inf
  else:
    sqrt(max(s - ssrHat, 0.0)) / sigma

proc pinnedTau[F, V, B](
    fn: F,
    x, y, w: seq[float],
    j: int,
    free: B,
    lower, upper: V,
    ssrHat, sigma: float,
    lm: LMOptions,
): ProfileTau[V] =
  ## τ for parameter j: each refit holds it at β.
  var pinned = free
  pinned[j] = false
  result = proc(beta: float, start: V, params: var V): float =
    var p = start
    p[j] = beta
    let r = levenbergMarquardt(fn, x, y, w, p, pinned, lower, upper, lm)
    params = r.params
    profileTauOf(r.ssr, ssrHat, sigma)

proc constrainedTau[F, G, V, B](
    fn: F,
    g: G,
    x, y, w: seq[float],
    free: B,
    lower, upper: V,
    ssrHat, sigma, scale: float,
    lm: LMOptions,
): ProfileTau[V] =
  ## τ for a function g(p) of the parameters (`g(p, grad)` returns it and
  ## writes its gradient), with standard error about `scale`: each refit
  ## minimises the SSR subject to g(p) = β, within the bounds. The constraint
  ## is imposed by the augmented Lagrangian method, as one extra observation
  ## (x = NaN) with response β − λ/μ, value g(p) and weight μ, refitting by
  ## Levenberg–Marquardt and updating λ ← λ + μ·(g − β) until g = β. μ makes
  ## the extra term ~1000× as curved as the SSR along g, so each round cuts
  ## the violation ~1000-fold. Inf if no fit within the bounds has g = β.
  let
    mu = 1e3 * (sigma / scale) ^ 2
    exact = 1e-10 * scale # stop here
    feasible = 1e-6 * scale # accept up to here (Levenberg–Marquardt's own precision)
  var xa = x
  xa.add NaN
  var wa = w
  wa.add mu
  let augmented = proc(xi: float, p: V, grad: var V): float =
    if xi.isNaN:
      g(p, grad)
    else:
      fn(xi, p, grad)
  result = proc(beta: float, start: V, params: var V): float =
    var
      ya = y
      p = start
      lambda = 0.0
      violation = Inf
      grad: V
    ya.add 0.0
    grad.ensureLen p.len
    for _ in 0 ..< 30:
      ya[^1] = beta - lambda / mu
      p = levenbergMarquardt(augmented, xa, ya, wa, p, free, lower, upper, lm).params
      let v = g(p, grad) - beta
      if not v.isFinite:
        return Inf
      let previous = violation
      violation = v
      lambda += mu * v
      if abs(v) <= exact or (abs(v) <= feasible and abs(v) > 0.1 * abs(previous)):
        break # converged, or at the limit of the refits' precision
    if not (abs(violation) <= feasible):
      return Inf
    params = p
    # S at g = β exactly, to first order: the profile's slope is −2λ.
    let s = weightedSsr(fn, x, y, w, p) + 2.0 * lambda * violation
    profileTauOf(s, ssrHat, sigma)

proc profileEnd[V](
    tau: ProfileTau[V],
    pHat: V,
    center, lower, upper, tq, step, direction: float,
    positive: bool,
): float =
  ## One end of a profile-likelihood interval: the value β, beyond `center`
  ## (the estimate) in `direction` (±1), at which τ(β) reaches `tq`. It steps
  ## outwards from `center` (by `step`, growing 1.5× each time) until τ passes
  ## `tq`, then finds the crossing by the Illinois method. Each refit starts
  ## from the nearest point inside the crossing (`pHat` at first), so the
  ## profile follows one valley. If τ stays below `tq` up to `lower` or
  ## `upper`, that bound is the end (0 for a `positive` quantity, i.e.
  ## Asymmetry, which is approached but never evaluated). As R's
  ## `profile.nls`, the profile is followed only to 10·`tq` steps from the
  ## estimate: far out, the refits are poorly conditioned (on a partial
  ## curve, a plateau runs off to infinity) and can stop short of the
  ## minimum, which would look like a crossing. If τ has not reached `tq` by
  ## then, the data are taken not to bound the quantity, and the end is the
  ## bound (±Inf if none) without a refit there.
  let open = positive and direction < 0.0 and not (lower > 0.0)
  let limit =
    if direction > 0.0:
      upper
    elif open:
      0.0
    else:
      lower
  let reach = 10.0 * tq * step
  var
    inner = center
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
    if abs(beta - center) > reach:
      return limit
    let t = tau(beta, innerP, p)
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
    let gc = tau(c, innerP, p) - tq
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

func profileStep(se, center: float): float =
  ## The first step out from the estimate: its standard error `se`, or
  ## 0.1·max(|estimate|, 1) if that is unusable.
  if se > 0.0 and se < Inf:
    se
  else:
    0.1 * max(abs(center), 1.0)

proc profileInterval[V](
    tau: ProfileTau[V],
    pHat: V,
    center, lower, upper, sigma, tq, se: float,
    positive = false,
): ConfInt =
  ## Both ends of the profile-likelihood interval (see `profileEnd`) of a
  ## quantity estimated at `center`, with standard error `se`, within
  ## `lower` .. `upper`. NaN without a usable σ.
  if not (sigma > 0.0 and sigma < Inf):
    return (NaN, NaN)
  let step = profileStep(se, center)
  for k, direction in [-1.0, 1.0]:
    let bound =
      profileEnd(tau, pHat, center, lower, upper, tq, step, direction, positive)
    if k == 0:
      result.lower = bound
    else:
      result.upper = bound

proc checkProfile(fit: DrcFit) =
  if fit.loss != Squared:
    raise
      newException(ValueError, "profile-likelihood intervals need a least-squares fit")
  if fit.jointView:
    raise newException(
      ValueError,
      "for a curve of a joint fit, use the MultiFit overloads " &
        "(e.g. `profileCI(fit, curve, param)`)",
    )

proc profileData(
    fit: DrcFit
): tuple[x: seq[float], free: array[drcParamCount, bool], box: Box] =
  result.x = logDoses(fit.doses)
  for param in DrcParam:
    result.free[ord(param)] = param notin fit.fixed
  result.box = (fit.lower.toVec, fit.upper.toVec)

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
  ## parameter within 10·t standard errors of the estimate (where R gives
  ## NA). A fixed parameter's interval is its value; NaN without residual
  ## degrees of freedom. `opts.lm` controls the refits. Least-squares fits
  ## only; for the curves of a joint fit, use `profileCI(MultiFit, …)`.
  fit.checkProfile
  let tq = fit.tQuantile(level)
  if param in fit.fixed:
    return (fit.params[param], fit.params[param])
  let (x, free, box) = fit.profileData
  let j = ord(param)
  profileInterval(
    pinnedTau(
      ll5, x, fit.responses, fit.weights, j, free, box.lower, box.upper, fit.ssr,
      fit.sigma, opts.lm,
    ),
    fit.params.toVec,
    fit.params[param],
    box.lower[j],
    box.upper[j],
    fit.sigma,
    tq,
    fit.se[param],
    param == Asymmetry,
  )

proc logEcxConstraint[V](
    slots: array[DrcParam, int], r: float
): proc(p: V, grad: var V): float =
  ## g(p) = the log10 dose at which the curve whose parameters are at `slots`
  ## of p is the fraction r = sᴬ of the way from Bottom to Top (see
  ## `zOffset`), with its gradient.
  result = proc(p: V, grad: var V): float =
    let
      h = p[slots[Hill]]
      (z, _, dzda) = zOffset(r, p[slots[Asymmetry]])
    grad.setZero()
    grad[slots[LogEC50]] += 1.0
    grad[slots[Hill]] += z / (h * h)
    grad[slots[Asymmetry]] -= dzda / h
    p[slots[LogEC50]] - z / h

proc constrainedInterval[F, G, V, B](
    fn: F,
    g: G,
    x, y, w: seq[float],
    pHat: V,
    free: B,
    lower, upper: V,
    ssrHat, sigma, tq, se: float,
    lm: LMOptions,
): ConfInt =
  ## The profile-likelihood interval for g(p) (see `constrainedTau`), or the
  ## point g(p̂) if g does not depend on any free parameter.
  var grad: V
  grad.ensureLen pHat.len
  let center = g(pHat, grad)
  var depends = false
  for j in 0 ..< pHat.len:
    depends = depends or (free[j] and grad[j] != 0.0)
  if not depends:
    return (center, center)
  if not (sigma > 0.0 and sigma < Inf):
    return (NaN, NaN)
  let tau = constrainedTau(
    fn, g, x, y, w, free, lower, upper, ssrHat, sigma, profileStep(se, center), lm
  )
  profileInterval(tau, pHat, center, -Inf, Inf, sigma, tq, se)

proc ecxProfileCI*(
    fit: DrcFit, pct: float, level = 0.95, opts = FitOptions()
): ConfInt =
  ## Profile-likelihood confidence interval for `ecx`: log10(ECx) is
  ## profiled as in `profileCI`, then back-transformed. An unbounded end is 0
  ## or Inf. The model is reparametrised with log10(ECx) in place of LogEC50
  ## or, if LogEC50 is fixed or bounded, each refit is constrained to the
  ## ECx. If ECx depends on no free parameter (the EC50 of a 4PL with
  ## LogEC50 fixed), the interval is the estimate.
  checkPercent(pct)
  fit.checkProfile
  let tq = fit.tQuantile(level)
  let r = fractionAt(fit.params[Hill], pct)
  let (x, free, box) = fit.profileData
  let (logEcx, grad) = fit.logEcxAndGradient(pct)
  let se = sqrt(fit.deltaVariance(grad, LogEC50 .. Asymmetry))
  if LogEC50 in fit.fixed or fit.lower[LogEC50] > -Inf or fit.upper[LogEC50] < Inf:
    var slots: array[DrcParam, int]
    for param in DrcParam:
      slots[param] = ord(param)
    let ci = constrainedInterval(
      ll5,
      logEcxConstraint[ParamVec](slots, r),
      x,
      fit.responses,
      fit.weights,
      fit.params.toVec,
      free,
      box.lower,
      box.upper,
      fit.ssr,
      fit.sigma,
      tq,
      se,
      opts.lm,
    )
    return (pow(10.0, ci.lower), pow(10.0, ci.upper))
  # q is the parameter vector with q[LogEC50] = log10(ECx), so
  # LogEC50 = log10(ECx) + z/Hill (see `zOffset`). z depends only on
  # Asymmetry, so it is kept between calls.
  var cachedA = NaN
  var z, dzda: float
  let model = proc(x: float, q: ParamVec, grad: var ParamVec): float =
    let (h, a) = (q[ord(Hill)], q[ord(Asymmetry)])
    if a != cachedA:
      let o = zOffset(r, a)
      (z, dzda, cachedA) = (o.z, o.dzda, a)
    var p = q
    p[ord(LogEC50)] = q[ord(LogEC50)] + z / h
    result = ll5(x, p, grad)
    let dl = grad[ord(LogEC50)]
    grad[ord(Hill)] -= dl * z / (h * h)
    grad[ord(Asymmetry)] += dl * dzda / h
  var q = fit.params.toVec
  q[ord(LogEC50)] = logEcx
  let j = ord(LogEC50)
  let ci = profileInterval(
    pinnedTau(
      model, x, fit.responses, fit.weights, j, free, box.lower, box.upper, fit.ssr,
      fit.sigma, opts.lm,
    ),
    q,
    logEcx,
    -Inf,
    Inf,
    fit.sigma,
    tq,
    se,
  )
  (pow(10.0, ci.lower), pow(10.0, ci.upper))

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

proc standardisedResiduals[F, V](fn: F, x, y, w: openArray[float], p: V): seq[float] =
  ## √wᵢ·(yᵢ − f(xᵢ)), and 0 for observations with zero weight.
  result = newSeq[float](x.len)
  var grad: V
  grad.ensureLen p.len
  for i in 0 ..< x.len:
    if w[i] > 0.0:
      result[i] = sqrt(w[i]) * (y[i] - fn(x[i], p, grad))

type IrlsResult[R] = object
  fit: R ## an `LMResult` (one curve) or `LMDynResult` (a joint fit)
  scale: float
  residuals: seq[float] ## standardised, at fit.params
  converged: bool

proc irls[F, V, B, R](
    fn: F,
    x, y, w: openArray[float],
    start: R,
    free: B,
    lower, upper: V,
    loss: Loss,
    c: float,
    scaleFloor: float,
    tol: float,
    opts: FitOptions,
): IrlsResult[R] =
  ## M-estimation by iteratively reweighted least squares, as R's
  ## robustbase::nlrob. Each round takes the robust scale s of the current
  ## residuals (at least `scaleFloor`), weights observation i by wᵢ·ψ(uᵢ)/uᵢ
  ## with uᵢ = rᵢ/s, and refits. It stops when a converged refit changes the
  ## residuals by less than `tol` (relative).
  result.fit = start
  result.residuals = standardisedResiduals(fn, x, y, w, start.params)
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
      fn, x, y, robustW, result.fit.params, free, lower, upper, opts.lm
    )
    if result.fit.status in {lmStalled, lmBadStart}:
      return # these weights cannot be fitted; unchanged residuals would look converged
    let residuals = standardisedResiduals(fn, x, y, w, result.fit.params)
    var change, size = 0.0
    for i in 0 ..< x.len:
      change += (residuals[i] - result.residuals[i]) ^ 2
      size += result.residuals[i] ^ 2
    result.residuals = residuals
    if result.fit.status.converged and sqrt(change / max(size, 1e-20)) <= tol:
      result.converged = true
      return

proc robustFit[F, V, B, R](
    fn: F,
    x, y, w: openArray[float],
    best: R,
    free: B,
    lower, upper: V,
    loss: Loss,
    tuning, sst: float,
    n: int,
    opts: FitOptions,
): IrlsResult[R] =
  ## The robust fit by IRLS (`irls`) from the least-squares fit `best`. `sst`
  ## is the responses' weighted total sum of squares over `n` observations.
  # 1 ppm of the responses' (weighted RMS) spread: far below any real scatter.
  let scaleFloor = 1e-6 * sqrt(sst / float(n))
  var start = best
  if loss == Tukey:
    # Tukey's redescending ψ can stop in a local minimum, so start from Huber,
    # converged only as far as nlrob's default tolerance.
    start = irls(
      fn,
      x,
      y,
      w,
      best,
      free,
      lower,
      upper,
      Huber,
      defaultTuning(Huber),
      scaleFloor,
      max(opts.irlsTol, 1e-6),
      opts,
    ).fit
  irls(
    fn, x, y, w, start, free, lower, upper, loss, tuning, scaleFloor, opts.irlsTol, opts
  )

proc robustUnits[R](
    r: IrlsResult[R], w: openArray[float], loss: Loss, tuning: float
): tuple[u, weights: seq[float]] =
  ## Each observation's residual in units of the robust scale, and its final
  ## robust weight ψ(u)/u (0 where the user's weight is 0).
  result.u = newSeq[float](w.len)
  result.weights = newSeq[float](w.len)
  for i in 0 ..< w.len:
    if w[i] > 0.0:
      if r.scale > 0.0:
        result.u[i] = r.residuals[i] / r.scale
      result.weights[i] = psiWeight(loss, result.u[i], tuning)

proc robustVarianceScale(
    loss: Loss, tuning: float, u, w: openArray[float], scale: float, nFree, n, df: int
): float =
  ## The factor that turns (JᵀWJ)⁻¹ into the covariance of a robust fit:
  ## Huber's correction, as MASS::rlm, κ²·S / mean(ψ')², where
  ## S = s²·Σψ(uᵢ)²/(n − p) and κ = 1 + p·var(ψ')/(n·mean(ψ')²). NaN without
  ## residual degrees of freedom.
  var
    sumPsi2 = 0.0
    derivative: RunningStat # of ψ'(uᵢ)
  for i in 0 ..< u.len:
    if w[i] > 0.0:
      sumPsi2 += psi(loss, u[i], tuning) ^ 2
      derivative.push(psiDerivative(loss, u[i], tuning))
  let
    meanD = derivative.mean
    kappa = 1.0 + float(nFree) * derivative.varianceS / (float(n) * meanD * meanD)
  if df > 0 and meanD > 0.0:
    (kappa / meanD) ^ 2 * scale * scale * sumPsi2 / float(df)
  else:
    NaN

proc resolveTuning(opts: FitOptions): float =
  ## The robust loss's tuning constant (0 for least squares, as before a
  ## robust fit sets it), after checking the robust options.
  if opts.loss == Squared:
    return 0.0
  result =
    if opts.tuning.isNaN:
      defaultTuning(opts.loss)
    else:
      opts.tuning
  if not (result > 0.0 and result < Inf):
    raise newException(
      ValueError, "tuning must be finite and positive, got " & $opts.tuning
    )
  if opts.irlsMaxIter < 1:
    raise newException(
      ValueError, "irlsMaxIter must be at least 1, got " & $opts.irlsMaxIter
    )
  if not (opts.irlsTol >= 0.0):
    raise newException(ValueError, "irlsTol must be non-negative, got " & $opts.irlsTol)

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

func freeMask(fixed: set[DrcParam]): array[drcParamCount, bool] =
  ## True for each parameter not in `fixed`.
  for param in DrcParam:
    result[ord(param)] = param notin fixed

func nonZeroLogDoseRange(doses, weights: openArray[float]): tuple[lower, upper: float] =
  ## The smallest and largest log10 of the non-zero doses with non-zero
  ## weight, or NaN if there are none.
  var
    lo = Inf
    hi = NegInf
  for i, dose in doses:
    let x = log10(dose)
    if weights[i] > 0.0 and x > NegInf:
      lo = min(lo, x)
      hi = max(hi, x)
  if lo <= hi:
    (lo, hi)
  else:
    (NaN, NaN)

func totalSumOfSquares(y, w: openArray[float]): float =
  ## Σ wᵢ·(yᵢ − ȳ)² about the weighted mean ȳ.
  var sumW, sumWY = 0.0
  for i in 0 ..< y.len:
    sumW += w[i]
    sumWY += w[i] * y[i]
  let mean = sumWY / sumW
  for i in 0 ..< y.len:
    result += w[i] * (y[i] - mean) ^ 2

proc likelihoodCriteria(
    ssr: float, n, nFree: int, weights: openArray[float]
): tuple[logLik, aic, bic: float] =
  ## The Gaussian log-likelihood of a least-squares fit with `nFree` free
  ## parameters to `n` observations with non-zero weight, and its AIC and BIC,
  ## which count σ as a parameter (as R and drc).
  var sumLogW = 0.0
  for w in weights:
    if w > 0.0:
      sumLogW += ln(w)
  result.logLik = gaussianLogLik(ssr, n, sumLogW)
  let k = float(nFree + 1)
  result.aic = -2.0 * result.logLik + 2.0 * k
  result.bic = -2.0 * result.logLik + ln(float(n)) * k

proc setData(fit: var DrcFit, doses, responses, weights: openArray[float]) =
  ## Checks the observations and stores them in `fit`, with weights of 1 if
  ## none are given, along with `n` and `logDoseRange`.
  if doses.len != responses.len:
    raise newException(ValueError, "doses and responses differ in length")
  if weights.len != 0 and weights.len != doses.len:
    raise newException(ValueError, "weights and doses differ in length")
  fit.doses = @doses
  fit.responses = @responses
  fit.weights = newSeq[float](doses.len)
  for i in 0 ..< doses.len:
    fit.weights[i] =
      if weights.len == 0:
        1.0
      else:
        weights[i]
    checkObservation(doses[i], responses[i], fit.weights[i])
    if fit.weights[i] > 0.0:
      inc fit.n
  fit.logDoseRange = nonZeroLogDoseRange(fit.doses, fit.weights)

proc guessStart(
    x, y, w: openArray[float], fixedValues: ParamVec, fixed: set[DrcParam]
): ParamVec =
  ## `drcGuess` from the observations with non-zero weight.
  var guessX, guessY: seq[float]
  for i in 0 ..< x.len:
    if w[i] > 0.0:
      guessX.add x[i]
      guessY.add y[i]
  drcGuess(guessX, guessY, fixedValues, fixed)

proc fallback(
    x, y, w: seq[float],
    start: ParamVec,
    free: array[drcParamCount, bool],
    box: Box,
    opts: FitOptions,
): tuple[fit: LMResult[drcParamCount], converged: bool] =
  ## Nelder–Mead from `start`, which needs no derivatives, then
  ## Levenberg–Marquardt to polish its result. Each stage only goes downhill,
  ## so from the best fit so far the polished fit is never worse. It has
  ## converged if either stage has.
  let nm = nelderMead(
    proc(p: ParamVec): float =
      weightedSsr(ll5, x, y, w, p),
    start,
    free,
    box.lower,
    box.upper,
    opts.nelderMead,
  )
  result.fit =
    levenbergMarquardt(ll5, x, y, w, nm.params, free, box.lower, box.upper, opts.lm)
  result.converged = result.fit.status.converged or nm.converged

proc flagBounds(fit: var DrcFit) =
  ## `atBound`, and its warning: the free parameters that are on a bound.
  for param in DrcParam:
    if param notin fit.fixed and
        fit.params[param] in [fit.lower[param], fit.upper[param]]:
      fit.atBound.incl param
  if fit.atBound != {}:
    fit.warnings.incl fwAtBound

proc flagEc50Range(fit: var DrcFit) =
  ## Warns if the EC50 is outside the tested doses.
  let logEC50 = log10(fit.ec50)
  if not (logEC50 >= fit.logDoseRange.lower and logEC50 <= fit.logDoseRange.upper):
    fit.warnings.incl fwEc50OutsideRange

proc setUncertainty(fit: var DrcFit, x: seq[float], robustScale: float, u: seq[float]) =
  ## `ssr`, `df`, `sigma`, `cov` and `se` at the final estimates, independent
  ## of which optimiser produced them: cov = varianceScale · (JᵀWJ)⁻¹ over the
  ## free parameters, with varianceScale = σ² for least squares. Parameters at
  ## a bound count as free, as in R's nls(algorithm = "port") and SciPy's
  ## curve_fit. For a robust fit, `robustScale` is its scale and `u` the
  ## residuals in units of it.
  let
    free = freeMask(fit.fixed)
    nFree = drcParamCount - card(fit.fixed)
  var
    jtwj, unscaled: array[drcParamCount, array[drcParamCount, float]]
    gradient: ParamVec
  # Under the user weights, so robust fits report their plain weighted SSR too.
  fit.ssr = normalEquations(
    ll5, x, fit.responses, fit.weights, fit.params.toVec, jtwj, gradient
  )
  fit.df = fit.n - nFree
  var varianceScale = NaN
  if fit.loss == Squared:
    fit.sigma =
      if fit.df > 0:
        sqrt(fit.ssr / float(fit.df))
      else:
        NaN
    varianceScale = fit.sigma ^ 2
  else:
    fit.sigma = robustScale
    varianceScale = robustVarianceScale(
      fit.loss, fit.tuning, u, fit.weights, robustScale, nFree, fit.n, fit.df
    )
  if not unscaledCovariance(jtwj, free, unscaled):
    fit.warnings.incl fwSingular
  if fit.df == 0:
    fit.warnings.incl fwNoResidualDf
  for i in DrcParam:
    for j in DrcParam:
      # Fixed entries are exactly 0; a NaN σ or (JᵀWJ)⁻¹ carries through to the rest.
      if i in fit.fixed or j in fit.fixed:
        fit.cov[i][j] = 0.0
      else:
        fit.cov[i][j] = varianceScale * unscaled[ord(i)][ord(j)]
  for param in DrcParam:
    fit.se[param] = sqrt(fit.cov[param][param])

proc setGoodnessOfFit(fit: var DrcFit, sst: float) =
  ## `rSquared` from `ssr` and the total sum of squares `sst`, and the
  ## likelihood criteria (NaN for a robust fit, which does not maximise the
  ## Gaussian likelihood).
  fit.rSquared = 1.0 - fit.ssr / sst
  (fit.logLik, fit.aic, fit.bic) =
    if fit.loss == Squared:
      likelihoodCriteria(fit.ssr, fit.n, drcParamCount - card(fit.fixed), fit.weights)
    else:
      (NaN, NaN, NaN)

proc fitCurveFrom(
    model: Model,
    doses, responses: openArray[float],
    weights: openArray[float],
    fixed: openArray[(DrcParam, float)],
    bounds: openArray[(DrcParam, Slice[float])],
    opts: FitOptions,
    start: ParamVec,
    warm: bool,
): DrcFit =
  ## `fitCurve`, or with `warm`, one Levenberg–Marquardt run from `start`
  ## in place of the guess and multi-start search (the fallback still
  ## applies); for refits of data close to an earlier fit's.
  result.model = model
  result.setData(doses, responses, weights)
  let (fixedValues, fixedParams) = resolveFixed(model, fixed)
  result.fixed = fixedParams
  let box = resolveBounds(model, fixedParams, bounds)
  for param in DrcParam:
    result.lower[param] = box.lower[ord(param)]
    result.upper[param] = box.upper[ord(param)]
  result.loss = opts.loss
  result.tuning = resolveTuning(opts)

  let nFree = drcParamCount - card(result.fixed)
  if result.n < max(nFree, 1):
    raise newException(
      ValueError,
      "need at least " & $max(nFree, 1) & " observations with " &
        "non-zero weight to fit " & $nFree & " parameters, got " & $result.n,
    )

  let
    x = logDoses(doses)
    y = result.responses
    w = result.weights
    free = freeMask(result.fixed)
  var p0 =
    if warm:
      start
    else:
      guessStart(x, y, w, fixedValues, result.fixed)
  for j in 0 ..< drcParamCount:
    p0[j] = clamp(p0[j], box.lower[j], box.upper[j])
  var best =
    if warm:
      levenbergMarquardt(ll5, x, y, w, p0, free, box.lower, box.upper, opts.lm)
    else:
      searchLM(x, y, w, p0, free, box, opts)
  result.converged = best.status.converged
  if opts.fallback and not result.converged:
    let fb =
      fallback(x, y, w, if best.ssr.isFinite: best.params else: p0, free, box, opts)
    best = fb.fit
    result.usedFallback = true
    result.fallbackConverged = fb.converged
    result.converged = fb.converged

  # A robust fit reweights from the least-squares fit.
  let sst = totalSumOfSquares(y, w)
  var
    robust: IrlsResult[LMResult[drcParamCount]]
    u: seq[float] # residuals in units of the robust scale
  if result.loss != Squared:
    robust = robustFit(
      ll5, x, y, w, best, free, box.lower, box.upper, result.loss, result.tuning, sst,
      result.n, opts,
    )
    best = robust.fit
    result.converged = best.status.converged
    if not robust.converged:
      result.warnings.incl fwRobustNotConverged
    (u, result.robustWeights) = robust.robustUnits(w, result.loss, result.tuning)

  result.params = best.params.toParams
  result.iterations = best.iterations
  result.status = best.status
  if not result.converged:
    result.warnings.incl fwNotConverged
  result.converged = result.warnings * {fwNotConverged, fwRobustNotConverged} == {}
  result.flagBounds
  result.flagEc50Range
  result.setUncertainty(x, robust.scale, u)
  result.setGoodnessOfFit(sst)

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
  fitCurveFrom(
    model, doses, responses, weights, fixed, bounds, opts, default(ParamVec), false
  )

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

func robustWeightNote(robustWeights: seq[float], n: int): string =
  ## How many observations a robust fit downweighted or rejected.
  var strong, rejected = 0
  for rw in robustWeights:
    if rw == 0.0:
      inc rejected
    elif rw < 0.5:
      inc strong
  # Observations the user excluded (weight 0) are not counted as rejected.
  rejected -= robustWeights.len - n
  &"Robust weight below 0.5: {strong} observations; 0 (rejected): {rejected}"

func fallbackNote(
    loss: Loss, usedFallback, fallbackConverged: bool, status: LMStatus, iterations: int
): string =
  ## How the fit stopped, for `summary`.
  if usedFallback:
    let outcome = if fallbackConverged: "converged" else: "did not converge"
    result =
      "Levenberg–Marquardt did not converge; the Nelder–Mead fallback " & outcome
    if loss == Squared:
      result.add &" (final polishing run: {status} after {iterations} iterations)"
    else:
      result.add &". Last reweighting round: {status} after {iterations} iterations"
  else:
    result = &"Stopped after {iterations} iterations: {status}"

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
    result.add &"Robust residual scale (MAD): {fit.sigma:.6g} on {fit.df} degrees of freedom\n"
    result.add robustWeightNote(fit.robustWeights, fit.n) & "\n"
    result.add &"R²: {fit.rSquared:.6g}\n"
  result.add fallbackNote(
    fit.loss, fit.usedFallback, fit.fallbackConverged, fit.status, fit.iterations
  )
  for warning in fit.warnings:
    result.add &"\nWarning: {warning}"

include drc/multi

include drc/bootstrap
