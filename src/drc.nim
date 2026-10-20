## Dose-response curve fitting.
##
## .. code-block:: nim
##   let fit = fitLL4(doses, responses, fixed = [(Bottom, 0.0)])
##   echo fit.ec50, " ", fit.ec50CI
##   echo fit.summary

import std/[math, strformat]
import drc/[linalg, lm, models, neldermead, robust, stats]

export lm, models, neldermead, robust, stats

type
  Model* = enum
    LL3 = "3PL"  ## 4PL with Bottom fixed
    LL4 = "4PL"
    LL5 = "5PL"

  FitWarning* = enum
    fwNotConverged = "the optimiser did not converge"
    fwSingular = "parameters are not identifiable (JᵀWJ is singular); standard errors unavailable"
    fwNoResidualDf = "no residual degrees of freedom; standard errors unavailable"
    fwEc50OutsideRange = "EC50 lies outside the tested dose range"
    fwRobustNotConverged = "the robust reweighting did not converge"

  FitOptions* = object
    lm*: LMOptions = LMOptions()
    nelderMead*: NMOptions = NMOptions()
    multiStart*: bool = true  ## also start from shifted LogEC50 and scaled Hill, keeping the best fit
    fallback*: bool = true    ## run Nelder–Mead when Levenberg–Marquardt does not converge
    loss*: Loss = Squared     ## Huber or Tukey for a robust fit that downweights outliers
    tuning*: float = NaN      ## the loss's tuning constant in units of the robust scale;
                              ## NaN uses the 95%-efficient default (1.345 Huber, 4.685 Tukey)
    irlsMaxIter*: int = 100   ## reweighting rounds for a robust fit
    irlsTol*: float = 1e-10   ## stop reweighting when the residuals change by less than this (relative)

  ConfInt* = tuple[lower, upper: float]

  Band* = tuple[fit, lower, upper: float]  ## a fitted response and its interval

  DrcFit* = object
    model*: Model
    params*: array[DrcParam, float]
    se*: array[DrcParam, float]   ## standard errors: 0 for fixed parameters, NaN if unavailable
    cov*: array[DrcParam, array[DrcParam, float]]  ## parameter covariance, σ²(JᵀWJ)⁻¹
    fixed*: set[DrcParam]         ## includes Asymmetry (= 1) for 3PL and 4PL fits
    n*: int              ## observations with non-zero weight
    df*: int             ## residual degrees of freedom: n − number of free parameters
    ssr*: float          ## weighted residual sum of squares
    sigma*: float        ## residual standard error, sqrt(ssr / df); for a robust fit, the
                         ## robust (MAD) scale of the residuals
    rSquared*: float     ## 1 − SSR / (weighted) total sum of squares
    logLik*: float       ## Gaussian log-likelihood, as R's logLik.nls (NaN for robust fits)
    aic*: float          ## AIC, counting σ as a parameter, as R and drc (NaN for robust fits)
    bic*: float          ## (NaN for robust fits)
    loss*: Loss
    tuning*: float              ## the loss's tuning constant (robust fits)
    robustWeights*: seq[float]  ## final ψ(u)/u per observation (robust fits; 0 where the
                                ## user's weight is 0); 1 means not downweighted
    logDoseRange*: tuple[lower, upper: float]  ## smallest and largest non-zero log10 dose
    converged*: bool
    usedFallback*: bool  ## Levenberg–Marquardt failed and Nelder–Mead was run
    iterations*: int     ## iterations of the final Levenberg–Marquardt run
    status*: LMStatus    ## how the final Levenberg–Marquardt run stopped; with the
                         ## fallback, `converged` may be true when this is not
    warnings*: set[FitWarning]

  ParamVec = array[drcParamCount, float]

func toVec(p: array[DrcParam, float]): ParamVec =
  for param in DrcParam: result[ord(param)] = p[param]

func toParams(v: ParamVec): array[DrcParam, float] =
  for param in DrcParam: result[param] = v[ord(param)]

proc predict*(fit: DrcFit; dose: float): float =
  var grad: ParamVec
  ll5(log10(dose), fit.params.toVec, grad)

proc tQuantile(fit: DrcFit; level: float): float =
  if not (level > 0.0 and level < 1.0):
    raise newException(ValueError, "confidence level must be in (0, 1), got " & $level)
  studentTQuantile(0.5 + 0.5 * level, float(fit.df))

func waldInterval(fit: DrcFit; param: DrcParam; tq: float): ConfInt =
  let half = tq * fit.se[param]
  (fit.params[param] - half, fit.params[param] + half)

proc confint*(fit: DrcFit; param: DrcParam; level = 0.95): ConfInt =
  ## Wald confidence interval using the t distribution with `fit.df` degrees of freedom.
  fit.waldInterval(param, fit.tQuantile(level))

proc checkPercent(pct: float) =
  if not (pct > 0.0 and pct < 100.0):
    raise newException(ValueError, "percentage must be in (0, 100), got " & $pct)

func logEcxAndGradient(fit: DrcFit; pct: float): (float, array[DrcParam, float]) =
  ## log10(ECx) and its gradient with respect to the parameters.
  ##
  ## The fraction of the change from the low-dose plateau is sᴬ for a rising
  ## curve and 1 − sᴬ for a falling one, where s = 1/(1 + 10^z) and
  ## z = (LogEC50 − x)·Hill. Solving for x gives x = LogEC50 − z/Hill with
  ## z = log10(r^(−1/A) − 1), r = pct/100 (rising) or 1 − pct/100 (falling).
  let
    h = fit.params[Hill]
    a = fit.params[Asymmetry]
    r = if h > 0.0: pct / 100.0 else: 1.0 - pct / 100.0
    v = pow(r, -1.0 / a)
    z = log10(v - 1.0)
  var grad: array[DrcParam, float]
  grad[LogEC50] = 1.0
  grad[Hill] = z / (h * h)
  grad[Asymmetry] = -v * ln(r) / (a * a * (v - 1.0) * ln(10.0) * h)
  (fit.params[LogEC50] - z / h, grad)

proc ecx*(fit: DrcFit; pct: float): float =
  ## Dose giving `pct`% of the change from the low-dose plateau to the
  ## high-dose plateau (drc's relative ED).
  checkPercent(pct)
  pow(10.0, fit.logEcxAndGradient(pct)[0])

proc ec50*(fit: DrcFit): float =
  ## The half-way dose. For 3PL and 4PL fits this is 10^LogEC50.
  fit.ecx(50.0)

func ecxInterval(fit: DrcFit; pct, tq: float): ConfInt =
  # Delta method on log10(ECx).
  let (center, grad) = fit.logEcxAndGradient(pct)
  var variance = 0.0
  for i in LogEC50 .. Asymmetry:  # ECx does not depend on Bottom or Top
    for j in LogEC50 .. Asymmetry:
      variance += grad[i] * fit.cov[i][j] * grad[j]
  let half = tq * sqrt(variance)
  (pow(10.0, center - half), pow(10.0, center + half))

proc ecxCI*(fit: DrcFit; pct: float; level = 0.95): ConfInt =
  ## Confidence interval for `ecx`, by the delta method on log10 dose and
  ## then back-transformed, so it is asymmetric on the dose scale.
  checkPercent(pct)
  fit.ecxInterval(pct, fit.tQuantile(level))

proc ec50CI*(fit: DrcFit; level = 0.95): ConfInt = fit.ecxCI(50.0, level)

proc bandAt(fit: DrcFit; dose, tq, extraVariance: float): Band =
  ## The fitted response at `dose` ± tq·sqrt(Var(f̂) + extraVariance), with
  ## Var(f̂) = gᵀ·cov·g by the delta method.
  if not (dose >= 0.0 and dose < Inf):
    raise newException(ValueError, "dose must be finite and non-negative, got " & $dose)
  var grad: ParamVec
  let value = ll5(log10(dose), fit.params.toVec, grad)
  var variance = extraVariance
  for i in DrcParam:
    for j in DrcParam:
      variance += grad[ord(i)] * fit.cov[i][j] * grad[ord(j)]
  let half = tq * sqrt(variance)
  (value, value - half, value + half)

proc confidenceBand*(fit: DrcFit; doses: openArray[float]; level = 0.95): seq[Band] =
  ## Pointwise confidence band for the fitted curve: where the true mean
  ## response lies at each dose.
  let tq = fit.tQuantile(level)
  for dose in doses: result.add fit.bandAt(dose, tq, 0.0)

proc confidenceBand*(fit: DrcFit; dose: float; level = 0.95): Band =
  fit.confidenceBand([dose], level)[0]

proc predictionBand*(fit: DrcFit; doses: openArray[float]; level = 0.95; weight = 1.0): seq[Band] =
  ## Pointwise prediction band: where a new observation at each dose is
  ## expected to fall. It adds the residual variance σ²/`weight` to the
  ## confidence band, with `weight` on the same relative scale as the fit's.
  if not (weight > 0.0 and weight < Inf):
    raise newException(ValueError, "weight must be finite and positive, got " & $weight)
  let tq = fit.tQuantile(level)
  for dose in doses: result.add fit.bandAt(dose, tq, fit.sigma ^ 2 / weight)

proc predictionBand*(fit: DrcFit; dose: float; level = 0.95; weight = 1.0): Band =
  fit.predictionBand([dose], level, weight)[0]

func better(a, b: LMResult[drcParamCount]): bool =
  ## A lower SSR wins, and any finite SSR beats NaN. Within rounding of each
  ## other, a converged fit beats an unconverged one; otherwise the earlier
  ## (default) start keeps its place.
  if b.ssr.isNaN: not a.ssr.isNaN
  elif a.ssr < b.ssr * (1.0 - 1e-10): true
  elif a.ssr <= b.ssr * (1.0 + 1e-10): a.status.converged and not b.status.converged
  else: false

proc multiStart(x, y, w: openArray[float]; p0: ParamVec; free: array[drcParamCount, bool];
                opts: FitOptions): LMResult[drcParamCount] =
  ## Levenberg–Marquardt from `p0` and, with `opts.multiStart`, from starts with
  ## LogEC50 shifted by ±1 and Hill scaled by 0.4 or 2.5; the best fit wins.
  result = levenbergMarquardt(ll5, x, y, w, p0, free, opts.lm)
  if not opts.multiStart: return
  for shift in [0.0, -1.0, 1.0]:
    for scale in [1.0, 0.4, 2.5]:
      if (shift == 0.0 and scale == 1.0) or (shift != 0.0 and not free[ord(LogEC50)]) or
         (scale != 1.0 and not free[ord(Hill)]):
        continue
      var p = p0
      p[ord(LogEC50)] += shift
      p[ord(Hill)] *= scale
      let r = levenbergMarquardt(ll5, x, y, w, p, free, opts.lm)
      if r.better(result): result = r

proc searchLM(x, y, w: openArray[float]; p0: ParamVec; free: array[drcParamCount, bool];
              opts: FitOptions): LMResult[drcParamCount] =
  if not free[ord(Asymmetry)]:
    return multiStart(x, y, w, p0, free, opts)
  # 5PL: fit the symmetric curve first, then release Asymmetry from there.
  var symmetric = free
  symmetric[ord(Asymmetry)] = false
  let base = multiStart(x, y, w, p0, symmetric, opts)
  result = base
  let released = levenbergMarquardt(ll5, x, y, w, base.params, free, opts.lm)
  if released.better(result): result = released
  if opts.multiStart:
    for a in [0.5, 2.0]:
      var p = base.params
      p[ord(Asymmetry)] = a
      let r = levenbergMarquardt(ll5, x, y, w, p, free, opts.lm)
      if r.better(result): result = r

proc standardisedResiduals(x, y, w: openArray[float]; p: ParamVec): seq[float] =
  ## √wᵢ·(yᵢ − f(xᵢ)), and 0 for observations with zero weight.
  result = newSeq[float](x.len)
  var grad: ParamVec
  for i in 0 ..< x.len:
    if w[i] > 0.0: result[i] = sqrt(w[i]) * (y[i] - ll5(x[i], p, grad))

type IrlsResult = object
  fit: LMResult[drcParamCount]
  scale: float
  residuals: seq[float]  ## standardised, at fit.params
  converged: bool

proc irls(x, y, w: openArray[float]; start: LMResult[drcParamCount];
          free: array[drcParamCount, bool]; loss: Loss; c: float; opts: FitOptions): IrlsResult =
  ## M-estimation by iteratively reweighted least squares, as R's
  ## robustbase::nlrob. Each round takes the robust scale s of the current
  ## residuals, weights observation i by wᵢ·ψ(uᵢ)/uᵢ with uᵢ = rᵢ/s, and refits.
  ## It stops when the residuals change by less than `opts.irlsTol` (relative).
  result.fit = start
  result.residuals = standardisedResiduals(x, y, w, start.params)
  var
    used: seq[float]
    robustW = newSeq[float](x.len)
  for _ in 1 .. opts.irlsMaxIter:
    used.setLen 0
    for i in 0 ..< x.len:
      if w[i] > 0.0: used.add result.residuals[i]
    result.scale = madScale(used)
    if not (result.scale > 0.0):
      return  # most residuals are exactly zero, so there is no scale to weight by
    for i in 0 ..< x.len:
      robustW[i] = w[i] * psiWeight(loss, result.residuals[i] / result.scale, c)
    # Copy the start: the call writes its result straight into result.fit, and
    # Nim passes the start array by reference, so it would read zeros.
    let current = result.fit.params
    result.fit = levenbergMarquardt(ll5, x, y, robustW, current, free, opts.lm)
    let residuals = standardisedResiduals(x, y, w, result.fit.params)
    var change, size = 0.0
    for i in 0 ..< x.len:
      change += (residuals[i] - result.residuals[i]) ^ 2
      size += result.residuals[i] ^ 2
    result.residuals = residuals
    if sqrt(change / max(size, 1e-20)) <= opts.irlsTol:
      result.converged = true
      return

proc fitCurve(model: Model; doses, responses, weights: openArray[float];
              fixed: openArray[(DrcParam, float)]; opts: FitOptions): DrcFit =
  result.model = model
  if doses.len != responses.len:
    raise newException(ValueError, "doses and responses differ in length")
  if weights.len != 0 and weights.len != doses.len:
    raise newException(ValueError, "weights and doses differ in length")

  var
    x, w = newSeq[float](doses.len)
    guessX, guessY: seq[float]  # observations with non-zero weight
    sumW, sumWY, sumLogW = 0.0
    lo = Inf
    hi = NegInf
  for i in 0 ..< doses.len:
    if not (doses[i] >= 0.0 and doses[i] < Inf):
      raise newException(ValueError, "dose must be finite and non-negative, got " & $doses[i])
    if not responses[i].isFinite:
      raise newException(ValueError, "response must be finite, got " & $responses[i])
    x[i] = log10(doses[i])
    w[i] = if weights.len == 0: 1.0 else: weights[i]
    if not (w[i] >= 0.0 and w[i] < Inf):
      raise newException(ValueError, "weight must be finite and non-negative, got " & $w[i])
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
  result.logDoseRange = if lo <= hi: (lo, hi) else: (NaN, NaN)

  var
    fixedValues: ParamVec
    free: array[drcParamCount, bool] = [true, true, true, true, true]
  for (param, value) in fixed:
    if param == Asymmetry and model != LL5:
      raise newException(ValueError, "Asymmetry can only be fixed in a 5PL fit")
    if not value.isFinite:
      raise newException(ValueError, "fixed " & $param & " must be finite, got " & $value)
    fixedValues[ord(param)] = value
    free[ord(param)] = false
    result.fixed.incl param
  if model != LL5:
    fixedValues[ord(Asymmetry)] = 1.0
    free[ord(Asymmetry)] = false
    result.fixed.incl Asymmetry
  if not (fixedValues[ord(Asymmetry)] > 0.0) and Asymmetry in result.fixed:
    raise newException(ValueError, "Asymmetry must be positive")

  result.loss = opts.loss
  if opts.loss != Squared:
    result.tuning = if opts.tuning.isNaN: defaultTuning(opts.loss) else: opts.tuning
    if not (result.tuning > 0.0 and result.tuning < Inf):
      raise newException(ValueError, "tuning must be finite and positive, got " & $opts.tuning)

  let nFree = drcParamCount - card(result.fixed)
  if result.n < max(nFree, 1):
    raise newException(ValueError, "need at least " & $max(nFree, 1) & " observations with " &
                       "non-zero weight to fit " & $nFree & " parameters, got " & $result.n)

  let p0 = drcGuess(guessX, guessY, fixedValues, result.fixed)
  var best = searchLM(x, responses, w, p0, free, opts)
  result.converged = best.status.converged
  if opts.fallback and not result.converged:
    # Nelder–Mead needs no derivatives; Levenberg–Marquardt then polishes its result.
    let y = @responses
    let start = if best.ssr.isFinite: best.params else: p0
    let nm = nelderMead(proc (p: ParamVec): float = weightedSsr(ll5, x, y, w, p),
                        start, free, opts.nelderMead)
    # Each stage only goes downhill from the best fit so far, so the polished
    # fit is never worse.
    best = levenbergMarquardt(ll5, x, y, w, nm.params, free, opts.lm)
    result.usedFallback = true
    result.converged = best.status.converged or nm.converged

  # A robust fit reweights from the least-squares fit.
  var robustFit: IrlsResult
  if result.loss != Squared:
    var start = best
    if result.loss == Tukey:
      # Tukey's redescending ψ can stop in a local minimum, so start from Huber.
      start = irls(x, responses, w, best, free, Huber, defaultTuning(Huber), opts).fit
    robustFit = irls(x, responses, w, start, free, result.loss, result.tuning, opts)
    best = robustFit.fit
    result.converged = best.status.converged
    if not robustFit.converged: result.warnings.incl fwRobustNotConverged
    result.robustWeights = newSeq[float](x.len)
    for i in 0 ..< x.len:
      if w[i] > 0.0:
        result.robustWeights[i] =
          psiWeight(result.loss, robustFit.residuals[i] / robustFit.scale, result.tuning)

  result.params = best.params.toParams
  result.ssr = weightedSsr(ll5, x, responses, w, best.params)
  result.iterations = best.iterations
  result.status = best.status
  result.df = result.n - nFree

  if not result.converged:
    result.warnings.incl fwNotConverged
  let logEC50 = log10(result.ec50)
  if not (logEC50 >= result.logDoseRange.lower and logEC50 <= result.logDoseRange.upper):
    result.warnings.incl fwEc50OutsideRange

  # Uncertainty: cov = varianceScale · (JᵀWJ)⁻¹ over the free parameters, with
  # varianceScale = σ² for least squares.
  var varianceScale = NaN
  if result.loss == Squared:
    result.sigma = if result.df > 0: sqrt(result.ssr / float(result.df)) else: NaN
    varianceScale = result.sigma ^ 2
  else:
    # Huber's correction, as MASS::rlm: κ²·S / mean(ψ')², where
    # S = s²·Σψ(uᵢ)²/(n − p) and κ = 1 + p·var(ψ')/(n·mean(ψ')²).
    let s = robustFit.scale
    result.sigma = s
    var sumPsi2, sumD, sumD2 = 0.0
    for i in 0 ..< x.len:
      if w[i] > 0.0:
        let u = robustFit.residuals[i] / s
        sumPsi2 += psi(result.loss, u, result.tuning) ^ 2
        let d = psiDerivative(result.loss, u, result.tuning)
        sumD += d
        sumD2 += d * d
    let
      n = float(result.n)
      meanD = sumD / n
      varD = (sumD2 - n * meanD * meanD) / (n - 1.0)
      kappa = 1.0 + float(nFree) * varD / (n * meanD * meanD)
    if result.df > 0 and meanD > 0.0:
      varianceScale = (kappa / meanD) ^ 2 * s * s * sumPsi2 / float(result.df)
  var
    jtwj, unscaled: array[drcParamCount, array[drcParamCount, float]]
    gradient: ParamVec
  discard normalEquations(ll5, x, responses, w, best.params, jtwj, gradient)
  if not unscaledCovariance(jtwj, free, unscaled): result.warnings.incl fwSingular
  if result.df == 0: result.warnings.incl fwNoResidualDf
  for i in DrcParam:
    for j in DrcParam:
      # Fixed entries are exactly 0; a NaN σ or (JᵀWJ)⁻¹ carries through to the rest.
      result.cov[i][j] =
        if i in result.fixed or j in result.fixed: 0.0
        else: varianceScale * unscaled[ord(i)][ord(j)]
  for param in DrcParam:
    result.se[param] = sqrt(result.cov[param][param])

  # Goodness of fit.
  let meanY = sumWY / sumW
  var sst = 0.0
  for i in 0 ..< responses.len:
    sst += w[i] * (responses[i] - meanY) ^ 2
  result.rSquared = 1.0 - result.ssr / sst
  if result.loss == Squared:
    result.logLik = gaussianLogLik(result.ssr, result.n, sumLogW)
    let k = float(nFree + 1)  # free parameters plus σ
    result.aic = -2.0 * result.logLik + 2.0 * k
    result.bic = -2.0 * result.logLik + ln(float(result.n)) * k
  else:
    # A robust fit does not maximise the Gaussian likelihood.
    (result.logLik, result.aic, result.bic) = (NaN, NaN, NaN)

proc fitLL4*(doses, responses: openArray[float];
             weights: openArray[float] = [];
             fixed: openArray[(DrcParam, float)] = [];
             opts = FitOptions()): DrcFit =
  ## Fits a four-parameter logistic to `responses` at `doses` (linear scale;
  ## zero doses are treated as untreated controls). Parameters in `fixed`
  ## are held at the given value. Weights are relative: the error variance
  ## of observation i is taken to be σ²/wᵢ, with σ² estimated from the fit.
  fitCurve(LL4, doses, responses, weights, fixed, opts)

proc fitLL3*(doses, responses: openArray[float];
             weights: openArray[float] = [];
             bottom = 0.0;
             opts = FitOptions()): DrcFit =
  ## Fits a three-parameter logistic: the 4PL with Bottom fixed at `bottom`.
  fitCurve(LL3, doses, responses, weights, [(Bottom, bottom)], opts)

proc fitLL5*(doses, responses: openArray[float];
             weights: openArray[float] = [];
             fixed: openArray[(DrcParam, float)] = [];
             opts = FitOptions()): DrcFit =
  ## Fits a five-parameter logistic, for curves that are not symmetric about
  ## their midpoint. See `ll5` for the parameterisation; `ec50` and `ecx`
  ## give the actual doses for a given response.
  fitCurve(LL5, doses, responses, weights, fixed, opts)

proc summary*(fit: DrcFit; level = 0.95): string =
  ## A printable table of estimates, standard errors and confidence intervals.
  let
    pct = &"{100.0 * level:g}%"
    tq = fit.tQuantile(level)
  result.add &"{fit.model} fit to {fit.n} observations"
  if fit.loss != Squared:
    result.add &", robust ({fit.loss}, tuning {fit.tuning:g})"
  result.add "\n\n"
  result.add &"""{"":<10}{"Estimate":>12}{"Std.Error":>12}{"Lower " & pct:>14}{"Upper " & pct:>14}""" & "\n"
  for param in DrcParam:
    if param == Asymmetry and fit.model != LL5:
      continue
    if param in fit.fixed:
      result.add &"{$param:<10}{fit.params[param]:>12.6g}{\"(fixed)\":>12}\n"
    else:
      let ci = fit.waldInterval(param, tq)
      result.add &"{$param:<10}{fit.params[param]:>12.6g}{fit.se[param]:>12.6g}{ci.lower:>14.6g}{ci.upper:>14.6g}\n"
  let ci = fit.ecxInterval(50.0, tq)
  result.add &"{\"EC50\":<10}{fit.ec50:>12.6g}{\"\":>12}{ci.lower:>14.6g}{ci.upper:>14.6g}\n\n"
  if fit.loss == Squared:
    result.add &"Residual standard error: {fit.sigma:.6g} on {fit.df} degrees of freedom\n"
    result.add &"R²: {fit.rSquared:.6g}   AIC: {fit.aic:.6g}   BIC: {fit.bic:.6g}\n"
  else:
    var strong, rejected = 0
    for rw in fit.robustWeights:
      if rw == 0.0: inc rejected
      elif rw < 0.5: inc strong
    # Observations the user excluded (weight 0) are not counted as rejected.
    rejected -= fit.robustWeights.len - fit.n
    result.add &"Robust residual scale (MAD): {fit.sigma:.6g} on {fit.df} degrees of freedom\n"
    result.add &"Robust weight below 0.5: {strong} observations; 0 (rejected): {rejected}\n"
    result.add &"R²: {fit.rSquared:.6g}\n"
  if fit.usedFallback:
    let outcome = if fit.converged: "converged" else: "did not converge"
    result.add "Levenberg–Marquardt did not converge; the Nelder–Mead fallback " & outcome &
               &" (final polishing run: {fit.status} after {fit.iterations} iterations)"
  else:
    result.add &"Stopped after {fit.iterations} iterations: {fit.status}"
  for warning in fit.warnings:
    result.add &"\nWarning: {warning}"
