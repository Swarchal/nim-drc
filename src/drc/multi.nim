## Joint fits of several dose-response curves, with parameters either shared
## by all curves or separate for each, as R drc's `curveid` with `pmodels`.
## Included by drc.nim, so it uses that module's private helpers.
##
## The joint parameter vector has one slot per shared parameter and one per
## curve for each unshared parameter, ordered by `DrcParam` and then by
## curve. Fixed parameters are shared by all curves and pinned.

import std/[algorithm, hashes, strutils]

type
  MultiFit* = object
    model*: Model
    shared*: set[DrcParam] ## one value for all curves; includes the fixed parameters
    fixed*: set[DrcParam] ## includes Asymmetry (= 1) for 3PL and 4PL fits
    lower*, upper*: array[DrcParam, float]
      ## the bounds on each parameter, for every curve: −Inf and Inf when unbounded
    atBound*: set[DrcParam] ## free parameters that ended on a bound in any curve
    curves*: seq[DrcFit] ## each curve as a single fit (see `fitLL4Multi`)
    index*: seq[array[DrcParam, int]]
      ## index[c][param]: slot of curve c's parameter in `params`
    params*: seq[float] ## the joint parameter vector
    free*: seq[bool] ## false for the slots of fixed parameters
    cov*: seq[seq[float]] ## joint covariance σ²(JᵀWJ)⁻¹; 0 for fixed slots
    n*: int ## observations with non-zero weight
    df*: int ## n − number of free joint parameters
    ssr*: float
    sigma*: float ## pooled residual standard error, sqrt(ssr / df)
    rSquared*: float
    logLik*: float ## Gaussian log-likelihood, as R's logLik.nls
    aic*: float
    bic*: float
    converged*: bool
    iterations*: int ## iterations of the Levenberg–Marquardt run that won
    status*: LMStatus
    warnings*: set[FitWarning] ## includes fwEc50OutsideRange if any curve has it
    fixedValues: ParamVec
    curveIds: seq[int]
    dataHash: Hash

  FTest* = tuple[f: float, df1, df2: int, p: float]

proc median(values: openArray[float]): float =
  var s = @values
  s.sort()
  let m = s.len div 2
  if s.len mod 2 == 1:
    s[m]
  else:
    0.5 * (s[m - 1] + s[m])

proc sharedStart(param: DrcParam, values: openArray[float]): float =
  ## One starting value for a shared parameter from the curves' own values.
  case param
  of Hill:
    # Majority direction, typical steepness.
    var rising = 0
    var steepness: seq[float]
    for v in values:
      if v > 0.0:
        inc rising
      steepness.add abs(v)
    (if 2 * rising >= values.len: 1.0 else: -1.0) * median(steepness)
  of Asymmetry:
    var logs: seq[float]
    for v in values:
      logs.add ln(v)
    exp(median(logs))
  else:
    median(values)

proc fitMulti(
    model: Model,
    doses, responses: openArray[float],
    curve: openArray[int],
    shared: set[DrcParam],
    weights: openArray[float],
    fixed: openArray[(DrcParam, float)],
    bounds: openArray[(DrcParam, Slice[float])],
    opts: FitOptions,
): MultiFit =
  if opts.loss != Squared:
    raise newException(ValueError, "robust losses are not supported for joint fits")
  if doses.len != responses.len or curve.len != doses.len:
    raise newException(ValueError, "doses, responses and curve differ in length")
  if weights.len != 0 and weights.len != doses.len:
    raise newException(ValueError, "weights and doses differ in length")
  if doses.len == 0:
    raise newException(ValueError, "no observations")

  let (fixedValues, fixedParams) = resolveFixed(model, fixed)
  result.model = model
  result.fixed = fixedParams
  result.shared = shared + fixedParams
  let box = resolveBounds(model, fixedParams, bounds)
  for param in DrcParam:
    result.lower[param] = box.lower[ord(param)]
    result.upper[param] = box.upper[ord(param)]
  result.fixedValues = fixedValues
  result.curveIds = @curve

  var m = 0
  for c in curve:
    if c < 0:
      raise newException(ValueError, "curve ids must be non-negative, got " & $c)
    m = max(m, c + 1)
  var
    x, w = newSeq[float](doses.len)
    nCurve = newSeq[int](m)
    sumW, sumWY, sumLogW = 0.0
  for i in 0 ..< doses.len:
    x[i] = log10(doses[i])
    w[i] =
      if weights.len == 0:
        1.0
      else:
        weights[i]
    checkObservation(doses[i], responses[i], w[i])
    if w[i] > 0.0:
      inc nCurve[curve[i]]
      inc result.n
      sumW += w[i]
      sumWY += w[i] * responses[i]
      sumLogW += ln(w[i])
  for c in 0 ..< m:
    if nCurve[c] == 0:
      raise newException(
        ValueError, "curve " & $c & " has no observations with non-zero weight"
      )
  result.dataHash = !$(hash(@doses) !& hash(@responses) !& hash(w))

  # Layout of the joint vector; fixed slots hold their value from the start.
  result.index = newSeq[array[DrcParam, int]](m)
  var lower, upper: seq[float] # bounds on each slot
  for param in DrcParam:
    let copies = if param in result.shared: 1 else: m
    for k in 0 ..< copies:
      for c in 0 ..< m:
        if copies == 1 or c == k:
          result.index[c][param] = result.params.len
      result.params.add fixedValues[ord(param)]
      result.free.add param notin fixedParams
      lower.add box.lower[ord(param)]
      upper.add box.upper[ord(param)]
  var nFree = 0
  for f in result.free:
    if f:
      inc nFree
  let unshared = drcParamCount - card(result.shared)
  if result.n < max(nFree, 1):
    raise newException(
      ValueError,
      "need at least " & $max(nFree, 1) & " observations with " &
        "non-zero weight to fit " & $nFree & " parameters, got " & $result.n,
    )
  for c in 0 ..< m:
    if nCurve[c] < unshared:
      raise newException(
        ValueError,
        "curve " & $c & " has " & $nCurve[c] &
          " observations with non-zero weight, fewer than its " & $unshared &
          " unshared parameters",
      )

  # Starting values: each curve fitted on its own where it has enough data,
  # with shared parameters started at the typical value over curves.
  var singleOpts = opts
  singleOpts.fallback = false
  let curveFree = drcParamCount - card(fixedParams)
  var single = newSeq[ParamVec](m)
  for c in 0 ..< m:
    var dc, rc, wc, gx, gy: seq[float]
    for i in 0 ..< doses.len:
      if curve[i] == c:
        dc.add doses[i]
        rc.add responses[i]
        wc.add w[i]
        if w[i] > 0.0:
          gx.add x[i]
          gy.add responses[i]
    var ok = false
    if nCurve[c] >= max(curveFree, 1):
      single[c] = fitCurve(model, dc, rc, wc, fixed, bounds, singleOpts).params.toVec
      ok = true
      for v in single[c]:
        ok = ok and v.isFinite
    if not ok:
      single[c] = drcGuess(gx, gy, fixedValues, fixedParams)
      for j in 0 ..< drcParamCount:
        single[c][j] = clamp(single[c][j], box.lower[j], box.upper[j])

  var startA = result.params
  for param in DrcParam:
    if param in fixedParams:
      continue
    if param in shared:
      var values: seq[float]
      for c in 0 ..< m:
        values.add single[c][ord(param)]
      startA[result.index[0][param]] =
        clamp(sharedStart(param, values), box.lower[ord(param)], box.upper[ord(param)])
    else:
      for c in 0 ..< m:
        startA[result.index[c][param]] = single[c][ord(param)]
  var starts = @[startA]
  if shared - fixedParams != {} and result.n >= max(curveFree, 1):
    # The pooled curve is a good guess for the shared parameters.
    let pooled = fitCurve(model, doses, responses, w, fixed, bounds, singleOpts).params
    var startB = startA
    for param in shared - fixedParams:
      if pooled[param].isFinite:
        startB[result.index[0][param]] = pooled[param]
    starts.add startB
  if opts.multiStart:
    # As the single-curve search: shift every LogEC50 by ±1 and scale every Hill.
    for shift in [0.0, -1.0, 1.0]:
      for scale in [1.0, 0.4, 2.5]:
        if (shift == 0.0 and scale == 1.0) or (shift != 0.0 and LogEC50 in fixedParams) or
            (scale != 1.0 and Hill in fixedParams):
          continue
        var p = startA
        for param in [LogEC50, Hill]:
          var done: seq[int]
          for c in 0 ..< m:
            let slot = result.index[c][param]
            if slot notin done:
              done.add slot
              if param == LogEC50:
                p[slot] += shift
              else:
                p[slot] *= scale
        starts.add p

  let
    obs = block:
      var s = newSeq[float](doses.len)
      for i in 0 ..< s.len:
        s[i] = float(i)
      s
    index = result.index
    curveIds = result.curveIds
    ys = @responses
  # Observation i's curve gets its parameters from the joint vector; LM sees
  # the observation index as x.
  let joint = proc(i: float, p: seq[float], grad: var seq[float]): float =
    let k = int(i)
    let slots = index[curveIds[k]]
    var q, g5: ParamVec
    for param in DrcParam:
      q[ord(param)] = p[slots[param]]
    result = ll5(x[k], q, g5)
    grad.setZero()
    for param in DrcParam:
      grad[slots[param]] += g5[ord(param)]

  proc search(starts: seq[seq[float]], free: seq[bool]): LMDynResult =
    for k, start in starts:
      let r = levenbergMarquardt(joint, obs, ys, w, start, free, lower, upper, opts.lm)
      if k == 0 or r.better(result):
        result = r

  var best: LMDynResult
  if Asymmetry notin fixedParams:
    # As for one curve: fit symmetric curves first, then release Asymmetry.
    var symFree = result.free
    var symStarts = starts
    for c in 0 ..< m:
      symFree[result.index[c][Asymmetry]] = false
    for s in symStarts.mitems:
      for c in 0 ..< m:
        s[result.index[c][Asymmetry]] =
          clamp(1.0, box.lower[ord(Asymmetry)], box.upper[ord(Asymmetry)])
    best = search(symStarts, symFree)
    let released = search(@[best.params, startA], result.free)
    if released.better(best):
      best = released
  else:
    best = search(starts, result.free)

  result.params = best.params
  result.iterations = best.iterations
  result.status = best.status
  result.converged = best.status.converged
  if not result.converged:
    result.warnings.incl fwNotConverged
  for c in 0 ..< m:
    for param in DrcParam:
      let slot = result.index[c][param]
      if result.free[slot] and result.params[slot] in [lower[slot], upper[slot]]:
        result.atBound.incl param
  if result.atBound != {}:
    result.warnings.incl fwAtBound

  # Uncertainty, as for a single curve, over the joint vector, with parameters
  # at a bound counted as free (as R's nls with algorithm = "port").
  var
    jtwj, unscaled: seq[seq[float]]
    gradient: seq[float]
  result.ssr = normalEquations(joint, obs, responses, w, result.params, jtwj, gradient)
  if not unscaledCovariance(jtwj, result.free, unscaled):
    result.warnings.incl fwSingular
  result.df = result.n - nFree
  if result.df == 0:
    result.warnings.incl fwNoResidualDf
  result.sigma =
    if result.df > 0:
      sqrt(result.ssr / float(result.df))
    else:
      NaN
  let np = result.params.len
  result.cov = newSeq[seq[float]](np)
  for i in 0 ..< np:
    result.cov[i] = newSeq[float](np)
    for j in 0 ..< np:
      if result.free[i] and result.free[j]:
        result.cov[i][j] = result.sigma ^ 2 * unscaled[i][j]

  let meanY = sumWY / sumW
  var sst = 0.0
  for i in 0 ..< doses.len:
    sst += w[i] * (responses[i] - meanY) ^ 2
  result.rSquared = 1.0 - result.ssr / sst
  result.logLik = gaussianLogLik(result.ssr, result.n, sumLogW)
  let k = float(nFree + 1) # free parameters plus σ
  result.aic = -2.0 * result.logLik + 2.0 * k
  result.bic = -2.0 * result.logLik + ln(float(result.n)) * k

  # Each curve as a DrcFit: its parameters and covariance block, with the
  # joint σ and df so intervals and bands account for the pooled fit.
  for c in 0 ..< m:
    var v: DrcFit
    v.model = model
    v.fixed = fixedParams
    (v.lower, v.upper) = (result.lower, result.upper)
    v.loss = Squared
    for param in DrcParam:
      v.params[param] = result.params[index[c][param]]
      for other in DrcParam:
        v.cov[param][other] = result.cov[index[c][param]][index[c][other]]
    for param in DrcParam:
      v.se[param] = sqrt(v.cov[param][param])
    v.n = nCurve[c]
    v.df = result.df
    v.sigma = result.sigma
    var
      cw, cwy = 0.0
      lo = Inf
      hi = NegInf
      grad: ParamVec
    for i in 0 ..< doses.len:
      if curve[i] == c and w[i] > 0.0:
        cw += w[i]
        cwy += w[i] * responses[i]
        v.ssr += w[i] * (responses[i] - ll5(x[i], v.params.toVec, grad)) ^ 2
        if x[i] > NegInf:
          lo = min(lo, x[i])
          hi = max(hi, x[i])
    var cst = 0.0
    for i in 0 ..< doses.len:
      if curve[i] == c and w[i] > 0.0:
        cst += w[i] * (responses[i] - cwy / cw) ^ 2
    v.rSquared = 1.0 - v.ssr / cst
    (v.logLik, v.aic, v.bic) = (NaN, NaN, NaN)
    v.logDoseRange =
      if lo <= hi:
        (lo, hi)
      else:
        (NaN, NaN)
    v.converged = result.converged
    v.iterations = result.iterations
    v.status = result.status
    v.warnings = result.warnings - {fwAtBound}
    for param in DrcParam:
      if param notin fixedParams and v.params[param] in [v.lower[param], v.upper[param]]:
        v.atBound.incl param
    if v.atBound != {}:
      v.warnings.incl fwAtBound
    let logEC50 = log10(v.ec50)
    if not (logEC50 >= v.logDoseRange.lower and logEC50 <= v.logDoseRange.upper):
      v.warnings.incl fwEc50OutsideRange
      result.warnings.incl fwEc50OutsideRange
    result.curves.add v

proc fitLL4Multi*(
    doses, responses: openArray[float],
    curve: openArray[int],
    shared: set[DrcParam] = {},
    weights: openArray[float] = [],
    fixed: openArray[(DrcParam, float)] = [],
    bounds: openArray[(DrcParam, Slice[float])] = [],
    opts = FitOptions(),
): MultiFit =
  ## Fits four-parameter logistics to several curves at once. Observation i
  ## belongs to curve `curve[i]` (ids 0 ..< m, each used). Parameters in
  ## `shared` take one value for all curves; the rest are fitted per curve.
  ## For example, `shared = {Bottom, Top, Hill}` fits parallel curves that
  ## differ only in EC50. Parameters in `fixed` are held at one value for all
  ## curves, and those in `bounds` (as for `fitLL4`) are kept within one
  ## range for all curves.
  ##
  ## `fit.curves[c]` is curve c as a `DrcFit`, so `ec50`, `ecxCI`, `predict`
  ## and the bands work on it. Its covariance is curve c's block of the joint
  ## covariance, and its `df` and `sigma` are the joint fit's; its `ssr`,
  ## `rSquared` and `n` are its own. Only least squares is supported.
  fitMulti(LL4, doses, responses, curve, shared, weights, fixed, bounds, opts)

proc fitLL3Multi*(
    doses, responses: openArray[float],
    curve: openArray[int],
    shared: set[DrcParam] = {},
    weights: openArray[float] = [],
    bottom = 0.0,
    bounds: openArray[(DrcParam, Slice[float])] = [],
    opts = FitOptions(),
): MultiFit =
  ## As `fitLL4Multi`, with Bottom fixed at `bottom` for every curve.
  fitMulti(
    LL3, doses, responses, curve, shared, weights, [(Bottom, bottom)], bounds, opts
  )

proc fitLL5Multi*(
    doses, responses: openArray[float],
    curve: openArray[int],
    shared: set[DrcParam] = {},
    weights: openArray[float] = [],
    fixed: openArray[(DrcParam, float)] = [],
    bounds: openArray[(DrcParam, Slice[float])] = [],
    opts = FitOptions(),
): MultiFit =
  ## As `fitLL4Multi` for five-parameter logistics.
  fitMulti(LL5, doses, responses, curve, shared, weights, fixed, bounds, opts)

proc relativePotency*(
    fit: MultiFit, a, b: int, pct = 50.0, level = 0.95
): tuple[ratio: float, ci: ConfInt] =
  ## ECx of curve `a` divided by ECx of curve `b`, with a confidence interval
  ## by the delta method on log10 of the ratio (so it is asymmetric), using
  ## the joint covariance including the parts shared by the two curves.
  ## R drc's `EDcomp` applies the delta method to the ratio itself by default,
  ## which gives a symmetric interval.
  checkPercent(pct)
  if a notin 0 ..< fit.curves.len or b notin 0 ..< fit.curves.len:
    raise newException(ValueError, "no such curve: " & $a & " or " & $b)
  let
    (logA, gradA) = fit.curves[a].logEcxAndGradient(pct)
    (logB, gradB) = fit.curves[b].logEcxAndGradient(pct)
  var g = newSeq[float](fit.params.len)
  for param in DrcParam:
    g[fit.index[a][param]] += gradA[param]
    g[fit.index[b][param]] -= gradB[param]
  var variance = 0.0
  for i in 0 ..< g.len:
    if g[i] != 0.0:
      for j in 0 ..< g.len:
        if g[j] != 0.0:
          variance += g[i] * fit.cov[i][j] * g[j]
  let
    center = logA - logB
    half = fit.curves[a].tQuantile(level) * sqrt(variance)
  (pow(10.0, center), (pow(10.0, center - half), pow(10.0, center + half)))

proc fTest*(reduced, full: MultiFit): FTest =
  ## F test of a reduced joint fit (e.g. parallel curves) against a fuller one
  ## (e.g. separate curves) on the same data, as R's `anova` for nested
  ## models. A small p means the extra parameters improve the fit.
  if reduced.dataHash != full.dataHash:
    raise newException(ValueError, "the fits are not to the same data")
  if reduced.curveIds == full.curveIds:
    if not (reduced.shared >= full.shared and reduced.fixed >= full.fixed):
      raise newException(
        ValueError,
        "the reduced fit must share and fix at least " &
          "the parameters the full fit does",
      )
    for param in full.fixed:
      if reduced.fixedValues[ord(param)] != full.fixedValues[ord(param)]:
        raise
          newException(ValueError, "the fits fix " & $param & " at different values")
  if not (reduced.df > full.df):
    raise newException(
      ValueError, "the reduced fit must have more residual degrees of freedom"
    )
  let
    df1 = reduced.df - full.df
    f = ((reduced.ssr - full.ssr) / float(df1)) / (full.ssr / float(full.df))
  (f, df1, full.df, fTail(f, float(df1), float(full.df)))

proc summary*(fit: MultiFit, level = 0.95): string =
  ## A printable table of the joint estimates and each curve's EC50.
  let pct = &"{100.0 * level:g}%"
  let tq = fit.curves[0].tQuantile(level) # every curve has the joint df
  result.add &"{fit.model} joint fit to {fit.n} observations in {fit.curves.len} curves"
  var sharedNames: seq[string]
  for param in fit.shared - fit.fixed:
    sharedNames.add $param
  if sharedNames.len > 0:
    result.add "; shared: " & sharedNames.join(", ")
  result.add "\n\n"
  result.add &"""{"":<12}{"Estimate":>12}{"Std.Error":>12}{"Lower " & pct:>14}{"Upper " & pct:>14}""" &
    "\n"
  for param in DrcParam:
    if param == Asymmetry and fit.model != LL5:
      continue
    for c in 0 ..< fit.curves.len:
      if c > 0 and param in fit.shared:
        break
      let
        slot = fit.index[c][param]
        name =
          if param in fit.shared:
            $param
          else:
            &"{param}[{c}]"
        value = fit.params[slot]
      if param in fit.fixed:
        result.add &"{name:<12}{value:>12.6g}{\"(fixed)\":>12}\n"
      else:
        let se = sqrt(fit.cov[slot][slot])
        result.add &"{name:<12}{value:>12.6g}{se:>12.6g}{value - tq * se:>14.6g}{value + tq * se:>14.6g}"
        result.add fit.curves[c].boundNote(param) & "\n"
  for c, curve in fit.curves:
    let ci = curve.ecxInterval(50.0, tq)
    result.add &"""{"EC50[" & $c & "]":<12}{curve.ec50:>12.6g}{"":>12}{ci.lower:>14.6g}{ci.upper:>14.6g}""" &
      "\n"
  result.add &"\nResidual standard error: {fit.sigma:.6g} on {fit.df} degrees of freedom\n"
  result.add &"R²: {fit.rSquared:.6g}   AIC: {fit.aic:.6g}   BIC: {fit.bic:.6g}\n"
  result.add &"Stopped after {fit.iterations} iterations: {fit.status}"
  for warning in fit.warnings:
    result.add &"\nWarning: {warning}"
