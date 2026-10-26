## Joint fits of several dose-response curves, with parameters either shared
## by all curves or separate for each, as R drc's `curveid` with `pmodels`.
## Included by drc.nim, so it uses that module's private helpers.
##
## The joint parameter vector (private; see `param` and `cov`) has one slot
## per shared parameter and one per curve for each unshared parameter,
## ordered by `DrcParam` and then by curve. Fixed parameters are shared by all
## curves and pinned.

import std/[algorithm, hashes, strutils]

type MultiFit* = object
  model*: Model
  shared*: set[DrcParam] ## one value for all curves; includes the fixed parameters
  fixed*: set[DrcParam] ## includes Asymmetry (= 1) for 3PL and 4PL fits
  lower*, upper*: array[DrcParam, float]
    ## the bounds on each parameter, for every curve: −Inf and Inf when unbounded
  atBound*: set[DrcParam] ## free parameters that ended on a bound in any curve
  curves*: seq[DrcFit] ## each curve as a single fit (see `fitMulti`)
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
  slots: seq[array[DrcParam, int]]
    ## slots[c][param]: where curve c's parameter is in `joint`
  joint: seq[float] ## the joint parameter vector
  jointFree: seq[bool] ## false for the slots of fixed parameters
  jointCov: seq[seq[float]] ## joint covariance σ²(JᵀWJ)⁻¹; 0 for fixed slots
  fixedValues: ParamVec
  curveIds: seq[int]
  dataHash: Hash
  data: tuple[doses, responses, weights: seq[float]] ## as given, in order

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

func observationIndex(n: int): seq[float] =
  ## 0, 1, …, n − 1: the x that LM sees for a joint fit.
  result = newSeq[float](n)
  for i in 0 ..< n:
    result[i] = float(i)

proc jointModel(
    x: seq[float], slots: seq[array[DrcParam, int]], curveIds: seq[int]
): proc(i: float, p: seq[float], grad: var seq[float]): float =
  ## The joint model: observation i's curve gets its parameters from the
  ## joint vector, and LM sees the observation index as x.
  result = proc(i: float, p: seq[float], grad: var seq[float]): float =
    let k = int(i)
    let curveSlots = slots[curveIds[k]]
    var q, g5: ParamVec
    for param in DrcParam:
      q[ord(param)] = p[curveSlots[param]]
    result = ll5(x[k], q, g5)
    grad.setZero()
    for param in DrcParam:
      grad[curveSlots[param]] += g5[ord(param)]

proc fitMultiFrom(
    model: Model,
    doses, responses: openArray[float],
    curve: openArray[int],
    shared: set[DrcParam],
    weights: openArray[float],
    fixed: openArray[(DrcParam, float)],
    bounds: openArray[(DrcParam, Slice[float])],
    opts: FitOptions,
    start: seq[float],
    warm: bool,
): MultiFit =
  ## `fitMulti`, or with `warm`, one Levenberg–Marquardt run from the joint
  ## vector `start` in place of the starting-value search; for refits of data
  ## close to an earlier fit's.
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
    x = logDoses(doses)
    w = newSeq[float](doses.len)
    nCurve = newSeq[int](m)
    sumW, sumWY, sumLogW = 0.0
  for i in 0 ..< doses.len:
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
  result.data = (@doses, @responses, w)

  # Layout of the joint vector; fixed slots hold their value from the start.
  result.slots = newSeq[array[DrcParam, int]](m)
  var lower, upper: seq[float] # bounds on each slot
  for param in DrcParam:
    let copies = if param in result.shared: 1 else: m
    for k in 0 ..< copies:
      for c in 0 ..< m:
        if copies == 1 or c == k:
          result.slots[c][param] = result.joint.len
      result.joint.add fixedValues[ord(param)]
      result.jointFree.add param notin fixedParams
      lower.add box.lower[ord(param)]
      upper.add box.upper[ord(param)]
  var nFree = 0
  for f in result.jointFree:
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

  # Each curve's own data, for its starting values and its view.
  var curveData = newSeq[tuple[doses, responses, weights: seq[float]]](m)
  for i, c in curve:
    curveData[c].doses.add doses[i]
    curveData[c].responses.add responses[i]
    curveData[c].weights.add w[i]

  let
    obs = observationIndex(doses.len)
    ys = @responses
    joint = jointModel(x, result.slots, result.curveIds)

  proc search(starts: seq[seq[float]], free: seq[bool]): LMDynResult =
    for k, start in starts:
      let r = levenbergMarquardt(joint, obs, ys, w, start, free, lower, upper, opts.lm)
      if k == 0 or r.better(result):
        result = r

  var best: LMDynResult
  if warm:
    best = search(@[start], result.jointFree)
  else:
    # Starting values: each curve fitted on its own where it has enough data,
    # with shared parameters started at the typical value over curves.
    var singleOpts = opts
    singleOpts.fallback = false
    let curveFree = drcParamCount - card(fixedParams)
    var single = newSeq[ParamVec](m)
    for c in 0 ..< m:
      let (dc, rc, wc) = curveData[c]
      var gx, gy: seq[float]
      for k in 0 ..< dc.len:
        if wc[k] > 0.0:
          gx.add log10(dc[k])
          gy.add rc[k]
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

    var startA = result.joint
    for param in DrcParam:
      if param in fixedParams:
        continue
      if param in shared:
        var values: seq[float]
        for c in 0 ..< m:
          values.add single[c][ord(param)]
        startA[result.slots[0][param]] = clamp(
          sharedStart(param, values), box.lower[ord(param)], box.upper[ord(param)]
        )
      else:
        for c in 0 ..< m:
          startA[result.slots[c][param]] = single[c][ord(param)]
    var starts = @[startA]
    if shared - fixedParams != {} and result.n >= max(curveFree, 1):
      # The pooled curve is a good guess for the shared parameters.
      let pooled =
        fitCurve(model, doses, responses, w, fixed, bounds, singleOpts).params
      var startB = startA
      for param in shared - fixedParams:
        if pooled[param].isFinite:
          startB[result.slots[0][param]] = pooled[param]
      starts.add startB
    if opts.multiStart:
      # As the single-curve search: shift every LogEC50 by ±1 and scale every Hill.
      for shift in [0.0, -1.0, 1.0]:
        for scale in [1.0, 0.4, 2.5]:
          if (shift == 0.0 and scale == 1.0) or (
            shift != 0.0 and LogEC50 in fixedParams
          ) or (scale != 1.0 and Hill in fixedParams):
            continue
          var p = startA
          for param in [LogEC50, Hill]:
            var done: seq[int]
            for c in 0 ..< m:
              let slot = result.slots[c][param]
              if slot notin done:
                done.add slot
                if param == LogEC50:
                  p[slot] += shift
                else:
                  p[slot] *= scale
          starts.add p

    if Asymmetry notin fixedParams:
      # As for one curve: fit symmetric curves first, then release Asymmetry.
      var symFree = result.jointFree
      var symStarts = starts
      for c in 0 ..< m:
        symFree[result.slots[c][Asymmetry]] = false
      for s in symStarts.mitems:
        for c in 0 ..< m:
          s[result.slots[c][Asymmetry]] =
            clamp(1.0, box.lower[ord(Asymmetry)], box.upper[ord(Asymmetry)])
      best = search(symStarts, symFree)
      let released = search(@[best.params, startA], result.jointFree)
      if released.better(best):
        best = released
    else:
      best = search(starts, result.jointFree)

  result.joint = best.params
  result.iterations = best.iterations
  result.status = best.status
  result.converged = best.status.converged
  if not result.converged:
    result.warnings.incl fwNotConverged
  for c in 0 ..< m:
    for param in DrcParam:
      let slot = result.slots[c][param]
      if result.jointFree[slot] and result.joint[slot] in [lower[slot], upper[slot]]:
        result.atBound.incl param
  if result.atBound != {}:
    result.warnings.incl fwAtBound

  # Uncertainty, as for a single curve, over the joint vector, with parameters
  # at a bound counted as free (as R's nls with algorithm = "port").
  var
    jtwj, unscaled: seq[seq[float]]
    gradient: seq[float]
  result.ssr = normalEquations(joint, obs, responses, w, result.joint, jtwj, gradient)
  if not unscaledCovariance(jtwj, result.jointFree, unscaled):
    result.warnings.incl fwSingular
  result.df = result.n - nFree
  if result.df == 0:
    result.warnings.incl fwNoResidualDf
  result.sigma =
    if result.df > 0:
      sqrt(result.ssr / float(result.df))
    else:
      NaN
  let np = result.joint.len
  result.jointCov = newSeq[seq[float]](np)
  for i in 0 ..< np:
    result.jointCov[i] = newSeq[float](np)
    for j in 0 ..< np:
      if result.jointFree[i] and result.jointFree[j]:
        result.jointCov[i][j] = result.sigma ^ 2 * unscaled[i][j]

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
      v.params[param] = result.joint[result.slots[c][param]]
      for other in DrcParam:
        v.cov[param][other] =
          result.jointCov[result.slots[c][param]][result.slots[c][other]]
    for param in DrcParam:
      v.se[param] = sqrt(v.cov[param][param])
    v.n = nCurve[c]
    v.jointView = true
    (v.doses, v.responses, v.weights) = curveData[c]
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

proc fitMulti*(
    model: Model,
    doses, responses: openArray[float],
    curve: openArray[int],
    shared: set[DrcParam] = {},
    weights: openArray[float] = [],
    fixed: openArray[(DrcParam, float)] = [],
    bounds: openArray[(DrcParam, Slice[float])] = [],
    opts = FitOptions(),
): MultiFit =
  ## Fits `model` to several curves at once. Observation i belongs to curve
  ## `curve[i]` (ids 0 ..< m, each used). Parameters in `shared` take one
  ## value for all curves; the rest are fitted per curve. For example,
  ## `shared = {Bottom, Top, Hill}` fits parallel curves that differ only in
  ## EC50. Parameters in `fixed` are held at one value for all curves (a 3PL
  ## fixes Bottom at 0 unless `fixed` gives it), and those in `bounds` (as
  ## for `fitCurve`) are kept within one range for all curves.
  ##
  ## `fit.curves[c]` is curve c as a `DrcFit`, so `ec50`, `ecxCI`, `predict`
  ## and the bands work on it. Its covariance is curve c's block of the joint
  ## covariance, and its `df` and `sigma` are the joint fit's; its `ssr`,
  ## `rSquared` and `n` are its own. `param` and `cov` give the joint
  ## estimates, including covariances between curves. Only least squares is
  ## supported.
  fitMultiFrom(
    model, doses, responses, curve, shared, weights, fixed, bounds, opts, @[], false
  )

proc fitLL4Multi*(
    doses, responses: openArray[float],
    curve: openArray[int],
    shared: set[DrcParam] = {},
    weights: openArray[float] = [],
    fixed: openArray[(DrcParam, float)] = [],
    bounds: openArray[(DrcParam, Slice[float])] = [],
    opts = FitOptions(),
): MultiFit =
  ## `fitMulti` with four-parameter logistics.
  fitMulti(LL4, doses, responses, curve, shared, weights, fixed, bounds, opts)

proc fitLL3Multi*(
    doses, responses: openArray[float],
    curve: openArray[int],
    shared: set[DrcParam] = {},
    weights: openArray[float] = [],
    fixed: openArray[(DrcParam, float)] = [],
    bounds: openArray[(DrcParam, Slice[float])] = [],
    opts = FitOptions(),
): MultiFit =
  ## `fitMulti` with three-parameter logistics: Bottom is fixed, at 0 unless
  ## `fixed` gives it.
  fitMulti(LL3, doses, responses, curve, shared, weights, fixed, bounds, opts)

proc fitLL5Multi*(
    doses, responses: openArray[float],
    curve: openArray[int],
    shared: set[DrcParam] = {},
    weights: openArray[float] = [],
    fixed: openArray[(DrcParam, float)] = [],
    bounds: openArray[(DrcParam, Slice[float])] = [],
    opts = FitOptions(),
): MultiFit =
  ## `fitMulti` with five-parameter logistics.
  fitMulti(LL5, doses, responses, curve, shared, weights, fixed, bounds, opts)

proc checkCurve(fit: MultiFit, c: int) =
  if c notin 0 ..< fit.curves.len:
    raise newException(ValueError, "no such curve: " & $c)

proc param*(fit: MultiFit, curve: int, param: DrcParam): float =
  ## Curve `curve`'s estimate of `param` (shared parameters have one value
  ## for every curve).
  fit.checkCurve(curve)
  fit.joint[fit.slots[curve][param]]

proc cov*(fit: MultiFit, a: int, pa: DrcParam, b: int, pb: DrcParam): float =
  ## The covariance of curve `a`'s `pa` with curve `b`'s `pb`, from the joint
  ## covariance. It is 0 if either is fixed, and for curves a and b that
  ## share `pa` = `pb` it is that parameter's variance.
  fit.checkCurve(a)
  fit.checkCurve(b)
  fit.jointCov[fit.slots[a][pa]][fit.slots[b][pb]]

proc logPotency(fit: MultiFit, a, b: int, pct: float): tuple[center, variance: float] =
  ## log10(ECx of curve a / ECx of curve b) and its delta-method variance
  ## from the joint covariance.
  checkPercent(pct)
  fit.checkCurve(a)
  fit.checkCurve(b)
  let
    (logA, gradA) = fit.curves[a].logEcxAndGradient(pct)
    (logB, gradB) = fit.curves[b].logEcxAndGradient(pct)
  var g = newSeq[float](fit.joint.len)
  for param in DrcParam:
    g[fit.slots[a][param]] += gradA[param]
    g[fit.slots[b][param]] -= gradB[param]
  for i in 0 ..< g.len:
    if g[i] != 0.0:
      for j in 0 ..< g.len:
        if g[j] != 0.0:
          result.variance += g[i] * fit.jointCov[i][j] * g[j]
  result.center = logA - logB

proc relativePotency*(
    fit: MultiFit, a, b: int, pct = 50.0, level = 0.95
): tuple[ratio: float, ci: ConfInt] =
  ## ECx of curve `a` divided by ECx of curve `b`, with a confidence interval
  ## by the delta method on log10 of the ratio (so it is asymmetric), using
  ## the joint covariance including the parts shared by the two curves.
  ## R drc's `EDcomp` applies the delta method to the ratio itself by default,
  ## which gives a symmetric interval. See also `relativePotencyProfileCI`.
  let (center, variance) = fit.logPotency(a, b, pct)
  let half = fit.curves[a].tQuantile(level) * sqrt(variance)
  (pow(10.0, center), (pow(10.0, center - half), pow(10.0, center + half)))

proc slotBounds(fit: MultiFit): tuple[lower, upper: seq[float]] =
  result.lower = newSeq[float](fit.joint.len)
  result.upper = newSeq[float](fit.joint.len)
  for slots in fit.slots:
    for param in DrcParam:
      result.lower[slots[param]] = fit.lower[param]
      result.upper[slots[param]] = fit.upper[param]

proc relativePotencyProfileCI*(
    fit: MultiFit, a, b: int, pct = 50.0, level = 0.95, opts = FitOptions()
): ConfInt =
  ## Profile-likelihood confidence interval for `relativePotency(a, b, pct)`.
  ## The joint model is rewritten with log10(ECx) of curve `b` and
  ## log10(ECx(a)/ECx(b)) in place of the two curves' LogEC50s, and the log
  ## ratio is profiled as in `profileCI`, then back-transformed. An end is 0
  ## or Inf if the data do not bound the ratio within 10·t standard errors.
  ## Needs LogEC50 to be fitted separately for each curve: not shared, fixed
  ## or bounded. `opts.lm` controls the refits.
  if a == b:
    raise newException(ValueError, "the two curves must differ")
  let (center, variance) = fit.logPotency(a, b, pct)
  let tq = fit.curves[a].tQuantile(level)
  if LogEC50 in fit.shared:
    raise newException(
      ValueError, "LogEC50 is shared or fixed, so the potency ratio cannot be profiled"
    )
  if fit.lower[LogEC50] > -Inf or fit.upper[LogEC50] < Inf:
    raise newException(
      ValueError, "LogEC50 is bounded, so the potency ratio cannot be profiled"
    )

  # q is the joint vector with q[sb] = log10 ECx(b) and q[sa] = the log ratio,
  # so LogEC50(b) = q[sb] + z_b/Hill_b and LogEC50(a) = q[sa] + q[sb] + z_a/Hill_a
  # (see `zOffset`). Each curve's z depends only on its Asymmetry, so it is
  # kept between calls.
  type Shift = object
    r, cachedA, z, dzda: float
    hill, asymmetry: int # slots

  proc shift(curve: int): Shift =
    Shift(
      r: fractionAt(fit.param(curve, Hill), pct),
      cachedA: NaN,
      hill: fit.slots[curve][Hill],
      asymmetry: fit.slots[curve][Asymmetry],
    )

  var shifts = [shift(a), shift(b)]
  let
    (sa, sb) = (fit.slots[a][LogEC50], fit.slots[b][LogEC50])
    x = logDoses(fit.data.doses)
    slots = fit.slots
    curveIds = fit.curveIds
  let model = proc(i: float, q: seq[float], grad: var seq[float]): float =
    let k = int(i)
    let c = curveIds[k]
    let curveSlots = slots[c]
    var p, g: ParamVec
    for param in DrcParam:
      p[ord(param)] = q[curveSlots[param]]
    let side =
      if c == a:
        0
      elif c == b:
        1
      else:
        -1
    if side >= 0:
      let h = q[shifts[side].hill]
      let asym = q[shifts[side].asymmetry]
      if asym != shifts[side].cachedA:
        let o = zOffset(shifts[side].r, asym)
        (shifts[side].z, shifts[side].dzda, shifts[side].cachedA) = (o.z, o.dzda, asym)
      p[ord(LogEC50)] = q[sb] + shifts[side].z / h + (if side == 0: q[sa]
      else: 0.0)
    result = ll5(x[k], p, g)
    grad.setZero()
    for param in DrcParam:
      grad[curveSlots[param]] += g[ord(param)]
    if side >= 0:
      # Chain rule through the shifted LogEC50 (grad[sa] or grad[sb] is already ∂f/∂LogEC50).
      let dl = g[ord(LogEC50)]
      let h = q[shifts[side].hill]
      if side == 0:
        grad[sb] += dl
      grad[shifts[side].hill] -= dl * shifts[side].z / (h * h)
      grad[shifts[side].asymmetry] += dl * shifts[side].dzda / h
  var q = fit.joint
  q[sb] = fit.curves[b].logEcxAndGradient(pct)[0]
  q[sa] = center
  let (lower, upper) = fit.slotBounds
  let ci = profileInterval(
    model,
    observationIndex(x.len),
    fit.data.responses,
    fit.data.weights,
    q,
    sa,
    fit.jointFree,
    lower,
    upper,
    fit.ssr,
    fit.sigma,
    tq,
    sqrt(variance),
    false,
    opts.lm,
  )
  (pow(10.0, ci.lower), pow(10.0, ci.upper))

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

proc lackOfFit*(fit: MultiFit): FTest =
  ## Lack-of-fit F test of the joint fit against one free mean per dose of
  ## each curve, as `lackOfFit` for a single fit. A small p for parallel
  ## curves that `fTest` does not reject still means some curve is misfitted.
  var pe: tuple[ss: float, groups: int]
  for curve in fit.curves:
    let p = pureError(curve.doses, curve.responses, curve.weights)
    pe.ss += p.ss
    pe.groups += p.groups
  lackOfFitTest(fit.ssr, fit.df, fit.n, pe)

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
        slot = fit.slots[c][param]
        name =
          if param in fit.shared:
            $param
          else:
            &"{param}[{c}]"
        value = fit.joint[slot]
      if param in fit.fixed:
        result.add &"{name:<12}{value:>12.6g}{\"(fixed)\":>12}\n"
      else:
        let se = sqrt(fit.jointCov[slot][slot])
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
