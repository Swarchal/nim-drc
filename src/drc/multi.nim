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
  sigma*: float
    ## pooled residual standard error, sqrt(ssr / df); for a robust fit, the
    ## robust (MAD) scale of all the residuals
  rSquared*: float
  logLik*: float ## Gaussian log-likelihood, as R's logLik.nls (NaN for robust fits)
  aic*: float ## (NaN for robust fits)
  bic*: float ## (NaN for robust fits)
  loss*: Loss
  tuning*: float ## the loss's tuning constant (robust fits)
  robustWeights*: seq[float]
    ## final ψ(u)/u per observation, in the order given (robust fits; 0
    ## where the user's weight is 0)
  converged*: bool
    ## false if the optimiser or, for a robust fit, the reweighting did not
    ## converge (see `warnings`)
  usedFallback*: bool ## Levenberg–Marquardt failed and Nelder–Mead was run
  fallbackConverged*: bool
    ## with `usedFallback`: Nelder–Mead or its polishing run converged
  iterations*: int ## iterations of the final Levenberg–Marquardt run
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

type
  JointModel = proc(i: float, p: seq[float], grad: var seq[float]): float

  CurveData = tuple[doses, responses, weights: seq[float]] ## one curve's observations

  JointProblem = object
    ## A joint fit as the optimisers see it: least squares over the joint
    ## vector, with the observation index as x.
    model: JointModel
    obs: seq[float] ## 0, 1, …, n − 1
    y, w: seq[float]
    free: seq[bool] ## false for the slots of fixed parameters
    lower, upper: seq[float] ## the bounds on each slot

func observationIndex(n: int): seq[float] =
  ## 0, 1, …, n − 1: the x that LM sees for a joint fit.
  result = newSeq[float](n)
  for i in 0 ..< n:
    result[i] = float(i)

proc jointModel(
    x: seq[float], slots: seq[array[DrcParam, int]], curveIds: seq[int]
): JointModel =
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

iterator slotsOf(fit: MultiFit, param: DrcParam): int =
  ## Each slot of `param` in the joint vector: one if it is shared, else one
  ## per curve.
  for c in 0 ..< (if param in fit.shared: 1 else: fit.slots.len):
    yield fit.slots[c][param]

func freeCount(fit: MultiFit): int =
  ## The number of free joint parameters.
  for free in fit.jointFree:
    if free:
      inc result

proc slotBounds(fit: MultiFit): tuple[lower, upper: seq[float]] =
  result.lower = newSeq[float](fit.joint.len)
  result.upper = newSeq[float](fit.joint.len)
  for slots in fit.slots:
    for param in DrcParam:
      result.lower[slots[param]] = fit.lower[param]
      result.upper[slots[param]] = fit.upper[param]

proc layoutJoint(fit: var MultiFit, curves: int) =
  ## Lays out the joint vector for `curves` curves (`slots`, `joint` and
  ## `jointFree`) from `shared` and `fixed`. Fixed slots hold their value
  ## from the start, and the rest 0.
  fit.slots = newSeq[array[DrcParam, int]](curves)
  for param in DrcParam:
    let copies = if param in fit.shared: 1 else: curves
    for k in 0 ..< copies:
      for c in 0 ..< curves:
        if copies == 1 or c == k:
          fit.slots[c][param] = fit.joint.len
      fit.joint.add fit.fixedValues[ord(param)]
      fit.jointFree.add param notin fit.fixed

proc problem(fit: MultiFit): JointProblem =
  ## The least-squares problem of the fit's data and layout.
  result = JointProblem(
    model: jointModel(logDoses(fit.data.doses), fit.slots, fit.curveIds),
    obs: observationIndex(fit.data.doses.len),
    y: fit.data.responses,
    w: fit.data.weights,
    free: fit.jointFree,
  )
  (result.lower, result.upper) = fit.slotBounds

proc solve(
    problem: JointProblem, start: seq[float], free: seq[bool], opts: FitOptions
): LMDynResult =
  ## One Levenberg–Marquardt run from `start`, over the slots in `free`.
  levenbergMarquardt(
    problem.model, problem.obs, problem.y, problem.w, start, free, problem.lower,
    problem.upper, opts.lm,
  )

proc search(
    problem: JointProblem, starts: seq[seq[float]], free: seq[bool], opts: FitOptions
): LMDynResult =
  ## The best of the Levenberg–Marquardt runs from each of `starts`.
  for k, start in starts:
    let r = problem.solve(start, free, opts)
    if k == 0 or r.better(result):
      result = r

func nonZero(data: CurveData): int =
  ## The curve's observations with non-zero weight.
  for w in data.weights:
    if w > 0.0:
      inc result

proc curveStart(
    fit: MultiFit,
    data: CurveData,
    fixed: openArray[(DrcParam, float)],
    bounds: openArray[(DrcParam, Slice[float])],
    opts: FitOptions,
): ParamVec =
  ## Starting values for one curve: its own fit where it has enough data and
  ## that fit is finite, else the guess, within the bounds.
  if data.nonZero >= max(drcParamCount - card(fit.fixed), 1):
    result = fitCurve(
      fit.model, data.doses, data.responses, data.weights, fixed, bounds, opts
    ).params.toVec
    var finite = true
    for v in result:
      finite = finite and v.isFinite
    if finite:
      return
  result = guessStart(
    logDoses(data.doses), data.responses, data.weights, fit.fixedValues, fit.fixed
  )
  for param in DrcParam:
    result[ord(param)] = clamp(result[ord(param)], fit.lower[param], fit.upper[param])

proc gridStarts(fit: MultiFit, start: seq[float]): seq[seq[float]] =
  ## As the single-curve search: `start` with every LogEC50 shifted by ±1 and
  ## every Hill scaled by 0.4 or 2.5.
  for shift in [0.0, -1.0, 1.0]:
    for scale in [1.0, 0.4, 2.5]:
      if (shift == 0.0 and scale == 1.0) or (shift != 0.0 and LogEC50 in fit.fixed) or
          (scale != 1.0 and Hill in fit.fixed):
        continue
      var p = start
      for slot in fit.slotsOf(LogEC50):
        p[slot] += shift
      for slot in fit.slotsOf(Hill):
        p[slot] *= scale
      result.add p

proc jointStarts(
    fit: MultiFit,
    curveData: seq[CurveData],
    fixed: openArray[(DrcParam, float)],
    bounds: openArray[(DrcParam, Slice[float])],
    opts: FitOptions,
): seq[seq[float]] =
  ## The joint vectors to start the search from. The first has each curve's
  ## own starting values (see `curveStart`), with shared parameters at the
  ## typical value over curves. Then comes one with the shared parameters from
  ## the pooled curve and, with `opts.multiStart`, the grid of `gridStarts`.
  var singleOpts = opts
  singleOpts.fallback = false
  singleOpts.loss = Squared # a robust fit reweights from the least-squares one
  var single: seq[ParamVec]
  for data in curveData:
    single.add fit.curveStart(data, fixed, bounds, singleOpts)

  let sharedFree = fit.shared - fit.fixed
  var first = fit.joint
  for param in DrcParam:
    if param in fit.fixed:
      continue
    if param in fit.shared:
      var values: seq[float]
      for s in single:
        values.add s[ord(param)]
      first[fit.slots[0][param]] =
        clamp(sharedStart(param, values), fit.lower[param], fit.upper[param])
    else:
      for c, s in single:
        first[fit.slots[c][param]] = s[ord(param)]
  result.add first
  if sharedFree != {} and fit.n >= max(drcParamCount - card(fit.fixed), 1):
    # The pooled curve is a good guess for the shared parameters.
    let pooled = fitCurve(
      fit.model, fit.data.doses, fit.data.responses, fit.data.weights, fixed, bounds,
      singleOpts,
    ).params
    var second = first
    for param in sharedFree:
      if pooled[param].isFinite:
        second[fit.slots[0][param]] = pooled[param]
    result.add second
  if opts.multiStart:
    result.add fit.gridStarts(first)

proc searchJoint(
    problem: JointProblem, fit: MultiFit, starts: seq[seq[float]], opts: FitOptions
): LMDynResult =
  ## The best fit from `starts` (the main one first). As for one curve, a 5PL
  ## fits symmetric curves first and then releases Asymmetry.
  if Asymmetry in fit.fixed:
    return problem.search(starts, problem.free, opts)
  let symmetric = clamp(1.0, fit.lower[Asymmetry], fit.upper[Asymmetry])
  var
    symFree = problem.free
    symStarts = starts
  for slot in fit.slotsOf(Asymmetry):
    symFree[slot] = false
    for s in symStarts.mitems:
      s[slot] = symmetric
  result = problem.search(symStarts, symFree, opts)
  let released = problem.search(@[result.params, starts[0]], problem.free, opts)
  if released.better(result):
    result = released

proc fallback(
    problem: JointProblem, start: seq[float], opts: FitOptions
): tuple[fit: LMDynResult, converged: bool] =
  ## As for one curve: Nelder–Mead from `start`, then Levenberg–Marquardt to
  ## polish. It has converged if either stage has.
  let nm = nelderMead(
    proc(p: seq[float]): float =
      weightedSsr(problem.model, problem.obs, problem.y, problem.w, p),
    start,
    problem.free,
    problem.lower,
    problem.upper,
    opts.nelderMead,
  )
  result.fit = problem.solve(nm.params, problem.free, opts)
  result.converged = result.fit.status.converged or nm.converged

proc setUncertainty(
    fit: var MultiFit, problem: JointProblem, robustScale: float, u: seq[float]
) =
  ## `ssr`, `df`, `sigma` and the joint covariance at the final estimates, as
  ## for a single curve, with parameters at a bound counted as free (as R's
  ## nls with algorithm = "port"). For a robust fit, `robustScale` is its
  ## scale and `u` the residuals in units of it.
  let nFree = fit.freeCount
  var
    jtwj, unscaled: seq[seq[float]]
    gradient: seq[float]
  fit.ssr = normalEquations(
    problem.model, problem.obs, problem.y, problem.w, fit.joint, jtwj, gradient
  )
  if not unscaledCovariance(jtwj, fit.jointFree, unscaled):
    fit.warnings.incl fwSingular
  fit.df = fit.n - nFree
  if fit.df == 0:
    fit.warnings.incl fwNoResidualDf
  var varianceScale: float
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
      fit.loss, fit.tuning, u, problem.w, robustScale, nFree, fit.n, fit.df
    )
  let np = fit.joint.len
  fit.jointCov = newSeq[seq[float]](np)
  for i in 0 ..< np:
    fit.jointCov[i] = newSeq[float](np)
    for j in 0 ..< np:
      if fit.jointFree[i] and fit.jointFree[j]:
        fit.jointCov[i][j] = varianceScale * unscaled[i][j]

proc setGoodnessOfFit(fit: var MultiFit, sst: float) =
  ## `rSquared` from `ssr` and the total sum of squares `sst`, and the
  ## likelihood criteria (NaN for a robust fit).
  fit.rSquared = 1.0 - fit.ssr / sst
  (fit.logLik, fit.aic, fit.bic) =
    if fit.loss == Squared:
      likelihoodCriteria(fit.ssr, fit.n, fit.freeCount, fit.data.weights)
    else:
      (NaN, NaN, NaN)

proc curveView(fit: MultiFit, c: int, data: CurveData): DrcFit =
  ## Curve `c` of the finished joint fit as a `DrcFit`: its parameters and
  ## covariance block, with the joint σ and df so intervals and bands account
  ## for the pooled fit. Its n, SSR, R², dose range and bound and EC50
  ## warnings are its own.
  result.model = fit.model
  result.fixed = fit.fixed
  (result.lower, result.upper) = (fit.lower, fit.upper)
  result.loss = fit.loss
  result.tuning = fit.tuning
  if fit.loss != Squared:
    for i, id in fit.curveIds:
      if id == c:
        result.robustWeights.add fit.robustWeights[i]
  for param in DrcParam:
    result.params[param] = fit.joint[fit.slots[c][param]]
    for other in DrcParam:
      result.cov[param][other] = fit.jointCov[fit.slots[c][param]][fit.slots[c][other]]
  for param in DrcParam:
    result.se[param] = sqrt(result.cov[param][param])
  result.jointView = true
  (result.doses, result.responses, result.weights) = data
  result.df = fit.df
  result.sigma = fit.sigma

  result.n = data.nonZero
  result.logDoseRange = nonZeroLogDoseRange(data.doses, data.weights)
  let p = result.params.toVec
  var grad: ParamVec
  for i, x in logDoses(data.doses):
    if data.weights[i] > 0.0:
      result.ssr += data.weights[i] * (data.responses[i] - ll5(x, p, grad)) ^ 2
  result.rSquared = 1.0 - result.ssr / totalSumOfSquares(data.responses, data.weights)
  (result.logLik, result.aic, result.bic) = (NaN, NaN, NaN)

  result.converged = fit.converged
  result.usedFallback = fit.usedFallback
  result.fallbackConverged = fit.fallbackConverged
  result.iterations = fit.iterations
  result.status = fit.status
  result.warnings = fit.warnings - {fwAtBound, fwEc50OutsideRange}
  result.flagBounds
  result.flagEc50Range

proc curveCount(curve: openArray[int]): int =
  ## The number of curves: one more than the largest id.
  for c in curve:
    if c < 0:
      raise newException(ValueError, "curve ids must be non-negative, got " & $c)
    result = max(result, c + 1)

proc setData(
    fit: var MultiFit,
    doses, responses: openArray[float],
    curve: openArray[int],
    weights: openArray[float],
): seq[CurveData] =
  ## Checks the observations and stores them in `fit`, with weights of 1 if
  ## none are given (`data`, `dataHash`, `curveIds`, `n`). Returns each
  ## curve's own observations.
  result = newSeq[CurveData](curveCount(curve))
  var w = newSeq[float](doses.len)
  for i, c in curve:
    w[i] =
      if weights.len == 0:
        1.0
      else:
        weights[i]
    checkObservation(doses[i], responses[i], w[i])
    if w[i] > 0.0:
      inc fit.n
    result[c].doses.add doses[i]
    result[c].responses.add responses[i]
    result[c].weights.add w[i]
  for c, data in result:
    if data.nonZero == 0:
      raise newException(
        ValueError, "curve " & $c & " has no observations with non-zero weight"
      )
  fit.curveIds = @curve
  fit.dataHash = !$(hash(@doses) !& hash(@responses) !& hash(w))
  fit.data = (@doses, @responses, w)

proc checkCounts(fit: MultiFit, curveData: seq[CurveData]) =
  ## Raises unless there are as many observations with non-zero weight as
  ## free joint parameters, and in each curve as many as its unshared ones.
  let
    nFree = fit.freeCount
    unshared = drcParamCount - card(fit.shared)
  if fit.n < max(nFree, 1):
    raise newException(
      ValueError,
      "need at least " & $max(nFree, 1) & " observations with " &
        "non-zero weight to fit " & $nFree & " parameters, got " & $fit.n,
    )
  for c, data in curveData:
    if data.nonZero < unshared:
      raise newException(
        ValueError,
        "curve " & $c & " has " & $data.nonZero &
          " observations with non-zero weight, fewer than its " & $unshared &
          " unshared parameters",
      )

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
  result.loss = opts.loss
  result.tuning = resolveTuning(opts)

  let curveData = result.setData(doses, responses, curve, weights)
  result.layoutJoint(curveData.len)
  result.checkCounts(curveData)

  let
    problem = result.problem
    starts =
      if warm:
        @[start]
      else:
        result.jointStarts(curveData, fixed, bounds, opts)
  var best =
    if warm:
      problem.solve(start, problem.free, opts)
    else:
      problem.searchJoint(result, starts, opts)
  result.converged = best.status.converged
  if opts.fallback and not result.converged:
    let fb = problem.fallback(
      if best.ssr.isFinite:
        best.params
      else:
        starts[0],
      opts,
    )
    best = fb.fit
    result.usedFallback = true
    result.fallbackConverged = fb.converged
    result.converged = fb.converged

  # A robust fit reweights from the least-squares fit, as for one curve.
  let sst = totalSumOfSquares(problem.y, problem.w)
  var
    robust: IrlsResult[LMDynResult]
    u: seq[float] # residuals in units of the robust scale
  if result.loss != Squared:
    robust = robustFit(
      problem.model, problem.obs, problem.y, problem.w, best, problem.free,
      problem.lower, problem.upper, result.loss, result.tuning, sst, result.n, opts,
    )
    best = robust.fit
    result.converged = best.status.converged
    if not robust.converged:
      result.warnings.incl fwRobustNotConverged
    (u, result.robustWeights) =
      robust.robustUnits(problem.w, result.loss, result.tuning)

  result.joint = best.params
  result.iterations = best.iterations
  result.status = best.status
  if not result.converged:
    result.warnings.incl fwNotConverged
  result.converged = result.warnings * {fwNotConverged, fwRobustNotConverged} == {}
  for param in DrcParam:
    if param notin result.fixed:
      for slot in result.slotsOf(param):
        if result.joint[slot] in [result.lower[param], result.upper[param]]:
          result.atBound.incl param
  if result.atBound != {}:
    result.warnings.incl fwAtBound

  result.setUncertainty(problem, robust.scale, u)
  result.setGoodnessOfFit(sst)
  for c, data in curveData:
    result.curves.add result.curveView(c, data)
  for view in result.curves:
    if fwEc50OutsideRange in view.warnings:
      result.warnings.incl fwEc50OutsideRange

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
  ## estimates, including covariances between curves.
  ##
  ## As for `fitCurve`, `opts.loss` = Huber or Tukey gives a robust fit, and
  ## Nelder–Mead takes over if Levenberg–Marquardt does not converge.
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

proc checkProfile(fit: MultiFit) =
  if fit.loss != Squared:
    raise
      newException(ValueError, "profile-likelihood intervals need a least-squares fit")

proc pinnedTau(
    fit: MultiFit, problem: JointProblem, slot: int, opts: FitOptions
): ProfileTau[seq[float]] =
  ## τ for slot `slot` of `problem`'s parameter vector: each refit holds it
  ## at β. `problem` is the fit's, or has a model that fits equally well.
  pinnedTau(
    problem.model, problem.obs, problem.y, problem.w, slot, problem.free, problem.lower,
    problem.upper, fit.ssr, fit.sigma, opts.lm,
  )

proc constrainedInterval[G](
    fit: MultiFit, g: G, tq, se: float, opts: FitOptions
): ConfInt =
  ## The profile-likelihood interval for g(joint vector), whose delta-method
  ## standard error is `se`: each refit of the joint model is constrained to
  ## g (see `constrainedInterval` for one curve).
  let problem = fit.problem
  constrainedInterval(
    problem.model, g, problem.obs, problem.y, problem.w, fit.joint, problem.free,
    problem.lower, problem.upper, fit.ssr, fit.sigma, tq, se, opts.lm,
  )

proc profileCI*(
    fit: MultiFit, curve: int, param: DrcParam, level = 0.95, opts = FitOptions()
): ConfInt =
  ## Profile-likelihood confidence interval for curve `curve`'s `param`, as
  ## `profileCI` for a single fit, refitting the whole joint model (a shared
  ## parameter has the same interval for every curve). Least-squares fits
  ## only.
  fit.checkCurve(curve)
  fit.checkProfile
  let tq = fit.curves[curve].tQuantile(level)
  let slot = fit.slots[curve][param]
  if not fit.jointFree[slot]:
    return (fit.joint[slot], fit.joint[slot])
  profileInterval(
    fit.pinnedTau(fit.problem, slot, opts),
    fit.joint,
    fit.joint[slot],
    fit.lower[param],
    fit.upper[param],
    fit.sigma,
    tq,
    sqrt(fit.jointCov[slot][slot]),
    param == Asymmetry,
  )

proc ecxProfileCI*(
    fit: MultiFit, curve: int, pct: float, level = 0.95, opts = FitOptions()
): ConfInt =
  ## Profile-likelihood confidence interval for curve `curve`'s `ecx`, as
  ## `ecxProfileCI` for a single fit: each refit of the joint model is
  ## constrained to the curve's log10(ECx). This works whether LogEC50 is
  ## separate, shared, fixed or bounded. If ECx depends on no free
  ## parameter, the interval is the estimate. Least-squares fits only.
  checkPercent(pct)
  fit.checkCurve(curve)
  fit.checkProfile
  let view = fit.curves[curve]
  let tq = view.tQuantile(level)
  let (_, grad) = view.logEcxAndGradient(pct)
  let ci = fit.constrainedInterval(
    logEcxConstraint[seq[float]](fit.slots[curve], fractionAt(view.params[Hill], pct)),
    tq,
    sqrt(view.deltaVariance(grad, LogEC50 .. Asymmetry)),
    opts,
  )
  antilog(ci)

proc ec50ProfileCI*(
    fit: MultiFit, curve: int, level = 0.95, opts = FitOptions()
): ConfInt =
  fit.ecxProfileCI(curve, 50.0, level, opts)

proc logPotencyConstraint(
    fit: MultiFit, a, b: int, pct: float
): proc(q: seq[float], grad: var seq[float]): float =
  ## g(joint vector) = log10(ECx of curve a / ECx of curve b), with its
  ## gradient.
  let
    ga = logEcxConstraint[seq[float]](fit.slots[a], fractionAt(fit.param(a, Hill), pct))
    gb = logEcxConstraint[seq[float]](fit.slots[b], fractionAt(fit.param(b, Hill), pct))
  result = proc(q: seq[float], grad: var seq[float]): float =
    var gradB = newSeq[float](q.len)
    result = ga(q, grad) - gb(q, gradB)
    for i in 0 ..< q.len:
      grad[i] -= gradB[i]

proc potencyModel(fit: MultiFit, a, b: int, pct: float): JointModel =
  ## The joint model rewritten for the potency of curve a relative to curve b
  ## at `pct`: in its parameter vector q, curve b's LogEC50 slot (sb) holds
  ## log10 ECx(b) and curve a's (sa) holds log10(ECx(a)/ECx(b)), so
  ## LogEC50(b) = q[sb] + z_b/Hill_b and LogEC50(a) = q[sa] + q[sb] + z_a/Hill_a
  ## (see `zOffset`). The curves must not share LogEC50.
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

  # Each curve's z depends only on its Asymmetry, so it is kept between calls.
  var shifts = [shift(a), shift(b)]
  let
    (sa, sb) = (fit.slots[a][LogEC50], fit.slots[b][LogEC50])
    x = logDoses(fit.data.doses)
    slots = fit.slots
    curveIds = fit.curveIds
  result = proc(i: float, q: seq[float], grad: var seq[float]): float =
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

proc relativePotencyProfileCI*(
    fit: MultiFit, a, b: int, pct = 50.0, level = 0.95, opts = FitOptions()
): ConfInt =
  ## Profile-likelihood confidence interval for `relativePotency(a, b, pct)`.
  ## The joint model is rewritten with log10(ECx) of curve `b` and
  ## log10(ECx(a)/ECx(b)) in place of the two curves' LogEC50s, and the log
  ## ratio is profiled as in `profileCI`, then back-transformed. An end is 0
  ## or Inf if the data do not bound the ratio within 10·t standard errors.
  ## If LogEC50 is shared, fixed or bounded, each refit is instead
  ## constrained to the log ratio. If the ratio depends on no free parameter
  ## (the EC50s of 4PL curves that share LogEC50), the interval is the
  ## estimate. Least-squares fits only. `opts.lm` controls the refits.
  if a == b:
    raise newException(ValueError, "the two curves must differ")
  let (center, variance) = fit.logPotency(a, b, pct)
  let tq = fit.curves[a].tQuantile(level)
  fit.checkProfile
  if LogEC50 in fit.shared or fit.lower[LogEC50] > -Inf or fit.upper[LogEC50] < Inf:
    return antilog(
      fit.constrainedInterval(
        fit.logPotencyConstraint(a, b, pct), tq, sqrt(variance), opts
      )
    )
  let (sa, sb) = (fit.slots[a][LogEC50], fit.slots[b][LogEC50])
  var problem = fit.problem
  problem.model = fit.potencyModel(a, b, pct)
  var q = fit.joint
  q[sb] = fit.curves[b].logEcxAndGradient(pct)[0]
  q[sa] = center
  antilog(
    profileInterval(
      fit.pinnedTau(problem, sa, opts),
      q,
      center,
      -Inf,
      Inf,
      fit.sigma,
      tq,
      sqrt(variance),
    )
  )

proc fTest*(reduced, full: MultiFit): FTest =
  ## F test of a reduced joint fit (e.g. parallel curves) against a fuller one
  ## (e.g. separate curves) on the same data, as R's `anova` for nested
  ## models. A small p means the extra parameters improve the fit.
  ## Least-squares fits only.
  if reduced.loss != Squared or full.loss != Squared:
    raise newException(ValueError, "an F test needs least-squares fits")
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
  ## Least-squares fits only.
  if fit.loss != Squared:
    raise newException(ValueError, "a lack-of-fit test needs a least-squares fit")
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
  if fit.loss != Squared:
    result.add &", robust ({fit.loss}, tuning {fit.tuning:g})"
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
  if fit.loss == Squared:
    result.add &"\nResidual standard error: {fit.sigma:.6g} on {fit.df} degrees of freedom\n"
    result.add &"R²: {fit.rSquared:.6g}   AIC: {fit.aic:.6g}   BIC: {fit.bic:.6g}\n"
  else:
    result.add &"\nRobust residual scale (MAD): {fit.sigma:.6g} on {fit.df} degrees of freedom\n"
    result.add robustWeightNote(fit.robustWeights, fit.n) & "\n"
    result.add &"R²: {fit.rSquared:.6g}\n"
  result.add fallbackNote(
    fit.loss, fit.usedFallback, fit.fallbackConverged, fit.status, fit.iterations
  )
  for warning in fit.warnings:
    result.add &"\nWarning: {warning}"
