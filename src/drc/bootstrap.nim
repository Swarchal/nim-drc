## Bootstrap intervals for single and joint fits: refit resampled data many
## times and take percentiles of any statistic over the refits.
##
## This file is included by drc.nim (after multi.nim), not imported, so it
## can use their private helpers.

import std/random

type
  Resample* = enum
    Residuals = "residuals"
      ## new responses: the fitted values plus residuals drawn with replacement
    Cases = "cases" ## observations drawn with replacement (within each curve)

  Bootstrap*[T] = object
    fit*: T ## the original fit
    replicates*: seq[T]
      ## the refits that converged; they do not keep their data
      ## (`doses`, `responses`, `weights`), so `residuals` and the tests
      ## need the original fit
    failed*: int ## resamples whose refit did not converge, left out
    resample*: Resample

proc refitArgs(
    model: Model,
    fixed: set[DrcParam],
    value: proc(param: DrcParam): float,
    lower, upper: array[DrcParam, float],
): tuple[fixed: seq[(DrcParam, float)], bounds: seq[(DrcParam, Slice[float])]] =
  ## The `fixed` and `bounds` arguments that reproduce a fit's settings.
  for param in fixed:
    if param != Asymmetry or model == LL5: # 3PL and 4PL fix it themselves
      result.fixed.add (param, value(param))
  for param in DrcParam:
    if param notin fixed and (lower[param] > -Inf or upper[param] < Inf):
      result.bounds.add (param, lower[param] .. upper[param])

proc residualPool(
    fitted, responses, weights: openArray[float], n, df: int
): seq[float] =
  ## The residuals √wᵢ·(yᵢ − f(xᵢ)) of observations with non-zero weight,
  ## centred and scaled by sqrt(n/df) so their spread matches σ (Davison and
  ## Hinkley's modified residuals, simplified for a nonlinear fit).
  if df < 1:
    raise
      newException(ValueError, "a residual bootstrap needs residual degrees of freedom")
  var total = 0.0
  for i in 0 ..< fitted.len:
    if weights[i] > 0.0:
      result.add sqrt(weights[i]) * (responses[i] - fitted[i])
      total += result[^1]
  let mean = total / float(result.len)
  let scale = sqrt(float(n) / float(df))
  for e in result.mitems:
    e = (e - mean) * scale

proc checkReplicates(n: int) =
  if n < 1:
    raise newException(ValueError, "need at least 1 replicate, got " & $n)

proc stripData(fit: var DrcFit) =
  fit.doses = @[]
  fit.responses = @[]
  fit.weights = @[]
  fit.robustWeights = @[]

proc bootstrap*(
    fit: DrcFit, n = 2000, resample = Residuals, seed = 1, opts = FitOptions()
): Bootstrap[DrcFit] =
  ## Refits `fit`'s model to `n` resampled data sets, with the same fixed
  ## parameters, bounds and loss (`opts` gives the other settings). With
  ## `Residuals`, the responses are the fitted values plus residuals
  ## √wⱼ·rⱼ drawn with replacement (centred and scaled by sqrt(n/df)),
  ## divided by √wᵢ; with `Cases`, observations are drawn with replacement.
  ## Zero-weight observations are left as they are. Resampling is
  ## reproducible for a given `seed`. Use `bootCI` for intervals.
  ##
  ## Residual resampling assumes the model and the weights are right; case
  ## resampling does not, but needs more data, since a resample can miss
  ## doses. Not available for the curves of a joint fit (bootstrap the
  ## `MultiFit` instead).
  checkReplicates(n)
  if fit.jointView:
    raise newException(
      ValueError, "for a curve of a joint fit, bootstrap the joint fit instead"
    )
  result.fit = fit
  result.resample = resample
  var refitOpts = opts
  refitOpts.loss = fit.loss
  if fit.loss != Squared:
    refitOpts.tuning = fit.tuning
  let (fixed, bounds) = refitArgs(
    fit.model,
    fit.fixed,
    proc(param: DrcParam): float =
      fit.params[param],
    fit.lower,
    fit.upper,
  )
  var rng = initRand(seed)
  let fitted = fit.fitted
  var pool: seq[float] # residuals, or the observations to draw from
  var cases: seq[int]
  if resample == Residuals:
    pool = residualPool(fitted, fit.responses, fit.weights, fit.n, fit.df)
  else:
    for i in 0 ..< fit.doses.len:
      if fit.weights[i] > 0.0:
        cases.add i
  for _ in 0 ..< n:
    var doses, responses, weights: seq[float]
    if resample == Residuals:
      (doses, responses, weights) = (fit.doses, fit.responses, fit.weights)
      for i in 0 ..< doses.len:
        if weights[i] > 0.0:
          responses[i] = fitted[i] + pool[rng.rand(pool.high)] / sqrt(weights[i])
    else:
      for _ in cases:
        let i = cases[rng.rand(cases.high)]
        doses.add fit.doses[i]
        responses.add fit.responses[i]
        weights.add fit.weights[i]
    try:
      var r = fitCurve(fit.model, doses, responses, weights, fixed, bounds, refitOpts)
      if r.converged:
        r.stripData
        result.replicates.add r
      else:
        inc result.failed
    except ValueError:
      inc result.failed

proc bootstrap*(
    fit: MultiFit, n = 2000, resample = Residuals, seed = 1, opts = FitOptions()
): Bootstrap[MultiFit] =
  ## As `bootstrap` for a single fit, refitting the joint model with the same
  ## shared and fixed parameters and bounds. Residuals are drawn from all
  ## curves together, since the joint fit has one σ; cases are drawn within
  ## each curve, so every curve keeps its number of observations. For
  ## example, an interval for the potency of curve 1 relative to curve 0:
  ##
  ## .. code-block:: nim
  ##   let boot = parallel.bootstrap()
  ##   echo boot.bootCI(proc(f: MultiFit): float = f.relativePotency(1, 0).ratio)
  checkReplicates(n)
  result.fit = fit
  result.resample = resample
  let (fixed, bounds) = refitArgs(
    fit.model,
    fit.fixed,
    proc(param: DrcParam): float =
      fit.fixedValues[ord(param)],
    fit.lower,
    fit.upper,
  )
  let shared = fit.shared - fit.fixed
  let data = fit.jointData
  var fitted = newSeq[float](data.doses.len)
  for i, c in fit.curveIds:
    fitted[i] = fit.curves[c].predict(data.doses[i])
  var rng = initRand(seed)
  var pool: seq[float]
  var cases = newSeq[seq[int]](fit.curves.len) # each curve's observations
  if resample == Residuals:
    pool = residualPool(fitted, data.responses, data.weights, fit.n, fit.df)
  else:
    for i, c in fit.curveIds:
      if data.weights[i] > 0.0:
        cases[c].add i
  for _ in 0 ..< n:
    var doses, responses, weights: seq[float]
    var curve: seq[int]
    if resample == Residuals:
      (doses, responses, weights) = (data.doses, data.responses, data.weights)
      curve = fit.curveIds
      for i in 0 ..< doses.len:
        if weights[i] > 0.0:
          responses[i] = fitted[i] + pool[rng.rand(pool.high)] / sqrt(weights[i])
    else:
      for c, rows in cases:
        for _ in rows:
          let i = rows[rng.rand(rows.high)]
          doses.add data.doses[i]
          responses.add data.responses[i]
          weights.add data.weights[i]
          curve.add c
    try:
      var r = fitMulti(
        fit.model, doses, responses, curve, shared, weights, fixed, bounds, opts
      )
      if r.converged:
        for c in r.curves.mitems:
          c.stripData
        result.replicates.add r
      else:
        inc result.failed
    except ValueError:
      inc result.failed

proc bootValues*[T](boot: Bootstrap[T], statistic: proc(fit: T): float): seq[float] =
  ## `statistic` of each replicate, e.g. for a histogram.
  for r in boot.replicates:
    result.add statistic(r)

proc percentile(sorted: openArray[float], alpha: float): float =
  ## The alpha quantile of sorted bootstrap values as R's `boot.ci` (its
  ## `norm.inter`): order statistic (R + 1)·alpha, interpolated on the normal
  ## quantile scale between neighbours, or the extreme value beyond them.
  let
    r = sorted.len
    rk = float(r + 1) * alpha
    k = int(rk)
  if float(k) == rk and k >= 1 and k <= r:
    sorted[k - 1]
  elif k <= 0:
    sorted[0]
  elif k >= r:
    sorted[r - 1]
  else:
    let
      z = normalQuantile(alpha)
      zk = normalQuantile(float(k) / float(r + 1))
      zk1 = normalQuantile(float(k + 1) / float(r + 1))
    sorted[k - 1] + (z - zk) / (zk1 - zk) * (sorted[k] - sorted[k - 1])

proc bootCI*[T](
    boot: Bootstrap[T], statistic: proc(fit: T): float, level = 0.95
): ConfInt =
  ## Percentile bootstrap confidence interval for `statistic`, any function
  ## of a fit: for example `boot.bootCI(ec50)`, or
  ## `boot.bootCI(proc(f: DrcFit): float = f.doseAt(50.0))`. The ends are
  ## the (1 ∓ level)/2 quantiles of the statistic over the replicates, as R's
  ## `boot.ci(type = "perc")`, so the interval respects any transformation
  ## (e.g. log dose). Replicates where the statistic is not finite are left
  ## out (see `bootValues`). With fewer than about 1000 replicates the ends
  ## of a 95% interval are noisy. NaN without any finite values.
  if not (level > 0.0 and level < 1.0):
    raise newException(ValueError, "confidence level must be in (0, 1), got " & $level)
  var values: seq[float]
  for v in boot.bootValues(statistic):
    if v.isFinite:
      values.add v
  if values.len == 0:
    return (NaN, NaN)
  values.sort()
  (values.percentile(0.5 - 0.5 * level), values.percentile(0.5 + 0.5 * level))
