## Levenberg–Marquardt for weighted nonlinear least squares:
##
##   minimise  SSR(p) = Σ wᵢ (yᵢ − f(xᵢ, p))²
##
## Each iteration solves the damped normal equations
##
##   (JᵀWJ + λ·diag(JᵀWJ)) δ = JᵀW r
##
## and adapts λ: smaller after a successful step (towards Gauss–Newton),
## larger after a failed one (towards scaled gradient descent).

import std/math
import linalg

type
  LMStatus* = enum
    lmSmallSsrChange = "relative change in SSR below ftol"
    lmSmallStep = "relative step size below xtol"
    lmStationary = "no downhill step found; the gradient is zero to working precision"
    lmStalled = "no downhill step found, but the gradient is not zero"
    lmMaxIter = "maximum iterations reached"
    lmBadStart = "SSR is not finite at the starting parameters"

  LMOptions* = object
    maxIter*: int = 200
    ftol*: float = 1e-12    ## stop when the relative SSR decrease falls below this
    xtol*: float = 1e-10    ## stop when every relative parameter step falls below this
    gtol*: float = 1e-6     ## a stuck fit is stationary if |cos(residuals, each Jacobian column)| ≤ gtol
    lambda0*: float = 1e-3  ## initial damping
    lambdaMax*: float = 1e16

  LMResult*[K: static int] = object
    params*: array[K, float]
    ssr*: float
    iterations*: int
    status*: LMStatus

func converged*(status: LMStatus): bool =
  status in {lmSmallSsrChange, lmSmallStep, lmStationary}

proc weightedSsr*[K: static int; F](fn: F; x, y, w: openArray[float]; p: array[K, float]): float =
  ## Σ wᵢ (yᵢ − f(xᵢ, p))²
  var grad: array[K, float]
  for i in 0 ..< x.len:
    let r = y[i] - fn(x[i], p, grad)
    result += w[i] * r * r

proc normalEquations*[K: static int; F](fn: F; x, y, w: openArray[float]; p: array[K, float];
                                        a: var array[K, array[K, float]]; g: var array[K, float]): float =
  ## Fills `a = JᵀWJ` and `g = JᵀWr` at `p`, and returns the SSR.
  a = default(array[K, array[K, float]])
  g = default(array[K, float])
  var grad: array[K, float]
  for i in 0 ..< x.len:
    let r = y[i] - fn(x[i], p, grad)
    result += w[i] * r * r
    for j in 0 ..< K:
      let wj = w[i] * grad[j]
      g[j] += wj * r
      for k in 0 .. j:
        a[j][k] += wj * grad[k]
  for j in 0 ..< K:
    for k in j + 1 ..< K:
      a[j][k] = a[k][j]

func stationary[K: static int](a: array[K, array[K, float]]; g: array[K, float]; ssr: float;
                               free: array[K, bool]; perfectSsr, gtol: float): bool =
  ## True at a minimum to working precision: an essentially perfect fit, or
  ## residuals orthogonal to every free Jacobian column, |gⱼ| ≤ gtol·‖Jⱼ‖·‖r‖.
  if ssr <= perfectSsr: return true
  for j in 0 ..< K:
    if free[j] and not (abs(g[j]) <= gtol * sqrt(a[j][j] * ssr)):  # NaN is not stationary
      return false
  true

proc levenbergMarquardt*[K: static int; F](fn: F; x, y, w: openArray[float];
                                           p0: array[K, float]; free: array[K, bool];
                                           opts = LMOptions()): LMResult[K] =
  ## Fits the parameters marked in `free`, holding the others at their value in `p0`.
  ## `fn(x: float; p: array[K, float]; grad: var array[K, float]): float` returns f(x, p) and
  ## writes ∂f/∂pⱼ into `grad`.
  assert x.len == y.len and x.len == w.len
  var
    p = p0
    a: array[K, array[K, float]]
    g: array[K, float]
    ssr = normalEquations(fn, x, y, w, p, a, g)
    lambda = opts.lambda0
    sumSquares = 0.0
  for i in 0 ..< y.len: sumSquares += w[i] * y[i] * y[i]
  # Residuals around 1e-10 of the data are rounding error: a perfect fit.
  let perfectSsr = 1e-20 * sumSquares

  result.status = if ssr.isFinite: lmMaxIter else: lmBadStart
  while result.status == lmMaxIter and result.iterations < opts.maxIter:
    inc result.iterations
    var maxDiag = 0.0
    for j in 0 ..< K:
      if free[j]: maxDiag = max(maxDiag, a[j][j])
    let diagFloor = max(maxDiag * 1e-12, 1e-300)
    var rhs = g
    for j in 0 ..< K:
      if not free[j]: rhs[j] = 0.0

    # Raise the damping until a step reduces the SSR. The first trial is
    # usually accepted, so it builds the next normal equations as it goes;
    # after a rejection, trials only compute the SSR.
    var
      accepted = false
      firstTrial = true
      aNew: array[K, array[K, float]]
      gNew: array[K, float]
    while not accepted and lambda <= opts.lambdaMax:
      var
        m = a
        delta: array[K, float]
      for j in 0 ..< K:
        if free[j]: m[j][j] += lambda * max(a[j][j], diagFloor)
      pinFixed(m, free)
      if solveSpd(m, rhs, delta):
        var pNew = p
        for j in 0 ..< K: pNew[j] += delta[j]
        let full = firstTrial
        firstTrial = false
        let ssrNew =
          if full: normalEquations(fn, x, y, w, pNew, aNew, gNew)
          else: weightedSsr(fn, x, y, w, pNew)
        if ssrNew.isFinite and ssrNew < ssr:
          accepted = true
          var smallStep = true
          for j in 0 ..< K:
            if abs(delta[j]) > opts.xtol * (abs(p[j]) + opts.xtol):
              smallStep = false
          if ssr - ssrNew <= opts.ftol * ssr:
            result.status = lmSmallSsrChange
          elif smallStep:
            result.status = lmSmallStep
          p = pNew
          if full: (ssr, a, g) = (ssrNew, aNew, gNew)
          else: ssr = normalEquations(fn, x, y, w, p, a, g)
          lambda = max(lambda / 10.0, 1e-15)
          continue
      lambda *= 10.0

    if not accepted:
      result.status =
        if stationary(a, g, ssr, free, perfectSsr, opts.gtol): lmStationary else: lmStalled

  result.params = p
  result.ssr = ssr
