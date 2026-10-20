## Distribution functions and least-squares fit statistics.

import std/math
import linalg

func betaContinuedFraction(a, b, x: float): float =
  ## Continued fraction for the incomplete beta function, by the modified
  ## Lentz method (Numerical Recipes §6.4).
  const
    maxIter = 1000
    eps = 1e-16
    tiny = 1e-300
  let
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
  var
    c = 1.0
    d = 1.0 - qab * x / qap
  if abs(d) < tiny: d = tiny
  d = 1.0 / d
  result = d
  for m in 1 .. maxIter:
    let
      fm = float(m)
      m2 = 2.0 * fm
    # Even step.
    var aa = fm * (b - fm) * x / ((qam + m2) * (a + m2))
    d = 1.0 + aa * d
    if abs(d) < tiny: d = tiny
    c = 1.0 + aa / c
    if abs(c) < tiny: c = tiny
    d = 1.0 / d
    result *= d * c
    # Odd step.
    aa = -(a + fm) * (qab + fm) * x / ((a + m2) * (qap + m2))
    d = 1.0 + aa * d
    if abs(d) < tiny: d = tiny
    c = 1.0 + aa / c
    if abs(c) < tiny: c = tiny
    d = 1.0 / d
    let delta = d * c
    result *= delta
    if abs(delta - 1.0) < eps:
      break

func incompleteBeta*(a, b, x: float; y = 1.0 - x): float =
  ## Regularised incomplete beta function Iₓ(a, b), for a, b > 0. Pass
  ## `y = 1 − x` when it is known more precisely than `x` (x close to 1).
  if x.isNaN or y.isNaN or not (a > 0.0 and b > 0.0): return NaN
  if x <= 0.0: return 0.0
  if y <= 0.0: return 1.0
  let front = exp(lgamma(a + b) - lgamma(a) - lgamma(b) + a * ln(x) + b * ln(y))
  # The continued fraction converges fastest for x < (a+1)/(a+b+2); otherwise
  # use the symmetry Iₓ(a, b) = 1 − I₁₋ₓ(b, a).
  if x < (a + 1.0) / (a + b + 2.0):
    front * betaContinuedFraction(a, b, x) / a
  else:
    1.0 - front * betaContinuedFraction(b, a, y) / b

func studentTTail(t, df: float): float =
  ## P(T > t) for t ≥ 0, computed directly so small tails keep their precision.
  let t2 = t * t
  0.5 * incompleteBeta(0.5 * df, 0.5, df / (df + t2), t2 / (df + t2))

func studentTPdf(t, df: float): float =
  exp(lgamma(0.5 * (df + 1.0)) - lgamma(0.5 * df) - 0.5 * ln(df * PI) -
      0.5 * (df + 1.0) * ln(1.0 + t * t / df))

func studentTCdf*(t, df: float): float =
  ## P(T ≤ t) for Student's t distribution with `df` degrees of freedom.
  if t.isNaN or not (df > 0.0): return NaN
  if t > 0.0: 1.0 - studentTTail(t, df) else: studentTTail(-t, df)

func studentTQuantile*(p, df: float): float =
  ## The t with P(T ≤ t) = p. NaN when `p` is outside [0, 1] or `df` ≤ 0.
  if p.isNaN or p < 0.0 or p > 1.0 or not (df > 0.0): return NaN
  if p == 0.0: return NegInf
  if p == 1.0: return Inf
  if p == 0.5: return 0.0
  if p < 0.5: return -studentTQuantile(1.0 - p, df)
  # Newton's method on the upper tail, falling back to bisection (or doubling
  # while there is no upper bound) whenever a step leaves the bracket.
  let q = 1.0 - p
  var
    lo = 0.0
    hi = Inf
    t = 1.0
  for _ in 0 ..< 200:
    let excess = studentTTail(t, df) - q
    if excess == 0.0: return t
    if excess > 0.0: lo = t else: hi = t
    var next = t + excess / studentTPdf(t, df)
    if not (next > lo and next < hi):
      next = if hi == Inf: 2.0 * t else: 0.5 * (lo + hi)
    if abs(next - t) <= 1e-15 * next: return next
    t = next
  t

proc unscaledCovariance*[K: static int](
    jtwj: array[K, array[K, float]];
    free: array[K, bool];
    cov: var array[K, array[K, float]]
): bool =
  ## (JᵀWJ)⁻¹ over the free parameters, with zero rows and columns for fixed
  ## ones. Multiply by σ² for the parameter covariance. If the free parameters
  ## are not identifiable (JᵀWJ is singular) their entries are NaN and the
  ## result is false.
  var m = jtwj
  pinFixed(m, free)
  result = invertSpd(m, cov)
  for j in 0 ..< K:
    for k in 0 ..< K:
      if not (free[j] and free[k]): cov[j][k] = 0.0
      elif not result: cov[j][k] = NaN

func gaussianLogLik*(ssr: float; n: int; sumLogWeights = 0.0): float =
  ## Maximised log-likelihood of a least-squares fit with normal errors,
  ## matching R's `logLik.nls`. `n` counts observations with non-zero weight
  ## and `sumLogWeights` is Σ ln(wᵢ) over them.
  let nf = float(n)
  -0.5 * nf * (ln(2.0 * PI) + 1.0 - ln(nf) - sumLogWeights / nf + ln(ssr))
