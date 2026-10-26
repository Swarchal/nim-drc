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
  if abs(d) < tiny:
    d = tiny
  d = 1.0 / d
  result = d
  for m in 1 .. maxIter:
    let
      fm = float(m)
      m2 = 2.0 * fm
    # Even step.
    var aa = fm * (b - fm) * x / ((qam + m2) * (a + m2))
    d = 1.0 + aa * d
    if abs(d) < tiny:
      d = tiny
    c = 1.0 + aa / c
    if abs(c) < tiny:
      c = tiny
    d = 1.0 / d
    result *= d * c
    # Odd step.
    aa = -(a + fm) * (qab + fm) * x / ((a + m2) * (qap + m2))
    d = 1.0 + aa * d
    if abs(d) < tiny:
      d = tiny
    c = 1.0 + aa / c
    if abs(c) < tiny:
      c = tiny
    d = 1.0 / d
    let delta = d * c
    result *= delta
    if abs(delta - 1.0) < eps:
      break

func incompleteBeta*(a, b, x: float, y = 1.0 - x): float =
  ## Regularised incomplete beta function Iₓ(a, b), for a, b > 0. Pass
  ## `y = 1 − x` when it is known more precisely than `x` (x close to 1).
  if x.isNaN or y.isNaN or not (a > 0.0 and b > 0.0):
    return NaN
  if x <= 0.0:
    return 0.0
  if y <= 0.0:
    return 1.0
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
  exp(
    lgamma(0.5 * (df + 1.0)) - lgamma(0.5 * df) - 0.5 * ln(df * PI) -
      0.5 * (df + 1.0) * ln(1.0 + t * t / df)
  )

func studentTCdf*(t, df: float): float =
  ## P(T ≤ t) for Student's t distribution with `df` degrees of freedom.
  if t.isNaN or not (df > 0.0):
    return NaN
  if t > 0.0:
    1.0 - studentTTail(t, df)
  else:
    studentTTail(-t, df)

func studentTQuantile*(p, df: float): float =
  ## The t with P(T ≤ t) = p. NaN when `p` is outside [0, 1] or `df` ≤ 0.
  if p.isNaN or p < 0.0 or p > 1.0 or not (df > 0.0):
    return NaN
  if p == 0.0:
    return NegInf
  if p == 1.0:
    return Inf
  if p == 0.5:
    return 0.0
  if p < 0.5:
    return -studentTQuantile(1.0 - p, df)
  # Newton's method on the upper tail, falling back to bisection (or doubling
  # while there is no upper bound) whenever a step leaves the bracket.
  let q = 1.0 - p
  var
    lo = 0.0
    hi = Inf
    t = 1.0
  for _ in 0 ..< 200:
    let excess = studentTTail(t, df) - q
    if excess == 0.0:
      return t
    if excess > 0.0:
      lo = t
    else:
      hi = t
    var next = t + excess / studentTPdf(t, df)
    if not (next > lo and next < hi):
      next =
        if hi == Inf:
          2.0 * t
        else:
          0.5 * (lo + hi)
    if abs(next - t) <= 1e-15 * next:
      return next
    t = next
  t

func normalQuantile*(p: float): float =
  ## The z with P(Z ≤ z) = p for a standard normal Z, by Wichura's algorithm
  ## AS 241 (PPND16), as R's `qnorm`: accurate to about 1e-16. NaN when `p`
  ## is outside [0, 1].
  if p.isNaN or p < 0.0 or p > 1.0:
    return NaN
  if p == 0.0:
    return NegInf
  if p == 1.0:
    return Inf
  let q = p - 0.5
  if abs(q) <= 0.425:
    let r = 0.180625 - q * q
    return
      q * (
        (
          (
            (
              (
                (
                  (r * 2509.0809287301226727 + 33430.575583588128105) * r +
                  67265.770927008700853
                ) * r + 45921.953931549871457
              ) * r + 13731.693765509461125
            ) * r + 1971.5909503065514427
          ) * r + 133.14166789178437745
        ) * r + 3.387132872796366608
      ) / (
        (
          (
            (
              (
                (
                  (r * 5226.495278852545925 + 28729.085735721942674) * r +
                  39307.89580009271061
                ) * r + 21213.794301586595867
              ) * r + 5394.1960214247511077
            ) * r + 687.1870074920579083
          ) * r + 42.313330701600911252
        ) * r + 1.0
      )
  var r = sqrt(
    -ln(
      if q < 0.0:
        p
      else:
        1.0 - p
    )
  )
  var z: float
  if r <= 5.0:
    r -= 1.6
    z =
      (
        (
          (
            (
              (
                (
                  (r * 7.7454501427834140764e-4 + 0.0227238449892691845833) * r +
                  0.24178072517745061177
                ) * r + 1.27045825245236838258
              ) * r + 3.64784832476320460504
            ) * r + 5.7694972214606914055
          ) * r + 4.6303378461565452959
        ) * r + 1.42343711074968357734
      ) / (
        (
          (
            (
              (
                (
                  (r * 1.05075007164441684324e-9 + 5.475938084995344946e-4) * r +
                  0.0151986665636164571966
                ) * r + 0.14810397642748007459
              ) * r + 0.68976733498510000455
            ) * r + 1.6763848301838038494
          ) * r + 2.05319162663775882187
        ) * r + 1.0
      )
  else:
    r -= 5.0
    z =
      (
        (
          (
            (
              (
                (
                  (r * 2.01033439929228813265e-7 + 2.71155556874348757815e-5) * r +
                  0.0012426609473880784386
                ) * r + 0.026532189526576123093
              ) * r + 0.29656057182850489123
            ) * r + 1.7848265399172913358
          ) * r + 5.4637849111641143699
        ) * r + 6.6579046435011037772
      ) / (
        (
          (
            (
              (
                (
                  (r * 2.04426310338993978564e-15 + 1.4215117583164458887e-7) * r +
                  1.8463183175100546818e-5
                ) * r + 7.868691311456132591e-4
              ) * r + 0.0148753612908506148525
            ) * r + 0.13692988092273580531
          ) * r + 0.59983220655588793769
        ) * r + 1.0
      )
  if q < 0.0:
    -z
  else:
    z

proc unscaledCovariance*[M, B](jtwj: M, free: B, cov: var M): bool =
  ## (JᵀWJ)⁻¹ over the free parameters, with zero rows and columns for fixed
  ## ones. Multiply by σ² for the parameter covariance. If the free parameters
  ## are not identifiable (JᵀWJ is singular) their entries are NaN and the
  ## result is false. Works on fixed-size arrays or seqs, as `linalg`.
  let n = jtwj.len
  var m = jtwj
  pinFixed(m, free)
  cov.ensureShape n
  result = invertSpd(m, cov)
  for j in 0 ..< n:
    for k in 0 ..< n:
      if not (free[j] and free[k]):
        cov[j][k] = 0.0
      elif not result:
        cov[j][k] = NaN

func fTail*(f, df1, df2: float): float =
  ## P(F > f) for the F distribution with `df1` and `df2` degrees of freedom.
  if f.isNaN or not (df1 > 0.0 and df2 > 0.0):
    return NaN
  if f <= 0.0:
    return 1.0
  if f == Inf:
    return 0.0
  let d = df2 + df1 * f
  incompleteBeta(0.5 * df2, 0.5 * df1, df2 / d, df1 * f / d)

func gaussianLogLik*(ssr: float, n: int, sumLogWeights = 0.0): float =
  ## Maximised log-likelihood of a least-squares fit with normal errors,
  ## matching R's `logLik.nls`. `n` counts observations with non-zero weight
  ## and `sumLogWeights` is Σ ln(wᵢ) over them.
  let nf = float(n)
  -0.5 * nf * (ln(2.0 * PI) + 1.0 - ln(nf) - sumLogWeights / nf + ln(ssr))
