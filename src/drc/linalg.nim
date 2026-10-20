## Small fixed-size dense linear algebra, enough for the normal equations
## of a nonlinear least-squares fit with a handful of parameters.

import std/math

func isFinite*(x: float): bool = x.classify notin {fcInf, fcNegInf, fcNan}

proc cholesky*[N: static int](a: array[N, array[N, float]]; l: var array[N, array[N, float]]): bool =
  ## Computes lower-triangular `l` with `a = l * lᵀ`. Only the lower triangle
  ## of `a` is read. Returns false if `a` is not (numerically) positive definite.
  l = default(array[N, array[N, float]])
  for j in 0 ..< N:
    var d = a[j][j]
    for k in 0 ..< j:
      d -= l[j][k] * l[j][k]
    if not (d > 0.0):  # also rejects NaN
      return false
    l[j][j] = sqrt(d)
    for i in j + 1 ..< N:
      var s = a[i][j]
      for k in 0 ..< j:
        s -= l[i][k] * l[j][k]
      l[i][j] = s / l[j][j]
  true

proc choleskySolve*[N: static int](l: array[N, array[N, float]]; b: array[N, float]): array[N, float] =
  ## Solves `l * lᵀ * x = b` given the Cholesky factor `l`.
  var z: array[N, float]
  for i in 0 ..< N:
    var s = b[i]
    for k in 0 ..< i:
      s -= l[i][k] * z[k]
    z[i] = s / l[i][i]
  for i in countdown(N - 1, 0):
    var s = z[i]
    for k in i + 1 ..< N:
      s -= l[k][i] * result[k]
    result[i] = s / l[i][i]

proc solveSpd*[N: static int](a: array[N, array[N, float]]; b: array[N, float]; x: var array[N, float]): bool =
  ## Solves `a * x = b` for symmetric positive definite `a`.
  var l: array[N, array[N, float]]
  if not cholesky(a, l):
    return false
  x = choleskySolve(l, b)
  true

proc invertSpd*[N: static int](a: array[N, array[N, float]]; inv: var array[N, array[N, float]]): bool =
  ## Inverts symmetric positive definite `a`. Returns false if it is not positive definite.
  var l: array[N, array[N, float]]
  if not cholesky(a, l):
    return false
  for j in 0 ..< N:
    var e: array[N, float]
    e[j] = 1.0
    let col = choleskySolve(l, e)
    # Fill the lower triangle from this column and mirror it, so the inverse
    # is exactly symmetric.
    for i in j ..< N:
      inv[i][j] = col[i]
      inv[j][i] = col[i]
  true

proc pinFixed*[N: static int](m: var array[N, array[N, float]]; free: array[N, bool]) =
  ## Replaces the rows and columns of fixed (non-free) parameters with those of
  ## the identity, so they decouple from the free ones in a linear solve.
  for j in 0 ..< N:
    if not free[j]:
      for k in 0 ..< N:
        m[j][k] = 0.0
        m[k][j] = 0.0
      m[j][j] = 1.0
