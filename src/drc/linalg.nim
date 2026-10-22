## Small dense linear algebra, enough for the normal equations of a nonlinear
## least-squares fit. Every routine takes either fixed-size arrays
## (`array[N, float]`, `array[N, array[N, float]]`) or their runtime-sized
## `seq` equivalents, so single-curve fits keep their stack-allocated arrays
## while joint fits size the problem at run time.

import std/math

func isFinite*(x: float): bool =
  x.classify notin {fcInf, fcNegInf, fcNan}

# Sizing helpers: no-ops for arrays, whose size is part of the type.

proc ensureLen*[N: static int](v: var array[N, float], n: int) {.inline.} =
  discard

proc ensureLen*(v: var seq[float], n: int) {.inline.} =
  v.setLen n

proc ensureShape*[N: static int](m: var array[N, array[N, float]], n: int) {.inline.} =
  discard

proc ensureShape*(m: var seq[seq[float]], n: int) =
  ## Makes `m` n×n (new entries are 0).
  m.setLen n
  for row in m.mitems:
    row.setLen n

proc setZero*[N: static int](v: var array[N, float]) {.inline.} =
  v = default(array[N, float])

proc setZero*(v: var seq[float]) {.inline.} =
  for x in v.mitems:
    x = 0.0

proc setZero*[N: static int](m: var array[N, array[N, float]]) {.inline.} =
  m = default(array[N, array[N, float]])

proc setZero*(m: var seq[seq[float]]) =
  for row in m.mitems:
    row.setZero()

proc cholesky*[M](a: M, l: var M): bool =
  ## Computes lower-triangular `l` with `a = l * lᵀ`. Only the lower triangle
  ## of `a` is read. Returns false if `a` is not (numerically) positive definite.
  let n = a.len
  l.ensureShape n
  l.setZero()
  for j in 0 ..< n:
    var d = a[j][j]
    for k in 0 ..< j:
      d -= l[j][k] * l[j][k]
    if not (d > 0.0): # also rejects NaN
      return false
    l[j][j] = sqrt(d)
    for i in j + 1 ..< n:
      var s = a[i][j]
      for k in 0 ..< j:
        s -= l[i][k] * l[j][k]
      l[i][j] = s / l[j][j]
  true

proc choleskySolve*[M, V](l: M, b: V): V =
  ## Solves `l * lᵀ * x = b` given the Cholesky factor `l`.
  let n = l.len
  var z: V
  z.ensureLen n
  result.ensureLen n
  for i in 0 ..< n:
    var s = b[i]
    for k in 0 ..< i:
      s -= l[i][k] * z[k]
    z[i] = s / l[i][i]
  for i in countdown(n - 1, 0):
    var s = z[i]
    for k in i + 1 ..< n:
      s -= l[k][i] * result[k]
    result[i] = s / l[i][i]

proc solveSpd*[M, V](a: M, b: V, x: var V): bool =
  ## Solves `a * x = b` for symmetric positive definite `a`.
  var l: M
  if not cholesky(a, l):
    return false
  x = choleskySolve(l, b)
  true

proc invertSpd*[M](a: M, inv: var M): bool =
  ## Inverts symmetric positive definite `a`. Returns false if it is not positive definite.
  let n = a.len
  var l: M
  if not cholesky(a, l):
    return false
  inv.ensureShape n
  for j in 0 ..< n:
    var e: typeof(a[0])
    e.ensureLen n
    e[j] = 1.0
    let col = choleskySolve(l, e)
    # Fill the lower triangle from this column and mirror it, so the inverse
    # is exactly symmetric.
    for i in j ..< n:
      inv[i][j] = col[i]
      inv[j][i] = col[i]
  true

proc pinFixed*[M, B](m: var M, free: B) =
  ## Replaces the rows and columns of fixed (non-free) parameters with those of
  ## the identity, so they decouple from the free ones in a linear solve.
  let n = m.len
  for j in 0 ..< n:
    if not free[j]:
      for k in 0 ..< n:
        m[j][k] = 0.0
        m[k][j] = 0.0
      m[j][j] = 1.0
