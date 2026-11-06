## Nelder–Mead downhill simplex minimisation. It needs no derivatives, so it
## is a robust fallback when Levenberg–Marquardt fails to converge.

import std/math

type
  NMOptions* = object
    maxEval*: int = 5000
    ftol*: float = 1e-12
      ## stop when the spread of values over the simplex is below ftol · |best|
    xtol*: float = 1e-8 ## ...and each vertex is within xtol (relative) of the best one

  NMResult*[K: static int] = object
    params*: array[K, float]
    value*: float
    evaluations*: int
    converged*: bool

  NMDynResult* = object ## `NMResult` for a parameter count known only at run time
    params*: seq[float]
    value*: float
    evaluations*: int
    converged*: bool

proc nmCore[V, B, R; F](f: F, p0: V, free: B, lower, upper: V, opts: NMOptions): R =
  ## The body of both `nelderMead` overloads: `V` is an `array[K, float]` or a
  ## `seq[float]`, and `R` the matching result type.
  var idx: seq[int]
  for j in 0 ..< p0.len:
    if free[j]:
      idx.add j
  let m = idx.len

  proc project(p: V): V =
    result = p
    for j in idx:
      result[j] = clamp(p[j], lower[j], upper[j])

  template eval(p: V): float =
    inc result.evaluations
    let v = f(p)
    if v.isNaN: Inf else: v

  let start = project(p0)
  if m == 0:
    # Nothing to optimise.
    result.params = start
    result.value = eval(start)
    result.converged = true
    return

  # Initial simplex: the start plus a 10% step (or 0.1 from zero) along each
  # free axis, stepping down instead where that would cross the upper bound.
  var
    simplex = newSeq[V](m + 1)
    values = newSeq[float](m + 1)
  for i in 0 .. m:
    simplex[i] = start
    if i > 0:
      let j = idx[i - 1]
      let step =
        if start[j] != 0.0:
          0.1 * abs(start[j])
        else:
          0.1
      simplex[i][j] =
        if start[j] + step <= upper[j]:
          start[j] + step
        else:
          max(start[j] - step, lower[j])
    values[i] = eval(simplex[i])

  proc along(a, b: V, t: float): V =
    ## a + t · (b − a) on the free coordinates, clamped into the bounds.
    result = a
    for j in idx:
      result[j] = clamp(a[j] + t * (b[j] - a[j]), lower[j], upper[j])

  while true:
    # Order the vertices: best first, worst last.
    for i in 1 .. m:
      var k = i
      while k > 0 and values[k] < values[k - 1]:
        swap(values[k], values[k - 1])
        swap(simplex[k], simplex[k - 1])
        dec k

    var small = values[m] - values[0] <= opts.ftol * abs(values[0])
    for i in 1 .. m:
      for j in idx:
        if abs(simplex[i][j] - simplex[0][j]) >
            opts.xtol * (abs(simplex[0][j]) + opts.xtol):
          small = false
    if small:
      result.converged = true
      break
    if result.evaluations >= opts.maxEval:
      break

    var centroid = simplex[0]
    for j in idx:
      var total = 0.0
      for i in 0 ..< m:
        total += simplex[i][j]
      centroid[j] = total / float(m)

    let
      reflected = along(centroid, simplex[m], -1.0)
      fr = eval(reflected)
    if fr < values[0]:
      let expanded = along(centroid, simplex[m], -2.0)
      let fe = eval(expanded)
      if fe < fr:
        (simplex[m], values[m]) = (expanded, fe)
      else:
        (simplex[m], values[m]) = (reflected, fr)
    elif fr < values[m - 1]:
      (simplex[m], values[m]) = (reflected, fr)
    else:
      # Contract towards the better of the reflected and worst points.
      let (target, ft) =
        if fr < values[m]:
          (reflected, fr)
        else:
          (simplex[m], values[m])
      let contracted = along(centroid, target, 0.5)
      let fc = eval(contracted)
      if fc < ft:
        (simplex[m], values[m]) = (contracted, fc)
      else:
        # Shrink everything towards the best vertex.
        for i in 1 .. m:
          simplex[i] = along(simplex[0], simplex[i], 0.5)
          values[i] = eval(simplex[i])

  result.params = simplex[0]
  result.value = values[0]

proc nelderMead*[K: static int, F](
    f: F,
    p0: array[K, float],
    free: array[K, bool],
    lower, upper: array[K, float],
    opts = NMOptions(),
): NMResult[K] =
  ## Minimises `f(p: array[K, float]): float` over the parameters marked in
  ## `free`, within `lower` ≤ p ≤ `upper` (±Inf for no bound), holding the
  ## others at their value in `p0`. NaN values count as +Inf. Every vertex is
  ## clamped into the bounds, as SciPy's bounded Nelder–Mead.
  nmCore[array[K, float], array[K, bool], NMResult[K], F](
    f, p0, free, lower, upper, opts
  )

proc nelderMead*[F](
    f: F, p0: seq[float], free: seq[bool], lower, upper: seq[float], opts = NMOptions()
): NMDynResult =
  ## As the array overload, for a parameter count known only at run time
  ## (`f(p: seq[float]): float`).
  nmCore[seq[float], seq[bool], NMDynResult, F](f, p0, free, lower, upper, opts)

proc nelderMead*[K: static int, F](
    f: F, p0: array[K, float], free: array[K, bool], opts = NMOptions()
): NMResult[K] =
  ## Without bounds.
  var lower, upper: array[K, float]
  for j in 0 ..< K:
    (lower[j], upper[j]) = (-Inf, Inf)
  nelderMead(f, p0, free, lower, upper, opts)
