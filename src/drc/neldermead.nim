## Nelder–Mead downhill simplex minimisation. It needs no derivatives, so it
## is a robust fallback when Levenberg–Marquardt fails to converge.

import std/math

type
  NMOptions* = object
    maxEval*: int = 5000
    ftol*: float = 1e-12  ## stop when the spread of values over the simplex is below ftol · |best|
    xtol*: float = 1e-8   ## ...and each vertex is within xtol (relative) of the best one

  NMResult*[K: static int] = object
    params*: array[K, float]
    value*: float
    evaluations*: int
    converged*: bool

proc nelderMead*[K: static int; F](f: F; p0: array[K, float]; free: array[K, bool];
                                   opts = NMOptions()): NMResult[K] =
  ## Minimises `f(p: array[K, float]): float` over the parameters marked in
  ## `free`, holding the others at their value in `p0`. NaN values count as +Inf.
  var idx: seq[int]
  for j in 0 ..< K:
    if free[j]: idx.add j
  let m = idx.len

  template eval(p: array[K, float]): float =
    inc result.evaluations
    let v = f(p)
    if v.isNaN: Inf else: v

  if m == 0:
    # Nothing to optimise.
    result.params = p0
    result.value = eval(p0)
    result.converged = true
    return

  # Initial simplex: p0 plus a 10% step (or 0.1 from zero) along each free axis.
  var
    simplex = newSeq[array[K, float]](m + 1)
    values = newSeq[float](m + 1)
  for i in 0 .. m:
    simplex[i] = p0
    if i > 0:
      let j = idx[i - 1]
      simplex[i][j] += (if p0[j] != 0.0: 0.1 * abs(p0[j]) else: 0.1)
    values[i] = eval(simplex[i])

  proc along(a, b: array[K, float]; t: float): array[K, float] =
    ## a + t · (b − a) on the free coordinates.
    result = a
    for j in idx: result[j] = a[j] + t * (b[j] - a[j])

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
        if abs(simplex[i][j] - simplex[0][j]) > opts.xtol * (abs(simplex[0][j]) + opts.xtol):
          small = false
    if small:
      result.converged = true
      break
    if result.evaluations >= opts.maxEval:
      break

    var centroid = simplex[0]
    for j in idx:
      var total = 0.0
      for i in 0 ..< m: total += simplex[i][j]
      centroid[j] = total / float(m)

    let
      reflected = along(centroid, simplex[m], -1.0)
      fr = eval(reflected)
    if fr < values[0]:
      let expanded = along(centroid, simplex[m], -2.0)
      let fe = eval(expanded)
      if fe < fr: (simplex[m], values[m]) = (expanded, fe)
      else: (simplex[m], values[m]) = (reflected, fr)
    elif fr < values[m - 1]:
      (simplex[m], values[m]) = (reflected, fr)
    else:
      # Contract towards the better of the reflected and worst points.
      let (target, ft) = if fr < values[m]: (reflected, fr) else: (simplex[m], values[m])
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
