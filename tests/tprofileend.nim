import std/[math, unittest]
include drc # white-box: profileEnd is private

# Toy models in which parameter 1 does not affect the fit near its estimate,
# so its profile is flat (τ = 0) there and the data do not bound it.

const
  x = [1.0, 2.0, 3.0, 4.0]
  y = [1.0, 2.0, 1.0, 2.0]
  w = [1.0, 1.0, 1.0, 1.0]
  free = [true, true, false, false, false]
  pHat = [1.5, 1.0, 0.0, 0.0, 0.0]
  ssrHat = 1.0 # at p[0] = mean(y)
  tq = 2.0

proc filled(value: float): ParamVec =
  for j in 0 ..< drcParamCount:
    result[j] = value

proc nanAtZero(x: float, p: ParamVec, grad: var ParamVec): float {.nimcall.} =
  ## Undefined for p[1] ≤ 0, as `ll5` is for Asymmetry ≤ 0.
  grad = default(ParamVec)
  grad[0] = 1.0
  p[0] + 0.0 * ln(p[1])

proc risesFarOut(x: float, p: ParamVec, grad: var ParamVec): float {.nimcall.} =
  ## Flat in p[1] down to −100, then a slope k·x with k = 1e-3·(−p[1] − 100)²
  ## that p[0] cannot absorb: like refits far out that stop short of the
  ## minimum and look like a crossing.
  let excess = max(0.0, -p[1] - 100.0)
  grad = default(ParamVec)
  grad[0] = 1.0
  grad[1] = -2e-3 * excess * x
  p[0] + 1e-3 * excess * excess * x

suite "profile ends beyond the reach":
  let upper = filled(Inf)

  test "a positive parameter's end is 0, never evaluated there":
    # The reach (10·tq·step = 0.02) is passed long before 0.
    let lower = filled(-Inf)
    check profileEnd(
      nanAtZero,
      x,
      y,
      w,
      pHat,
      1,
      free,
      lower,
      upper,
      ssrHat,
      1.0,
      tq,
      1e-3,
      -1.0,
      true,
      LMOptions(),
    ) == 0.0
    # Within the reach, 0 is approached by halving.
    check profileEnd(
      nanAtZero,
      x,
      y,
      w,
      pHat,
      1,
      free,
      lower,
      upper,
      ssrHat,
      1.0,
      tq,
      0.1,
      -1.0,
      true,
      LMOptions(),
    ) == 0.0

  test "a loose bound beyond the reach is the end, without a refit there":
    var lower = filled(-Inf)
    lower[1] = -1e6
    check profileEnd(
      risesFarOut,
      x,
      y,
      w,
      pHat,
      1,
      free,
      lower,
      upper,
      ssrHat,
      1.0,
      tq,
      1e-3,
      -1.0,
      false,
      LMOptions(),
    ) == -1e6
    check profileEnd(
      risesFarOut,
      x,
      y,
      w,
      pHat,
      1,
      free,
      filled(-Inf),
      upper,
      ssrHat,
      1.0,
      tq,
      1e-3,
      -1.0,
      false,
      LMOptions(),
    ) == -Inf

  test "within the reach, a crossing is found as before":
    # With a large step the reach (2000) covers the rise past −100. With the
    # best p[0], SSR(k) = 1 − 2k + 5k², which reaches Ŝ + tq²·σ² = 5 at
    # k = (2 + √84)/10, i.e. β = −100 − √(1000·k).
    var lower = filled(-Inf)
    lower[1] = -1e3
    let e = profileEnd(
      risesFarOut,
      x,
      y,
      w,
      pHat,
      1,
      free,
      lower,
      upper,
      ssrHat,
      1.0,
      tq,
      100.0,
      -1.0,
      false,
      LMOptions(),
    )
    check abs(e - (-100.0 - sqrt(1000.0 * (2.0 + sqrt(84.0)) / 10.0))) < 1e-6
