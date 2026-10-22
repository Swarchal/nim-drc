import std/[math, strutils, unittest]
import drc
import ryegrass

proc close(a, b: float, rel: float): bool =
  abs(a - b) <= rel * max(1e-3, abs(b))

proc near(a, b, se: float, k = 2e-5): bool =
  ## Estimates agree to a small fraction of their standard error (see tmulti).
  abs(a - b) <= k * se

proc expDecay(
    x: float, p: array[2, float], grad: var array[2, float]
): float {.nimcall.} =
  let e = exp(-p[1] * x)
  grad = [e, -p[0] * x * e]
  p[0] * e

suite "bounded optimisers":
  # Noiseless exponential decay with A = 3, k = 0.7.
  var x, y, w: seq[float]
  var grad: array[2, float]
  for k in 0 ..< 20:
    x.add 0.25 * float(k)
    y.add expDecay(x[^1], [3.0, 0.7], grad)
    w.add 1.0
  # With k held at 0.5, the best A is linear least squares on e^(−0.5x).
  var sxy, sxx = 0.0
  for i in 0 ..< x.len:
    sxy += y[i] * exp(-0.5 * x[i])
    sxx += exp(-x[i])
  let bestA = sxy / sxx

  test "Levenberg–Marquardt stops on an active bound":
    let r = levenbergMarquardt(
      expDecay, x, y, w, [1.0, 0.1], [true, true], [-Inf, -Inf], [Inf, 0.5]
    )
    check r.status.converged
    check r.params[1] == 0.5
    check close(r.params[0], bestA, 1e-9)

  test "a start outside the bounds is clamped into them":
    let r = levenbergMarquardt(
      expDecay, x, y, w, [1.0, 2.0], [true, true], [-Inf, -Inf], [Inf, 0.5]
    )
    check r.status.converged
    check r.params[1] == 0.5
    check close(r.params[0], bestA, 1e-9)

  test "an inactive bound changes nothing":
    let r = levenbergMarquardt(
      expDecay, x, y, w, [1.0, 0.1], [true, true], [0.0, 0.0], [10.0, 10.0]
    )
    check r.status.converged
    check close(r.params[0], 3.0, 1e-8)
    check close(r.params[1], 0.7, 1e-8)

  test "the runtime-sized overload takes bounds":
    proc decaySeq(x: float, p: seq[float], grad: var seq[float]): float =
      let e = exp(-p[1] * x)
      grad[0] = e
      grad[1] = -p[0] * x * e
      p[0] * e

    let r = levenbergMarquardt(
      decaySeq, x, y, w, @[1.0, 0.1], @[true, true], @[-Inf, -Inf], @[Inf, 0.5]
    )
    check r.status.converged
    check r.params[1] == 0.5
    check close(r.params[0], bestA, 1e-9)

  test "Nelder–Mead stays in the box":
    let r = nelderMead(
      proc(p: array[2, float]): float =
        (p[0] - 3.0) ^ 2 + (p[1] + 1.0) ^ 2,
      [0.0, 2.0],
      [true, true],
      [-Inf, 0.0],
      [2.0, Inf],
    )
    check r.converged
    check abs(r.params[0] - 2.0) < 1e-6
    check abs(r.params[1]) < 1e-6
    check r.params[0] <= 2.0 and r.params[1] >= 0.0

suite "fitLL4 with bounds: ryegrass":
  # References from R: nls(y ~ b + (t-b)/(1+10^((l-x)*h)), algorithm = "port",
  # lower = ..., upper = ...). Its standard errors treat parameters at a
  # bound as free, as SciPy's curve_fit does (which agrees to 7 digits).
  test "Hill and Top both at a bound":
    let fit = fitLL4(
      ryegrassDose,
      ryegrassRootLength,
      bounds = [(Hill, -2.5 .. Inf), (Top, -Inf .. 7.5)],
    )
    check fit.converged
    check fit.params[Hill] == -2.5
    check fit.params[Top] == 7.5
    check fit.atBound == {Top, Hill}
    check fwAtBound in fit.warnings
    check fit.df == 20
    check near(fit.params[Bottom], 0.317975155252540, 0.2709739602447437)
    check near(fit.params[LogEC50], 0.515137630966396, 0.0351880061587067)
    check close(fit.se[Bottom], 0.2709739602447437, 1e-5)
    check close(fit.se[Top], 0.2210532003248256, 1e-5)
    check close(fit.se[LogEC50], 0.0351880061587067, 1e-5)
    check close(fit.se[Hill], 0.4436246902454261, 1e-5)
    check close(fit.ssr, 7.01168204339883, 1e-9)
    check close(fit.sigma, 0.592101428954484, 1e-9)
    check close(fit.aic, 48.5776205789046, 1e-9)

  test "Hill at a bound, the rest free":
    let fit = fitLL4(ryegrassDose, ryegrassRootLength, bounds = [(Hill, -2.5 .. Inf)])
    check fit.converged
    check fit.atBound == {Hill}
    check fit.params[Hill] == -2.5
    check near(fit.params[Bottom], 0.351300738754152, 0.2390451584497744)
    check near(fit.params[Top], 7.877861671538866, 0.2022667256380237)
    check near(fit.params[LogEC50], 0.488961424703306, 0.0303990282365878)
    check close(fit.se[Hill], 0.3768188403397014, 1e-5)
    check close(fit.ssr, 5.76261374862384, 1e-9)

  test "a lower bound on Bottom":
    let fit = fitLL4(ryegrassDose, ryegrassRootLength, bounds = [(Bottom, 0.6 .. Inf)])
    check fit.converged
    check fit.atBound == {Bottom}
    check fit.params[Bottom] == 0.6
    check near(fit.params[Top], 7.781382090537693, 0.1902619724458674)
    check near(fit.params[LogEC50], 0.479186889639322, 0.0264241508264739)
    check near(fit.params[Hill], -3.101131048440974, 0.4875977669438525)
    check close(fit.se[Bottom], 0.2063521070433216, 1e-5)
    check close(fit.ssr, 5.48669241681981, 1e-9)

  test "bounds the fit does not reach give the unbounded fit":
    let free = fitLL4(ryegrassDose, ryegrassRootLength)
    let bounded = fitLL4(
      ryegrassDose,
      ryegrassRootLength,
      bounds = [(Bottom, 0.0 .. 1.0), (Top, 5.0 .. 10.0), (Hill, -10.0 .. -0.1)],
    )
    check bounded.atBound == {}
    check fwAtBound notin bounded.warnings
    for param in [Bottom, Top, LogEC50, Hill]:
      check near(bounded.params[param], free.params[param], free.se[param], 1e-6)
      check close(bounded.se[param], free.se[param], 1e-6)

  test "the summary marks parameters at a bound":
    let fit = fitLL4(
      ryegrassDose,
      ryegrassRootLength,
      bounds = [(Hill, -2.5 .. Inf), (Top, -Inf .. 7.5)],
    )
    let s = fit.summary
    check "(at upper bound)" in s
    check "(at lower bound)" in s
    check "Warning: a parameter is at a bound" in s

  test "the Nelder–Mead fallback respects the bounds":
    let fit = fitLL4(
      ryegrassDose,
      ryegrassRootLength,
      bounds = [(Hill, -2.5 .. Inf), (Top, -Inf .. 7.5)],
      opts = FitOptions(lm: LMOptions(maxIter: 2), multiStart: false),
    )
    check fit.usedFallback
    check fit.params[Hill] >= -2.5
    check fit.params[Top] <= 7.5

  test "3PL fits take bounds":
    # Unbounded, the 3PL's Hill is −2.47.
    let fit = fitLL3(ryegrassDose, ryegrassRootLength, bounds = [(Hill, -2.0 .. Inf)])
    check fit.converged
    check fit.params[Bottom] == 0.0
    check fit.params[Hill] == -2.0
    check fit.atBound == {Hill}
    check near(fit.params[Top], 7.973695046936829, 0.2285289575699544)
    check near(fit.params[LogEC50], 0.510169303319511, 0.0322015004574352)
    check close(fit.se[Hill], 0.2444988142639670, 1e-5)
    check close(fit.ssr, 7.45022737860643, 1e-9)

  test "invalid bounds raise":
    expect ValueError:
      discard fitLL4(ryegrassDose, ryegrassRootLength, bounds = [(Hill, 1.0 .. 1.0)])
    expect ValueError:
      discard fitLL4(ryegrassDose, ryegrassRootLength, bounds = [(Hill, 2.0 .. 1.0)])
    expect ValueError:
      discard fitLL4(ryegrassDose, ryegrassRootLength, bounds = [(Hill, NaN .. 1.0)])
    expect ValueError: # fixed and bounded
      discard fitLL4(
        ryegrassDose,
        ryegrassRootLength,
        fixed = [(Top, 8.0)],
        bounds = [(Top, 0.0 .. 9.0)],
      )
    expect ValueError: # Bottom is fixed in a 3PL
      discard fitLL3(ryegrassDose, ryegrassRootLength, bounds = [(Bottom, 0.0 .. 1.0)])
    expect ValueError:
      discard
        fitLL4(ryegrassDose, ryegrassRootLength, bounds = [(Asymmetry, 0.5 .. 2.0)])

suite "fitLL5 with bounds":
  test "Asymmetry at its lower bound matches nls port":
    # The README's asymmetric data; unbounded, Asymmetry is 0.392.
    var doses = @[0.0, 0.0]
    for k in 0 .. 12:
      let dose = pow(10.0, -9.0 + 0.5 * float(k))
      doses.add [dose, dose]
    let responses = [
      5.0, 5.6, 8.34, 7.1, 10.65, 9.57, 16.18, 18.74, 22.64, 22.38, 37.01, 36.75, 54.87,
      52.8, 75.55, 77.0, 86.17, 87.94, 89.67, 90.89, 90.97, 94.18, 92.39, 95.46, 95.3,
      94.61, 89.96, 93.92,
    ]
    let fit = fitLL5(doses, responses, bounds = [(Asymmetry, 0.5 .. Inf)])
    check fit.converged
    check fit.atBound == {Asymmetry}
    check fit.params[Asymmetry] == 0.5
    check near(fit.params[Bottom], 5.53390919615597, 1.035424243281960)
    check near(fit.params[Top], 93.54482040029316, 0.622968385782773)
    check near(fit.params[LogEC50], -6.19947116181473, 0.114909017774703)
    check near(fit.params[Hill], 1.06494747050679, 0.118797307766521)
    check close(fit.se[Asymmetry], 0.104195765402131, 1e-5)
    check close(fit.se[Hill], 0.118797307766521, 1e-5)
    check close(fit.ssr, 65.050503369584, 1e-9)

  test "bounds that exclude Asymmetry = 1 hold the symmetric stage inside them":
    var doses = @[0.0, 0.0]
    for k in 0 .. 12:
      let dose = pow(10.0, -9.0 + 0.5 * float(k))
      doses.add [dose, dose]
    let responses = [
      5.0, 5.6, 8.34, 7.1, 10.65, 9.57, 16.18, 18.74, 22.64, 22.38, 37.01, 36.75, 54.87,
      52.8, 75.55, 77.0, 86.17, 87.94, 89.67, 90.89, 90.97, 94.18, 92.39, 95.46, 95.3,
      94.61, 89.96, 93.92,
    ]
    let free = fitLL5(doses, responses)
    let fit = fitLL5(doses, responses, bounds = [(Asymmetry, 0.1 .. 0.8)])
    check fit.atBound == {}
    check near(fit.params[Asymmetry], free.params[Asymmetry], free.se[Asymmetry], 1e-5)

suite "robust fits with bounds":
  test "Huber with Hill at a bound matches nlrob port":
    # robustbase::nlrob(..., psi = .Mwgt.psi1("huber", 1.345), algorithm = "port",
    # lower = c(h = -1)), tol 1e-8; port's own precision limits the reference to
    # about 3e-7 (relative).
    const
      dose = [
        1e-09, 1e-09, 1e-09, 3.1622776601683795e-09, 3.1622776601683795e-09,
        3.1622776601683795e-09, 1e-08, 1e-08, 1e-08, 3.162277660168379e-08,
        3.162277660168379e-08, 3.162277660168379e-08, 1e-07, 1e-07, 1e-07,
        3.162277660168379e-07, 3.162277660168379e-07, 3.162277660168379e-07, 1e-06,
        1e-06, 1e-06, 3.162277660168379e-06, 3.162277660168379e-06,
        3.162277660168379e-06, 1e-05, 1e-05, 1e-05, 3.1622776601683795e-05,
        3.1622776601683795e-05, 3.1622776601683795e-05, 0.0001, 0.0001, 0.0001,
      ]
      response = [
        104.03, 90.24, 99.16, 95.92, 136.26, 96.97, 90.44, 95.81, 93.91, 102.27, 92.98,
        91.24, 77.88, 76.72, 75.56, 48.83, 16.45, 49.28, 24.15, 20.67, 21.35, 12.34,
        9.33, 6.18, 2.95, 35.12, 9.3, 1.57, 1.65, 5.39, -0.56, 1.22, 4.74,
      ]
    let fit = fitLL4(
      dose, response, bounds = [(Hill, -1.0 .. Inf)], opts = FitOptions(loss: Huber)
    )
    check fit.converged
    check fit.atBound == {Hill}
    check fit.params[Hill] == -1.0
    check close(fit.params[Bottom], 2.03171458855001, 1e-6)
    check close(fit.params[Top], 100.083282865875, 1e-6)
    check close(fit.params[LogEC50], -6.5467662623357, 1e-6)
    check close(fit.sigma, 3.60364962020509, 1e-6)
    check close(fit.robustWeights[4], 0.130095033154273, 1e-5)
    check close(fit.robustWeights[16], 0.151606080488306, 1e-5)

suite "joint fits with bounds":
  test "parallel DNase curves with Top at a bound match nls port":
    # nls(density ~ b + (t-b)/(1+10^((l[c]-x)*h)), algorithm = "port",
    # upper = c(Inf, 2, Inf, Inf, Inf, Inf)); unbounded, Top is 2.296.
    const
      conc = [0.04882812, 0.1953125, 0.390625, 0.78125, 1.5625, 3.125, 6.25, 12.5]
      density = [
        0.017, 0.018, 0.121, 0.124, 0.206, 0.215, 0.377, 0.374, 0.614, 0.609, 1.019,
        1.001, 1.334, 1.364, 1.73, 1.71, 0.011, 0.016, 0.118, 0.108, 0.2, 0.206, 0.364,
        0.36, 0.62, 0.64, 0.979, 0.973, 1.424, 1.399, 1.74, 1.732, 0.054, 0.054, 0.152,
        0.148, 0.226, 0.222, 0.392, 0.383, 0.658, 0.644, 1.043, 1.002, 1.466, 1.381,
        1.743, 1.724,
      ]
    var doses: seq[float]
    var curve: seq[int]
    for c in 0 .. 2:
      for d in conc:
        doses.add [d, d]
        curve.add [c, c]
    let fit = fitLL4Multi(
      doses, density, curve, shared = {Bottom, Top, Hill}, bounds = [(Top, -Inf .. 2.0)]
    )
    check fit.converged
    check fit.atBound == {Top}
    check fwAtBound in fit.warnings
    check fit.params[fit.index[0][Top]] == 2.0
    const
      expected = [
        0.0482560924153271, 2.0, 0.5090639705275717, 0.5011658447255355,
        0.4776673662131824, 1.2053495896159236,
      ]
      se = [
        0.0119088734986862, 0.0458727820739963, 0.0212789838371984, 0.0213136532771704,
        0.0214177940751630, 0.0484449304973926,
      ]
    for i in 0 ..< expected.len:
      check near(fit.params[i], expected[i], se[i])
      check close(sqrt(fit.cov[i][i]), se[i], 1e-5)
    check close(fit.ssr, 0.0400669504426697, 1e-9)
    check fit.df == 42
    for c in 0 .. 2:
      check fit.curves[c].atBound == {Top}
      check fwAtBound in fit.curves[c].warnings
    check "(at upper bound)" in fit.summary
