import std/[math, random, unittest]
import drc
import ryegrass

proc close(a, b: float; rel = 1e-6): bool = abs(a - b) <= rel * max(1.0, abs(b))

proc expDecay(x: float; p: array[2, float]; grad: var array[2, float]): float {.nimcall.} =
  let e = exp(-p[1] * x)
  grad = [e, -p[0] * x * e]
  p[0] * e

suite "Levenberg-Marquardt":
  # Noiseless exponential decay with A = 3, k = 0.7.
  var x, y, w: seq[float]
  var grad: array[2, float]
  for k in 0 ..< 20:
    x.add 0.25 * float(k)
    y.add expDecay(x[^1], [3.0, 0.7], grad)
    w.add 1.0

  test "fits a generic two-parameter model":
    let r = levenbergMarquardt(expDecay, x, y, w, [1.0, 0.1], [true, true])
    check r.status.converged
    check close(r.params[0], 3.0, 1e-8)
    check close(r.params[1], 0.7, 1e-8)

  test "fixed parameters do not move":
    let r = levenbergMarquardt(expDecay, x, y, w, [2.0, 0.1], [false, true])
    check r.params[0] == 2.0
    check r.status.converged

  test "a stall away from a minimum is not reported as converged":
    # A finite model value with a NaN gradient makes every damped solve fail.
    proc broken(x: float; p: array[2, float]; grad: var array[2, float]): float =
      grad = [NaN, NaN]
      p[0] + p[1] * x
    let r = levenbergMarquardt(broken, x, y, w, [1.0, 0.1], [true, true])
    check r.status == lmStalled
    check not r.status.converged

  test "normalEquations returns JᵀWJ for all parameters":
    var a: array[2, array[2, float]]
    var g: array[2, float]
    discard normalEquations(expDecay, x, y, w, [3.0, 0.7], a, g)
    var expected = 0.0
    for xi in x: expected += exp(-1.4 * xi)  # (∂f/∂A)² = e^(−2kx)
    check abs(a[0][0] - expected) <= 1e-12 * expected
    check a[0][1] == a[1][0]

  test "non-finite start is reported":
    let r = levenbergMarquardt(expDecay, [1.0], [NaN], [1.0], [1.0, 1.0], [true, true])
    check r.status == lmBadStart
    check not r.status.converged

suite "fitLL4: ryegrass reference":
  # Reference values agree between scipy.optimize.curve_fit, R nls() and the
  # published drc::drm(rootl ~ conc, fct = LL.4()) estimates
  # (b = 2.98222, c = 0.48141, d = 7.79296, e = 3.05795; drc's b = -Hill).
  test "all observations, including zero-dose controls":
    let fit = fitLL4(ryegrassDose, ryegrassRootLength)
    check fit.converged
    check close(fit.params[Bottom], 0.48140994064097836)
    check close(fit.params[Top], 7.792962329581427)
    check close(fit.params[LogEC50], 0.48543111513116)
    check close(fit.params[Hill], -2.9822286375650706)
    check close(fit.ec50, 3.0579551757658026)
    check close(fit.ssr, 5.400214650161578, 1e-10)

  test "zero doses dropped":
    var dose, resp: seq[float]
    for i in 0 ..< ryegrassDose.len:
      if ryegrassDose[i] > 0.0:
        dose.add ryegrassDose[i]
        resp.add ryegrassRootLength[i]
    let fit = fitLL4(dose, resp)
    check fit.converged
    check close(fit.params[Bottom], 0.4560145651203796)
    check close(fit.params[Top], 7.975508629556834)
    check close(fit.params[LogEC50], 0.4760056732498162)
    check close(fit.params[Hill], -2.809326109902281)
    check close(fit.ssr, 4.480875822250005, 1e-10)

  test "bottom fixed at zero":
    let fit = fitLL4(ryegrassDose, ryegrassRootLength, fixed = [(Bottom, 0.0)])
    check fit.converged
    check fit.params[Bottom] == 0.0
    check fit.fixed == {Bottom, Asymmetry}  # Asymmetry is always fixed (at 1) in a 4PL
    check close(fit.params[Top], 7.85542760957709)
    check close(fit.params[LogEC50], 0.5136644098261173)
    check close(fit.params[Hill], -2.4703242877044502)
    check close(fit.ssr, 6.622819869524972, 1e-10)

suite "fitLL4: synthetic data":
  test "recovers exact parameters from noiseless data":
    # Truth is given with Bottom < Top: the 4PL is unchanged by swapping Bottom
    # and Top and negating Hill, and fits use the Bottom < Top form.
    for truth in [[0.0, 100.0, -6.0, 1.2], [0.0, 100.0, -7.5, -0.8], [10.0, 60.0, 1.0, -3.0]]:
      var dose, resp: seq[float]
      for k in 0 .. 9:
        let x = truth[2] - 2.25 + 0.5 * float(k)
        for _ in 1 .. 3:
          dose.add pow(10.0, x)
          resp.add truth[0] + (truth[1] - truth[0]) / (1.0 + pow(10.0, (truth[2] - x) * truth[3]))
      let fit = fitLL4(dose, resp)
      check fit.converged
      for param in Bottom .. Hill:
        check close(fit.params[param], truth[ord(param)], 1e-6)

  test "the mirrored parameterisation describes the same curve":
    let dose = [1e-9, 1e-8, 1e-7, 1e-6, 1e-5]
    var resp: seq[float]
    for d in dose:
      resp.add 100.0 + (0.0 - 100.0) / (1.0 + pow(10.0, (-7.5 - log10(d)) * 0.8))
    let fit = fitLL4(dose, resp)
    check fit.params[Bottom] < fit.params[Top]
    check close(fit.params[Hill], -0.8, 1e-6)
    for i in 0 ..< dose.len:
      check close(fit.predict(dose[i]), resp[i], 1e-9)

  test "recovers parameters from noisy replicates":
    var rng = initRand(42)
    let truth = [2.0, 98.0, -6.3, -1.4]
    var dose, resp: seq[float]
    for k in 0 .. 11:
      let x = -9.0 + 0.5 * float(k)
      for _ in 1 .. 4:
        dose.add pow(10.0, x)
        resp.add truth[0] + (truth[1] - truth[0]) / (1.0 + pow(10.0, (truth[2] - x) * truth[3])) +
                 rng.gauss(0.0, 3.0)
    let fit = fitLL4(dose, resp)
    check fit.converged
    check abs(fit.params[LogEC50] - truth[2]) < 0.05
    check abs(fit.params[Hill] - truth[3]) < 0.2
    check abs(fit.params[Bottom] - truth[0]) < 3.0
    check abs(fit.params[Top] - truth[1]) < 3.0

  test "weights of zero ignore observations":
    var dose = @ryegrassDose
    var resp = @ryegrassRootLength
    var weights = newSeq[float](dose.len)
    for w in weights.mitems: w = 1.0
    dose.add 3.0
    resp.add 1000.0  # wild outlier, given zero weight
    weights.add 0.0
    let fit = fitLL4(dose, resp, weights)
    check close(fit.params[LogEC50], 0.48543111513116)

  test "zero-weight observations do not affect the starting values":
    var dose = @ryegrassDose
    var resp = @ryegrassRootLength
    var weights = newSeq[float](dose.len)
    for w in weights.mitems: w = 1.0
    let clean = fitLL4(dose, resp, weights)
    dose.add 3.0
    resp.add 1e6
    weights.add 0.0
    let withOutlier = fitLL4(dose, resp, weights)
    check withOutlier.params == clean.params
    check withOutlier.iterations == clean.iterations

  test "invalid input is rejected":
    expect ValueError: discard fitLL4([1.0, 2.0], [1.0])
    expect ValueError: discard fitLL4([-1.0, 1.0, 2.0, 3.0], [1.0, 2.0, 3.0, 4.0])
    expect ValueError: discard fitLL4([1.0, 2.0, 3.0, 4.0], [1.0, NaN, 3.0, 4.0])
    expect ValueError: discard fitLL4([1.0, 2.0, 3.0], [1.0, 2.0, 3.0])
    check fitLL4([1.0, 2.0, 3.0], [1.0, 2.0, 3.0], fixed = [(Bottom, 0.0)]).n == 3
    expect ValueError: discard fitLL4([1.0, 2.0, 3.0, 4.0], [1.0, 2.0, 3.0, 4.0], fixed = [(Bottom, NaN)])
    expect ValueError: discard fitLL4([1.0, 2.0, 3.0, 4.0], [1.0, 2.0, 3.0, 4.0], fixed = [(Top, Inf)])
    expect ValueError: discard fitLL4([1.0, 2.0, 3.0, 4.0, 5.0], [1.0, 2.0, 3.0, 4.0, 5.0],
                                      weights = [0.0, 0.0, 0.0, 0.0, 0.0],
                                      fixed = [(Bottom, 0.0), (Top, 9.0), (LogEC50, 0.0), (Hill, 1.0)])
