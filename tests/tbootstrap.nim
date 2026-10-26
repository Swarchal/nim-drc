import std/[math, sequtils, unittest]
import drc
import ryegrass

# Bootstrap draws cannot match R's random numbers, so the resampling is
# checked against the same scheme written independently in R with 20,000
# replicates (10,000 for the joint fit; nls() fails on 0.04–1% of them),
# allowing for Monte Carlo error: with 2000 replicates here, a percentile
# end has a standard error of about 0.06 bootstrap SDs, and an SD about 1.6%.
# The percentile rule itself is checked exactly against boot:::norm.inter,
# which boot.ci(type = "perc") uses.

proc close(a, b: float, rel = 1e-12): bool =
  abs(a - b) <= rel * abs(b)

proc sd(v: openArray[float]): float =
  let m = sum(v) / float(v.len)
  for x in v:
    result += (x - m) ^ 2
  sqrt(result / float(v.len - 1))

proc hill(f: DrcFit): float =
  f.params[Hill]

proc logEC50(f: DrcFit): float =
  f.params[LogEC50]

suite "bootstrap percentile intervals":
  proc withHills(values: seq[float]): Bootstrap[DrcFit] =
    for v in values:
      var f: DrcFit
      f.params[Hill] = v
      result.replicates.add f

  test "the percentile rule matches boot.ci":
    for (r, lo95, hi95, lo90, hi90) in [
      (
        39, -2.7485676957842871, 3.8913861588522636, -2.7179193376008284,
        3.8458894264756478,
      ),
      (
        50, -2.7788912705139435, 3.6674537493982466, -2.6396958149077956,
        3.6394784498404285,
      ),
      (
        10, -2.3767728239894153, 3.7680747398701451, -2.3767728239894153,
        3.7680747398701451,
      ),
    ]:
      var values: seq[float]
      for k in 1 .. r:
        values.add sin(float(k)) * 3.0 + float(k) / float(r)
      let boot = withHills(values)
      let ci95 = boot.bootCI(hill)
      let ci90 = boot.bootCI(hill, 0.9)
      check close(ci95.lower, lo95)
      check close(ci95.upper, hi95)
      check close(ci90.lower, lo90)
      check close(ci90.upper, hi90)

  test "non-finite values are left out":
    var values = @[NaN, Inf]
    for k in 1 .. 50:
      values.add sin(float(k)) * 3.0 + float(k) / 50.0
    let ci = withHills(values).bootCI(hill)
    check close(ci.lower, -2.7788912705139435)
    check close(ci.upper, 3.6674537493982466)
    check withHills(@[NaN]).bootCI(hill).lower.isNaN

  test "invalid levels are rejected":
    expect ValueError:
      discard withHills(@[1.0, 2.0]).bootCI(hill, 1.0)

suite "bootstrap: ryegrass 4PL":
  let fit = fitLL4(ryegrassDose, ryegrassRootLength)
  let boot = fit.bootstrap()

  test "residual resampling matches R within Monte Carlo error":
    check boot.replicates.len + boot.failed == 2000
    check boot.failed <= 10
    check boot.resample == Residuals
    let sdL = boot.bootValues(logEC50).sd
    check abs(sdL / 0.0267645 - 1.0) < 0.08
    check abs(boot.bootValues(hill).sd / 0.514065 - 1.0) < 0.08
    let ec = boot.bootCI(ec50)
    let sdEC = 0.189369
    check abs(ec.lower - 2.715317817) < 0.3 * sdEC
    check abs(ec.upper - 3.45906975) < 0.3 * sdEC
    let h = boot.bootCI(hill)
    check abs(h.lower - -4.232110773) < 0.3 * 0.514065
    check abs(h.upper - -2.226128184) < 0.3 * 0.514065

  test "case resampling matches R within Monte Carlo error":
    let cases = fit.bootstrap(resample = Cases)
    check cases.failed <= 40
    check abs(cases.bootValues(logEC50).sd / 0.0414714 - 1.0) < 0.08
    let ec = cases.bootCI(ec50)
    let sdEC = 0.298446
    check abs(ec.lower - 2.594911296) < 0.3 * sdEC
    check abs(ec.upper - 3.752459472) < 0.3 * sdEC

  test "a seed reproduces the resamples":
    let a = fit.bootstrap(n = 20, seed = 7)
    let b = fit.bootstrap(n = 20, seed = 7)
    let c = fit.bootstrap(n = 20, seed = 8)
    check a.bootValues(ec50) == b.bootValues(ec50)
    check a.bootValues(ec50) != c.bootValues(ec50)

  test "replicates keep the settings and drop the data":
    let fixed = fitLL4(
      ryegrassDose,
      ryegrassRootLength,
      fixed = [(Bottom, 0.0)],
      bounds = [(Hill, -10.0 .. -0.5)],
    )
    let b = fixed.bootstrap(n = 20)
    for r in b.replicates:
      check r.params[Bottom] == 0.0
      check r.model == LL4
      check r.lower[Hill] == -10.0
      check r.doses.len == 0
    check b.fit.doses == @ryegrassDose

  test "robust fits refit with their loss":
    let huber = fitLL4(ryegrassDose, ryegrassRootLength, opts = FitOptions(loss: Huber))
    let b = huber.bootstrap(n = 20)
    check b.replicates.len > 0
    for r in b.replicates:
      check r.loss == Huber

  test "zero-weight observations are left as they are":
    var weights = newSeqWith(24, 1.0)
    weights[3] = 0.0
    let w = fitLL4(ryegrassDose, ryegrassRootLength, weights)
    check w.bootstrap(n = 20).replicates.allIt(it.n == 23)
    check w.bootstrap(n = 20, resample = Cases).replicates.allIt(it.n == 23)

  test "invalid arguments are rejected":
    expect ValueError:
      discard fit.bootstrap(n = 0)
    expect ValueError:
      discard fitLL4([0.1, 1.0, 10.0, 100.0], [1.0, 3.0, 8.0, 9.0]).bootstrap(n = 5)

suite "bootstrap: joint fits":
  let response = [
    0.037, 0.07, 0.059, 0.157, 0.207, 0.161, 0.498, 0.508, 1.048, 1.013, 1.632, 1.587,
    1.832, 1.768, 2.002, 1.956, 1.987, 2.026, 0.138, 0.129, 0.269, 0.263, 0.574, 0.492,
    1.175, 1.148, 1.645, 1.593, 1.868, 1.904, 2.021, 1.963, 2.006, 1.988, 1.942, 1.981,
    0.036, 0.05, 0.101, 0.088, 0.069, 0.066, 0.167, 0.162, 0.306, 0.305, 0.799, 0.815,
    1.375, 1.415, 1.789, 1.749, 1.944, 1.885,
  ]
  var doses: seq[float]
  var curve: seq[int]
  for c in 0 .. 2:
    for k in 0 .. 8:
      let dose = pow(10.0, -2.0 + 0.5 * float(k))
      doses.add [dose, dose]
      curve.add [c, c]
  let fit = fitLL4Multi(doses, response, curve, shared = {Bottom, Top, Hill})
  proc potency(f: MultiFit): float =
    f.relativePotency(0, 1).ratio

  test "potency by residual resampling matches R within Monte Carlo error":
    let boot = fit.bootstrap()
    check boot.failed <= 20
    check abs(boot.bootValues(potency).sd / 0.167058 - 1.0) < 0.08
    let ci = boot.bootCI(potency)
    check abs(ci.lower - 3.493127288) < 0.3 * 0.167058
    check abs(ci.upper - 4.146510796) < 0.3 * 0.167058

  test "potency by case resampling matches R within Monte Carlo error":
    let boot = fit.bootstrap(resample = Cases)
    check boot.failed <= 40
    for r in boot.replicates:
      for c in 0 .. 2:
        check r.curves[c].n == 18
    let ci = boot.bootCI(potency)
    check abs(ci.lower - 3.457358658) < 0.3 * 0.188875
    check abs(ci.upper - 4.191799379) < 0.3 * 0.188875

  test "replicates keep the shared parameters":
    let boot = fit.bootstrap(n = 10)
    for r in boot.replicates:
      check r.shared == fit.shared
      check r.curves.len == 3

  test "curves of a joint fit cannot be bootstrapped alone":
    expect ValueError:
      discard fit.curves[0].bootstrap(n = 5)
