import std/[math, sequtils, unittest]
import drc
import ryegrass
import fixtures

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

  test "replicates keep the settings":
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
    check b.fit.doses == @ryegrassDose

  test "replicates are ordinary fits to their resampled data":
    # Data-dependent results on a replicate match a fresh fit to its data.
    for resample in [Residuals, Cases]:
      for r in fit.bootstrap(n = 5, resample = resample).replicates:
        check r.doses.len == 24
        if resample == Residuals:
          check r.doses == @ryegrassDose
          check r.responses != @ryegrassRootLength
        let fresh = fitLL4(r.doses, r.responses, r.weights)
        check close(fresh.ssr, r.ssr, 1e-8)
        check close(r.residuals.mapIt(it * it).sum, r.ssr, 1e-10)
        let (a, b) = (r.profileCI(Hill), fresh.profileCI(Hill))
        check abs(a.lower) < Inf and abs(a.upper) < Inf
        check close(a.lower, b.lower, 1e-6)
        check close(a.upper, b.upper, 1e-6)
        check close(r.lackOfFit.f, fresh.lackOfFit.f, 1e-6)

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
  let (doses, curve) = potencyDesign()
  const response = potencyResponse
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

  test "joint replicates are ordinary fits to their resampled data":
    for r in fit.bootstrap(n = 3).replicates:
      # Residual resampling keeps the curve ids, here in curve order.
      var d, y: seq[float]
      for c in r.curves:
        d.add c.doses
        y.add c.responses
      check d == doses
      let fresh = fitLL4Multi(d, y, curve, shared = {Bottom, Top, Hill})
      let (a, b) =
        (r.relativePotencyProfileCI(0, 1), fresh.relativePotencyProfileCI(0, 1))
      check abs(a.lower) < Inf and abs(a.upper) < Inf
      check close(a.lower, b.lower, 1e-6)
      check close(a.upper, b.upper, 1e-6)
      check close(r.lackOfFit.f, fresh.lackOfFit.f, 1e-6)

  test "replicates keep the shared parameters":
    let boot = fit.bootstrap(n = 10)
    for r in boot.replicates:
      check r.shared == fit.shared
      check r.curves.len == 3

  test "robust joint fits refit with their loss and tuning":
    let tukey = fitLL4Multi(
      doses,
      response,
      curve,
      shared = {Bottom, Top, Hill},
      opts = FitOptions(loss: Tukey, tuning: 4.0),
    )
    let b = tukey.bootstrap(n = 10)
    check b.replicates.len > 0
    for r in b.replicates:
      check r.loss == Tukey and r.tuning == 4.0

  test "curves of a joint fit are bootstrapped through the joint fit":
    expect ValueError:
      discard fit.curves[0].bootstrap(n = 5)
    let boot = fit.bootstrap(n = 50)
    let ci = boot.bootCI(
      proc(f: MultiFit): float =
        f.curves[2].ec50
    )
    check ci.lower < fit.curves[2].ec50 and fit.curves[2].ec50 < ci.upper
