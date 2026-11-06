import std/[math, strutils, unittest]
import drc
import drc/linalg
import ryegrass
import fixtures

proc close(a, b: float, rel: float): bool =
  abs(a - b) <= rel * max(1e-3, abs(b))

proc near(a, b, se: float, k = 2e-5): bool =
  ## Estimates agree to a small fraction of their standard error. Optimisers,
  ## R's nls included, stop short of the exact optimum by about this much
  ## along poorly determined directions (e.g. a Bottom near 0).
  abs(a - b) <= k * se

proc ll4Seq(x: float, p: seq[float], grad: var seq[float]): float =
  ## ll5 with Asymmetry = 1, through the runtime-sized interface.
  var q, g: array[drcParamCount, float]
  for j in 0 ..< 4:
    q[j] = p[j]
  q[ord(Asymmetry)] = 1.0
  result = ll5(x, q, g)
  for j in 0 ..< 4:
    grad[j] = g[j]

suite "runtime-sized linear algebra and Levenberg–Marquardt":
  test "seq and array routines agree exactly":
    let a = [[4.0, 1.0, 0.5], [1.0, 3.0, 0.2], [0.5, 0.2, 2.0]]
    var invA: array[3, array[3, float]]
    check invertSpd(a, invA)
    let s = @[@[4.0, 1.0, 0.5], @[1.0, 3.0, 0.2], @[0.5, 0.2, 2.0]]
    var invS: seq[seq[float]]
    check invertSpd(s, invS)
    for i in 0 ..< 3:
      for j in 0 ..< 3:
        check invS[i][j] == invA[i][j]
    var xa: array[3, float]
    var xs: seq[float]
    check solveSpd(a, [1.0, 2.0, 3.0], xa)
    check solveSpd(s, @[1.0, 2.0, 3.0], xs)
    for i in 0 ..< 3:
      check xs[i] == xa[i]

  test "a seq fit takes the same steps as the array fit":
    var x, w: seq[float]
    for d in ryegrassDose:
      x.add log10(d)
      w.add 1.0
    let start = [0.5, 7.8, 0.5, -3.0, 1.0]
    let fixed = levenbergMarquardt(
      ll5, x, ryegrassRootLength, w, start, [true, true, true, true, false]
    )
    let dynamic = levenbergMarquardt(
      ll4Seq, x, ryegrassRootLength, w, start[0 ..< 4], @[true, true, true, true]
    )
    check dynamic.iterations == fixed.iterations
    check dynamic.status == fixed.status
    check dynamic.ssr == fixed.ssr
    for j in 0 ..< 4:
      check dynamic.params[j] == fixed.params[j]

  test "the start may alias the variable receiving the result":
    var x, w: seq[float]
    for d in ryegrassDose:
      x.add log10(d)
      w.add 1.0
    var r: LMDynResult
    r.params = @[0.5, 7.8, 0.5, -3.0]
    r = levenbergMarquardt(
      ll4Seq, x, ryegrassRootLength, w, r.params, @[true, true, true, true]
    )
    check r.status.converged
    check close(r.ssr, 5.400214650161578, 1e-9)

suite "F distribution":
  test "upper tail matches R's pf":
    check close(fTail(2.13058, 6, 36), 0.073545187427460246, 1e-10)
    check close(fTail(0.5, 1, 10), 0.49564750438312, 1e-10)
    check close(fTail(10, 3, 200), 3.5948204444886642e-06, 1e-8)
    check close(fTail(40, 2, 5), 0.00083922361604267493, 1e-9)
    check close(fTail(1e-3, 4, 4), 0.999997007985024, 1e-12)
    check fTail(0.0, 3, 4) == 1.0
    check fTail(-1.0, 3, 4) == 1.0
    check fTail(Inf, 3, 4) == 0.0
    check fTail(1.0, 0, 4).isNaN

suite "joint fits match R nls (DNase)":
  # nls(density ~ b + (t-b)/(1+10^((l[c]-x)*h))) and variants, x = log10(conc).
  let (doses, curve) = dnase()

  test "parallel curves: shared Bottom, Top and Hill":
    let fit = fitLL4Multi(doses, dnaseDensity, curve, shared = {Bottom, Top, Hill})
    check fit.converged
    check fit.warnings == {}
    check fit.n == 48 and fit.df == 42
    check fit.n - fit.df == 6 # Bottom, Top, 3 × LogEC50, Hill
    let
      expected = [
        0.012153083145364535, 2.2960210356810506, 1.0022379964196355,
        0.62687500591840939, 0.6176929850472922, 0.59580357999667788,
      ]
      se = [
        0.0095714474069814141, 0.054224708409959126, 0.030384658860899121,
        0.022528818528393384, 0.022572490292757402, 0.022679810945941136,
      ]
      vcov = [
        9.1612605464611231e-05, -0.00031578470058738232, 0.00022747314762914495,
        -9.711885393838455e-05, -9.7327030488189039e-05, -9.784555424271014e-05,
        -0.00031578470058738232, 0.0029403190021450917, -0.0015127839420201962,
        0.0011610676508598259, 0.001163657218109401, 0.0011699882688050098,
        0.00022747314762914495, -0.0015127839420201962, 0.00092322749409321555,
        -0.00057269800789644387, -0.00057376870066595874, -0.00057642141848427016,
        -9.711885393838455e-05, 0.0011610676508598259, -0.00057269800789644387,
        0.00050754766428528109, 0.00047327285564801404, 0.00047581224280361515,
        -9.7327030488189039e-05, 0.001163657218109401, -0.00057376870066595874,
        0.00047327285564801404, 0.00050951731801662716, 0.00047686158454947526,
        -9.784555424271014e-05, 0.0011699882688050098, -0.00057642141848427016,
        0.00047581224280361515, 0.00047686158454947526, 0.00051437382454363144,
      ]
      # R's coefficient order b, t, h, l1, l2, l3 as (curve, parameter)
      coefs =
        [(0, Bottom), (0, Top), (0, Hill), (0, LogEC50), (1, LogEC50), (2, LogEC50)]
    for i, (c, param) in coefs:
      check near(fit.param(c, param), expected[i], se[i])
      check close(sqrt(fit.cov(c, param, c, param)), se[i], 1e-5)
      for j, (d, other) in coefs:
        let v = fit.cov(c, param, d, other)
        check abs(v - vcov[6 * i + j]) <= 1e-5 * abs(vcov[6 * i + j])
    check close(fit.ssr, 0.017877581773117615, 1e-8)
    check close(fit.aic, -228.76152210258857, 1e-8)
    check close(fit.bic, -215.66311502623333, 1e-8)
    # The views carry the shared values and their uncertainty.
    for c in 0 .. 2:
      check fit.curves[c].params[Hill] == fit.param(0, Hill)
      check close(fit.curves[c].params[LogEC50], expected[3 + c], 1e-6)
      check close(fit.curves[c].se[Top], se[1], 1e-5)
      check fit.curves[c].df == 42
      check fit.curves[c].sigma == fit.sigma
      check fit.curves[c].n == 16

  test "F test of parallel against separate curves matches anova":
    let parallel = fitLL4Multi(doses, dnaseDensity, curve, shared = {Bottom, Top, Hill})
    let separate = fitLL4Multi(doses, dnaseDensity, curve)
    check separate.converged
    check separate.df == 36
    check close(separate.ssr, 0.013192846017506888, 1e-8)
    check close(separate.aic, -231.34741644049552, 1e-8)
    let t = fTest(parallel, separate)
    check t.df1 == 6 and t.df2 == 36
    check close(t.f, 2.1305800504579935, 1e-6)
    check close(t.p, 0.073545181288113118, 1e-6)

  test "relative potency uses the covariance shared through Bottom and Top":
    # Bottom and Top shared, Hill and LogEC50 per curve. R: ECx found by root
    # finding on the fitted curves, gradient by central differences.
    let fit = fitLL4Multi(doses, dnaseDensity, curve, shared = {Bottom, Top})
    check fit.converged
    check near(fit.param(0, Bottom), 0.011673323618921537, 0.0092397605852578713)
    check close(fit.curves[1].params[Hill], 1.0270339163214497, 1e-6)
    check close(fit.curves[2].params[LogEC50], 0.5982597126486634, 1e-6)
    check close(fit.curves[2].se[Hill], 0.030988720959629807, 1e-5)
    let rp50 = fit.relativePotency(0, 2)
    check close(rp50.ratio, 1.0735395579650393, 1e-6)
    check close(rp50.ci.lower, 1.0326652380824977, 1e-6)
    check close(rp50.ci.upper, 1.1160317400203823, 1e-6)
    let rp20 = fit.relativePotency(0, 2, pct = 20.0)
    check close(rp20.ratio, 1.0801053914170189, 1e-6)
    check close(rp20.ci.lower, 1.0158737488843668, 1e-6)
    check close(rp20.ci.upper, 1.1483982707982194, 1e-6)
    # The ratio of a curve to itself is exactly 1 with no uncertainty.
    let same = fit.relativePotency(1, 1)
    check same.ratio == 1.0 and same.ci.lower == 1.0 and same.ci.upper == 1.0

  test "unequal curve lengths, weights and zero weights":
    # Three rows of run 4 dropped; R fitted without the zero-weight rows.
    const
      wConc = [
        0.04882812, 0.04882812, 0.1953125, 0.1953125, 0.390625, 0.390625, 0.78125,
        0.78125, 1.5625, 1.5625, 3.125, 3.125, 6.25, 6.25, 12.5, 12.5, 0.04882812,
        0.1953125, 0.390625, 0.390625, 0.78125, 0.78125, 1.5625, 3.125, 3.125, 6.25,
        6.25, 12.5, 12.5, 0.04882812, 0.04882812, 0.1953125, 0.1953125, 0.390625,
        0.390625, 0.78125, 0.78125, 1.5625, 1.5625, 3.125, 3.125, 6.25, 6.25, 12.5, 12.5,
      ]
      wDensity = [
        0.017, 0.018, 0.121, 0.124, 0.206, 0.215, 0.377, 0.374, 0.614, 0.609, 1.019,
        1.001, 1.334, 1.364, 1.73, 1.71, 0.016, 0.118, 0.2, 0.206, 0.364, 0.36, 0.64,
        0.979, 0.973, 1.424, 1.399, 1.74, 1.732, 0.054, 0.054, 0.152, 0.148, 0.226,
        0.222, 0.392, 0.383, 0.658, 0.644, 1.043, 1.002, 1.466, 1.381, 1.743, 1.724,
      ]
      wCurve = [
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
      ]
      wWeights = [
        1.0, 2, 0.5, 1, 0, 0.5, 1, 2, 0.5, 1, 2, 0.5, 1, 2, 0.5, 1, 2, 0.5, 1, 2, 0.5,
        1, 2, 0.5, 1, 2, 0.5, 1, 2, 0, 1, 2, 0.5, 1, 2, 0.5, 1, 2, 0.5, 1, 2, 0.5, 1, 2,
        0.5,
      ]
    let fit = fitLL4Multi(
      wConc, wDensity, wCurve, shared = {Bottom, Top, Hill}, weights = wWeights
    )
    check fit.n == 43 and fit.df == 37
    check fit.curves[1].n == 13
    check near(fit.param(0, Bottom), 0.0059287617775354053, 0.010642068841095056)
    check close(fit.param(0, Top), 2.3135156470532325, 1e-6)
    check close(fit.param(0, Hill), 0.98505516143174987, 1e-6)
    check close(fit.curves[0].params[LogEC50], 0.63013919556039821, 1e-6)
    check close(fit.curves[1].params[LogEC50], 0.61663361491202162, 1e-6)
    check close(fit.curves[2].params[LogEC50], 0.60400232702416579, 1e-6)
    check close(fit.curves[1].se[LogEC50], 0.02423132344325853, 1e-5)
    check close(fit.curves[0].se[Top], 0.057269565417712109, 1e-5)
    check close(fit.ssr, 0.018301583379243987, 1e-8)
    check close(fit.aic, -197.73590219435633, 1e-8)

  test "3PL: Bottom fixed at 0":
    let fit = fitLL3Multi(doses, dnaseDensity, curve, shared = {Top, Hill})
    check fit.converged
    check fit.df == 43
    check Bottom in fit.fixed and Bottom in fit.shared
    check fit.param(0, Bottom) == 0.0
    check close(fit.param(0, Top), 2.3404023237518992, 1e-6)
    check close(fit.param(0, Hill), 0.97194346869682302, 1e-6)
    check close(fit.curves[1].params[LogEC50], 0.6316735935066734, 1e-6)
    check close(fit.curves[2].se[LogEC50], 0.021523777346403494, 1e-5)
    check fit.curves[0].se[Bottom] == 0.0
    check close(fit.ssr, 0.018515062405032871, 1e-8)

suite "joint fits match R nls (simulated)":
  test "zero-dose controls with shared Bottom, Top and Hill":
    # Two falling curves, controls at dose 0 (x = -Inf in R as here).
    const y = [
      98.12, 100.55, 97.49, 103.52, 99.72, 93.2, 97.12, 88.42, 87.94, 63.36, 68.82,
      36.44, 33.41, 9.42, 19.44, 8.27, 8.35, 8.82, 8.45, 7.06, 8.04, 102.35, 100.22,
      94.03, 101.58, 99.55, 98.55, 94.6, 95.16, 97.85, 93.01, 88.63, 70.89, 69.56,
      36.59, 39.47, 17.61, 18.61, 12.64, 11.63, 5.77, 5.5,
    ]
    var doses: seq[float]
    var curve: seq[int]
    for c in 0 .. 1:
      for _ in 1 .. 3:
        doses.add 0.0
        curve.add c
      for k in 0 .. 8: # 10^seq(-8, -4, by = 0.5), duplicates
        for _ in 1 .. 2:
          doses.add pow(10.0, -8.0 + 0.5 * float(k))
          curve.add c
    let fit = fitLL4Multi(doses, y, curve, shared = {Bottom, Top, Hill})
    check fit.converged
    check close(fit.param(0, Bottom), 7.13590098008161, 1e-6)
    check close(fit.param(0, Top), 99.251720675476605, 1e-6)
    check close(fit.param(0, Hill), -1.2430191739060219, 1e-6)
    check close(fit.curves[0].params[LogEC50], -6.3016157818636929, 1e-6)
    check close(fit.curves[1].params[LogEC50], -5.727616082086584, 1e-6)
    check close(fit.curves[0].se[Bottom], 0.87537809766560548, 1e-5)
    check close(fit.curves[1].se[LogEC50], 0.022493025100516997, 1e-5)
    check close(fit.ssr, 223.74695123800646, 1e-8)
    # Controls sit on the low-dose plateau of each curve.
    check close(fit.curves[1].predict(0.0), fit.param(0, Top), 1e-12)

  test "5PL with shared Asymmetry":
    const y = [
      4.75, 6.37, 10.49, 6.41, 11.92, 12.24, 20.71, 19.29, 36.76, 33.58, 56.72, 57.56,
      78.13, 77.16, 93.04, 86.9, 95.21, 93.95, 96.27, 95.4, 98.08, 93.14, 97.37, 97.92,
      95, 91.32, 5.65, 4.04, 7.02, 6.26, 8.7, 8.07, 12.64, 10.6, 16.55, 16.82, 27.7,
      28.93, 50.31, 50.78, 74.32, 71.96, 87.66, 91.78, 94.45, 96.5, 94.2, 94.53, 94.65,
      93.13, 93.73, 98.08,
    ]
    var doses: seq[float]
    var curve: seq[int]
    for c in 0 .. 1:
      for k in 0 .. 12: # 10^seq(-9, -3, by = 0.5), duplicates
        for _ in 1 .. 2:
          doses.add pow(10.0, -9.0 + 0.5 * float(k))
          curve.add c
    let fit = fitLL5Multi(doses, y, curve, shared = {Bottom, Top, Hill, Asymmetry})
    check fit.converged
    check fit.n - fit.df == 6 # 4 shared, LogEC50 for each curve
    check close(fit.param(0, Bottom), 4.3933639193057008, 1e-5)
    check close(fit.param(0, Top), 95.575022100491481, 1e-6)
    check close(fit.param(0, Hill), 1.2468468905757522, 1e-5)
    check close(fit.param(0, Asymmetry), 0.479746422442341, 1e-5)
    check close(fit.curves[0].params[LogEC50], -6.2267271894299414, 1e-6)
    check close(fit.curves[1].params[LogEC50], -5.588832883366563, 1e-6)
    check close(fit.curves[0].se[Asymmetry], 0.084518932668916705, 1e-4)
    check close(fit.curves[1].se[LogEC50], 0.083226295801133041, 1e-4)
    check close(fit.ssr, 137.01299476247908, 1e-8)
    check close(fit.aic, 211.94887431481337, 1e-8)
    # Is the asymmetry real? anova(4PL, 5PL) in R.
    let symmetric = fitLL4Multi(doses, y, curve, shared = {Bottom, Top, Hill})
    check close(symmetric.ssr, 180.30455549333564, 1e-8)
    let t = fTest(symmetric, fit)
    check t.df1 == 1 and t.df2 == 46
    check close(t.f, 14.534473880172047, 1e-6)
    check close(t.p, 0.00040762775655688092, 1e-5)

suite "joint fits: consistency":
  let (doses, curve) = dnase()

  test "one curve is the single-curve fit":
    var ids = newSeq[int](ryegrassDose.len)
    let joint = fitLL4Multi(ryegrassDose, ryegrassRootLength, ids)
    let single = fitLL4(ryegrassDose, ryegrassRootLength)
    for param in Bottom .. Hill:
      check near(
        joint.curves[0].params[param], single.params[param], single.se[param], 1e-5
      )
      check close(joint.curves[0].se[param], single.se[param], 1e-6)
    check close(joint.ssr, single.ssr, 1e-12)
    check close(joint.aic, single.aic, 1e-12)

  test "sharing every parameter is the pooled fit":
    let joint =
      fitLL4Multi(doses, dnaseDensity, curve, shared = {Bottom, Top, LogEC50, Hill})
    let pooled = fitLL4(doses, dnaseDensity)
    check joint.n - joint.df == 4
    for param in Bottom .. Hill:
      check near(
        joint.curves[2].params[param], pooled.params[param], pooled.se[param], 1e-5
      )
      check close(joint.curves[2].se[param], pooled.se[param], 1e-6)
    check close(joint.ssr, pooled.ssr, 1e-12)

  test "sharing nothing gives each curve's own fit, with a pooled σ":
    let joint = fitLL4Multi(doses, dnaseDensity, curve)
    var totalSsr = 0.0
    for c in 0 .. 2:
      let single =
        fitLL4(doses[16 * c ..< 16 * (c + 1)], dnaseDensity[16 * c ..< 16 * (c + 1)])
      for param in Bottom .. Hill:
        check near(
          joint.curves[c].params[param], single.params[param], single.se[param], 1e-5
        )
      check close(joint.curves[c].ssr, single.ssr, 1e-8)
      totalSsr += single.ssr
    check close(joint.ssr, totalSsr, 1e-8)
    check close(joint.sigma, sqrt(totalSsr / 36.0), 1e-8)

  test "the per-curve views support single-curve inference":
    let fit = fitLL4Multi(doses, dnaseDensity, curve, shared = {Bottom, Top, Hill})
    let view = fit.curves[1]
    let ci = view.ec50CI
    check ci.lower < view.ec50 and view.ec50 < ci.upper
    check close(view.ec50, pow(10.0, view.params[LogEC50]), 1e-12)
    let band = view.confidenceBand(1.0)
    check band.fit == view.predict(1.0)
    check band.lower < band.fit and band.fit < band.upper
    check view.logDoseRange.lower == log10(0.04882812)

  test "summary":
    let fit = fitLL4Multi(doses, dnaseDensity, curve, shared = {Bottom, Top, Hill})
    let text = fit.summary
    check "4PL joint fit to 48 observations in 3 curves; shared: Bottom, Top, Hill" in
      text
    check "LogEC50[2]" in text
    check "EC50[0]" in text
    check "Residual standard error" in text
    check "Asymmetry" notin text

suite "joint fits: input checks":
  let (doses, curve) = dnase()

  test "invalid input is rejected":
    expect ValueError:
      discard fitLL4Multi(doses, dnaseDensity, curve[0 ..< 10])
    expect ValueError:
      discard fitLL4Multi([1.0, 2.0], [1.0, 2.0], [0, -1])
    var gap = curve
    for c in gap.mitems:
      if c == 1:
        c = 3 # ids 0, 2, 3: curve 1 is missing
    expect ValueError:
      discard fitLL4Multi(doses, dnaseDensity, gap)
    var w = newSeq[float](doses.len)
    for i in 0 ..< w.len:
      w[i] = (if curve[i] == 2: 0.0 else: 1.0)
    expect ValueError:
      discard fitLL4Multi(doses, dnaseDensity, curve, weights = w)
    expect ValueError:
      discard fitLL4Multi(
        doses, dnaseDensity, curve, opts = FitOptions(loss: Huber, tuning: -1.0)
      )
    expect ValueError:
      discard fitLL4Multi(doses, dnaseDensity, curve, fixed = [(Asymmetry, 2.0)])

  test "each curve needs as many points as its own parameters":
    let few = [0, 0, 0, 0, 0, 0, 1, 1, 1]
    let d = [0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 0.1, 1.0, 10.0]
    let r = [1.0, 1.2, 2.0, 3.5, 4.0, 4.1, 1.1, 2.2, 4.0]
    expect ValueError:
      discard fitLL4Multi(d, r, few) # curve 1: 3 points, 4 parameters
    let fit = fitLL4Multi(d, r, few, shared = {Bottom, Top})
    check fit.df == 9 - 6

  test "a curve of controls only leaves its LogEC50 unidentifiable":
    var d = doses
    var ids = curve
    var r = @dnaseDensity
    for _ in 1 .. 3:
      d.add 0.0
      ids.add 3
      r.add 0.02
    let fit = fitLL4Multi(d, r, ids, shared = {Bottom, Top, Hill})
    check fwSingular in fit.warnings
    check fit.curves[3].se[LogEC50].isNaN

  test "F test checks that the fits are nested on the same data":
    let parallel = fitLL4Multi(doses, dnaseDensity, curve, shared = {Bottom, Top, Hill})
    let separate = fitLL4Multi(doses, dnaseDensity, curve)
    expect ValueError:
      discard fTest(separate, parallel)
    let common = fitLL4Multi(doses, dnaseDensity, curve, shared = {Top, Hill})
    let other =
      fitLL4Multi(doses, dnaseDensity, curve, shared = {Bottom, Hill, LogEC50})
    expect ValueError:
      discard fTest(other, common) # does not share Top
    var shifted = @dnaseDensity
    shifted[0] += 0.001
    let elsewhere = fitLL4Multi(doses, shifted, curve)
    expect ValueError:
      discard fTest(parallel, elsewhere)
    # Pooled data against parallel curves: different curve ids, same data.
    let pooled = fitLL4Multi(doses, dnaseDensity, newSeq[int](doses.len))
    let t = fTest(pooled, parallel)
    check t.df1 == 2 and t.p < 0.05

suite "MultiFit accessors":
  let (doses, curve) = dnase()
  let fit = fitLL4Multi(doses, dnaseDensity, curve, shared = {Bottom, Top, Hill})

  test "fitMulti takes the model as an argument":
    let direct = fitMulti(LL4, doses, dnaseDensity, curve, shared = {Bottom, Top, Hill})
    check direct.ssr == fit.ssr
    for c in 0 .. 2:
      check direct.curves[c].params == fit.curves[c].params

  test "param and cov agree with the per-curve views":
    for c in 0 .. 2:
      for p in DrcParam:
        check fit.param(c, p) == fit.curves[c].params[p]
        for q in DrcParam:
          check fit.cov(c, p, c, q) == fit.curves[c].cov[p][q]

  test "cov spans curves":
    # Shared parameters are one slot, so across curves their covariance is the variance.
    check fit.cov(0, Hill, 2, Hill) == fit.cov(1, Hill, 1, Hill)
    check fit.cov(0, LogEC50, 1, LogEC50) != 0.0
    check fit.cov(0, LogEC50, 1, LogEC50) == fit.cov(1, LogEC50, 0, LogEC50)
    check fit.cov(0, Asymmetry, 1, Hill) == 0.0 # fixed

  test "a curve that does not exist raises":
    expect ValueError:
      discard fit.param(3, Hill)
    expect ValueError:
      discard fit.cov(0, Hill, -1, Hill)

suite "robust joint fits match robustbase::nlrob":
  # The README's potency data with gross errors added at observations 5, 22
  # and 47; parallel curves. R: nlrob() with the indexed formula written out
  # (b + (t-b)/(1+10^((l0*(c==0)+l1*(c==1)+l2*(c==2)-x)*h))), started from
  # nls() for Huber and from a Huber fit at nlrob's default tolerance for
  # bisquare, both run to tol = 1e-12. Standard errors follow MASS::rlm's
  # covariance, computed in R from a finite-difference Jacobian at the
  # nlrob fit (as tirls.nim).
  let (doses, curve) = potencyDesign()
  var response = @potencyResponse
  response[5] += 0.6
  response[22] -= 0.5
  response[47] += 0.6
  let
    huber = fitLL4Multi(
      doses,
      response,
      curve,
      shared = {Bottom, Top, Hill},
      opts = FitOptions(loss: Huber),
    )
    tukey = fitLL4Multi(
      doses,
      response,
      curve,
      shared = {Bottom, Top, Hill},
      opts = FitOptions(loss: Tukey),
    )
    # R's coefficient order b, t, h, l0, l1, l2 as (curve, parameter)
    coefs = [(0, Bottom), (0, Top), (0, Hill), (0, LogEC50), (1, LogEC50), (2, LogEC50)]

  test "Huber":
    check huber.converged
    check huber.warnings == {}
    check huber.loss == Huber and huber.tuning == 1.345
    let
      expected = [
        0.062090094158708059, 1.9976807979794629, 1.1042505388247708,
        -0.014497824016666246, -0.5813732251833873, 0.69134185114528857,
      ]
      se = [
        0.012114911583969381, 0.012460799742694534, 0.032703034698909401,
        0.017972175667133758, 0.017943567860478588, 0.017928420067593331,
      ]
    for i, (c, param) in coefs:
      check close(huber.param(c, param), expected[i], 2e-5)
      check close(sqrt(huber.cov(c, param, c, param)), se[i], 2e-5)
    check close(huber.cov(0, LogEC50, 1, LogEC50), 8.3784727733646073e-05, 2e-5)
    check close(huber.sigma, 0.034384829058765283, 2e-5)
    for (i, rw) in [
      (3, 0.89632501842134193),
      (5, 0.083651502066434541),
      (22, 0.095453255431882197),
      (47, 0.075077240610943508),
    ]:
      check close(huber.robustWeights[i], rw, 2e-5)

  test "Tukey rejects the gross errors":
    check tukey.converged
    let
      expected = [
        0.060824822915289144, 1.99836551884469, 1.1026026211419602,
        -0.011923882480351143, -0.59272372873542711, 0.70115412862418891,
      ]
      se = [
        0.011224534946907349, 0.011547616490340406, 0.030194275100864858,
        0.016651858014354825, 0.01662432277325401, 0.01660974894374281,
      ]
    for i, (c, param) in coefs:
      check close(tukey.param(c, param), expected[i], 1e-5)
      check close(sqrt(tukey.cov(c, param, c, param)), se[i], 1e-5)
    check close(tukey.sigma, 0.032876363476044078, 1e-5)
    var rejected: seq[int]
    for i, rw in tukey.robustWeights:
      if rw == 0.0:
        rejected.add i
    check rejected == @[5, 22, 47]
    check close(tukey.robustWeights[13], 0.40872327465844011, 1e-5)

  test "the curves carry the robust scale, weights and covariance":
    for c in 0 .. 2:
      let view = huber.curves[c]
      check view.loss == Huber and view.tuning == 1.345
      check view.sigma == huber.sigma
      check view.robustWeights == huber.robustWeights[18 * c ..< 18 * (c + 1)]
      check view.se[LogEC50] == sqrt(huber.cov(c, LogEC50, c, LogEC50))
    let ci = huber.curves[1].ec50CI
    check ci.lower < huber.curves[1].ec50 and huber.curves[1].ec50 < ci.upper

  test "no likelihood-based criteria or F tests":
    check huber.logLik.isNaN and huber.aic.isNaN and huber.bic.isNaN
    let separate = fitLL4Multi(doses, response, curve, opts = FitOptions(loss: Huber))
    expect ValueError:
      discard fTest(huber, separate)
    expect ValueError:
      discard huber.lackOfFit

  test "summary":
    let text = tukey.summary
    check "robust (Tukey bisquare, tuning 4.685); shared: Bottom, Top, Hill" in text
    check "Robust residual scale (MAD)" in text
    check "0 (rejected): 3" in text
    check "AIC" notin text

suite "joint fits: Nelder–Mead fallback":
  let (doses, curve) = potencyDesign()
  let reference =
    fitLL4Multi(doses, potencyResponse, curve, shared = {Bottom, Top, Hill})

  test "recovers when Levenberg–Marquardt does not converge":
    let starved = FitOptions(lm: LMOptions(maxIter: 2), multiStart: false)
    let fit = fitLL4Multi(
      doses, potencyResponse, curve, shared = {Bottom, Top, Hill}, opts = starved
    )
    check fit.usedFallback and fit.fallbackConverged and fit.converged
    check fit.curves[0].usedFallback
    check close(fit.ssr, reference.ssr, 1e-10)
    for c in 0 .. 2:
      for param in Bottom .. Hill:
        check near(
          fit.param(c, param),
          reference.param(c, param),
          sqrt(reference.cov(c, param, c, param)),
          1e-5,
        )
    check "Nelder–Mead fallback converged" in fit.summary

  test "without the fallback, non-convergence is reported":
    var starved = FitOptions(lm: LMOptions(maxIter: 2), multiStart: false)
    starved.fallback = false
    let fit = fitLL4Multi(
      doses, potencyResponse, curve, shared = {Bottom, Top, Hill}, opts = starved
    )
    check not fit.usedFallback
    check not fit.converged
    check fwNotConverged in fit.warnings

  test "is not used when Levenberg–Marquardt converges":
    check not reference.usedFallback
