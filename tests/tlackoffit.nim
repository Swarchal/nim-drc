import std/[math, sequtils, unittest]
import drc
import ryegrass
import fixtures

# References from R: nls() fits in x = log10(dose), residuals(fit) and
# residuals(fit, type = "pearson"), and lack-of-fit tests as
# anova(fit, lm(y ~ factor(dose))) (weighted lm for the weighted fit, and
# lm(y ~ factor(curve):factor(dose)) against an indexed nls for the joint fit).

proc close(a, b: float, rel = 1e-5): bool =
  abs(a - b) <= rel * abs(b)

proc near(a, b, sigma: float): bool =
  ## Residuals agree to 1e-5·σ: R's nls stops ~1e-7 short of the optimum,
  ## which moves a small residual by more than a relative tolerance allows.
  abs(a - b) <= 1e-5 * sigma

proc checkTest(t: FTest, f: float, df1, df2: int, p: float) =
  check close(t.f, f)
  check t.df1 == df1
  check t.df2 == df2
  check close(t.p, p)

suite "stored data and residuals":
  let fit = fitLL4(ryegrassDose, ryegrassRootLength)

  test "the fit keeps its data":
    check fit.doses == @ryegrassDose
    check fit.responses == @ryegrassRootLength
    check fit.weights == newSeqWith(24, 1.0)

  test "residuals match R":
    let r = fit.residuals
    check r.len == 24
    check near(r[0], -0.21296230697004681, fit.sigma)
    check near(r[6], 0.77321870708744633, fit.sigma)
    check near(r[12], -1.1331767273421434, fit.sigma)
    check near(r[23], -0.049463954895298867, fit.sigma)
    let f = fit.fitted
    for i in 0 ..< 24:
      check close(f[i] + r[i], ryegrassRootLength[i], 1e-12)
      check close(f[i], fit.predict(ryegrassDose[i]), 1e-12)

  test "Pearson residuals match R":
    let r = fit.pearsonResiduals
    check near(r[0], -0.40983800539792026, 1.0)
    check near(r[6], 1.4880305212586236, 1.0)
    check near(r[12], -2.1807562864285068, 1.0)
    check near(r[23], -0.09519153366531323, 1.0)
    var ss = 0.0
    for v in r:
      ss += v * v
    check close(ss, float(fit.df), 1e-10)

  test "Pearson residuals without a usable σ":
    # As many free parameters as observations with weight: σ is NaN.
    let exact = fitLL4(
      [0.1, 1.0, 10.0, 100.0, 1000.0],
      [1.0, 3.0, 8.0, 9.0, 50.0],
      [1.0, 1.0, 1.0, 1.0, 0.0],
    )
    check exact.df == 0
    var r = exact.pearsonResiduals
    check r[4] == 0.0
    for v in r[0 .. 3]:
      check v.isNaN
    # An exact fit, σ = 0.
    var weights = newSeqWith(24, 1.0)
    weights[3] = 0.0
    var zero = fitLL4(ryegrassDose, ryegrassRootLength, weights)
    zero.sigma = 0.0
    r = zero.pearsonResiduals
    check r[3] == 0.0
    check r[0].isNaN

  test "weighted Pearson residuals match R, and zero weights give 0":
    var weights: seq[float]
    for i in 0 ..< 8:
      weights.add [1.0, 2.0, 0.5]
    let wfit = fitLL4(ryegrassDose, ryegrassRootLength, weights)
    check wfit.weights == weights
    let r = wfit.pearsonResiduals
    check near(r[0], -0.17313946885483877, 1.0)
    check near(r[6], 1.7122366411325041, 1.0)
    check near(r[12], -1.9375873485852111, 1.0)
    check near(r[23], -0.045249903392545303, 1.0)
    weights[3] = 0.0
    let zfit = fitLL4(ryegrassDose, ryegrassRootLength, weights)
    check zfit.pearsonResiduals[3] == 0.0
    check zfit.residuals[3] != 0.0

  test "curves of a joint fit keep their own data":
    let (doses, resp, curve) = ryegrassCurves(2.0)
    let joint = fitLL4Multi(doses, resp, curve, shared = {Bottom, Top, Hill})
    check joint.curves[1].doses == doses[24 .. 47]
    check joint.curves[1].responses == @ryegrassRootLength
    var ssr = 0.0
    for c in joint.curves:
      for r in c.residuals:
        ssr += r * r
    check close(ssr, joint.ssr, 1e-10)

suite "lack-of-fit test":
  test "ryegrass 4PL matches R":
    checkTest(
      fitLL4(ryegrassDose, ryegrassRootLength).lackOfFit,
      0.24106272392832281,
      3,
      17,
      0.86648302477503492,
    )

  test "weighted fit matches R":
    var weights: seq[float]
    for i in 0 ..< 8:
      weights.add [1.0, 2.0, 0.5]
    checkTest(
      fitLL4(ryegrassDose, ryegrassRootLength, weights).lackOfFit,
      0.41769202229753633,
      3,
      17,
      0.74257974897135515,
    )

  test "fixed parameters free a degree of freedom for lack of fit":
    checkTest(
      fitLL3(ryegrassDose, ryegrassRootLength).lackOfFit,
      1.1839267223897403,
      4,
      17,
      0.35316785450538896,
    )

  test "asymmetric data: the 4PL lacks fit, the 5PL does not":
    checkTest(
      fitLL4(asymDose, asymResponse).lackOfFit,
      3.9860458464630013,
      10,
      14,
      0.0095048095151010949,
    )
    checkTest(
      fitLL5(asymDose, asymResponse).lackOfFit,
      1.8022426165285323,
      9,
      14,
      0.15587556031686667,
    )

  test "zero-weight observations are left out":
    var dose = @ryegrassDose
    var resp = @ryegrassRootLength
    var w = newSeqWith(24, 1.0)
    dose.add [3.0, 3.0]
    resp.add [1000.0, -1000.0]
    w.add [0.0, 0.0]
    let t = fitLL4(dose, resp, w).lackOfFit
    checkTest(t, 0.24106272392832281, 3, 17, 0.86648302477503492)

  test "joint fit matches R":
    let (doses, curve) = potencyDesign()
    const response = potencyResponse
    let fit = fitLL4Multi(doses, response, curve, shared = {Bottom, Top, Hill})
    checkTest(fit.lackOfFit, 1.6517296515743194, 21, 27, 0.10926207189636181)

  test "identical replicates give no pure error":
    # Responses read to the instrument's resolution: each dose's triplicates
    # are equal, so there is no scatter to test the curve against.
    var doses, resp: seq[float]
    for (d, y) in [(0.1, 9.1), (1.0, 8.3), (10.0, 5.7), (100.0, 2.3), (1000.0, 1.1)]:
      doses.add [d, d, d]
      resp.add [y, y, y]
    expect ValueError:
      discard fitLL4(doses, resp).lackOfFit
    # Scatter at one dose is enough, and the identical doses add exactly 0.
    resp[1] = 9.2
    let t = fitLL4(doses, resp).lackOfFit
    check t.df1 == 1
    check t.df2 == 10
    check t.f > 0.0 and t.f < Inf
    # Joint fits pool the pure error over curves in the same way.
    let joint = fitLL4Multi(
      doses & doses, resp & resp, newSeqWith(15, 0) & newSeqWith(15, 1), shared = {Hill}
    )
    check joint.lackOfFit.f < Inf
    resp[1] = 9.1
    expect ValueError:
      discard fitLL4Multi(
        doses & doses,
        resp & resp,
        newSeqWith(15, 0) & newSeqWith(15, 1),
        shared = {Hill},
      ).lackOfFit

  test "unsuitable fits are rejected":
    # No replicates.
    expect ValueError:
      discard
        fitLL4([0.1, 1.0, 10.0, 100.0, 1000.0], [9.0, 8.0, 5.0, 2.0, 1.0]).lackOfFit
    # Replicates, but no more doses than parameters.
    expect ValueError:
      discard fitLL4(
        [0.1, 0.1, 1.0, 1.0, 10.0, 10.0, 100.0, 100.0],
        [9.0, 9.2, 7.0, 7.3, 3.0, 3.1, 1.0, 1.2],
      ).lackOfFit
    expect ValueError:
      discard fitLL4(ryegrassDose, ryegrassRootLength, opts = FitOptions(loss: Huber)).lackOfFit
    let (doses, resp, curve) = ryegrassCurves(2.0)
    let joint = fitLL4Multi(doses, resp, curve, shared = {Bottom, Top, Hill})
    expect ValueError:
      discard joint.curves[0].lackOfFit
