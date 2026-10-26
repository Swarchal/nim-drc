import std/[math, sequtils, unittest]
import drc
import ryegrass

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

const
  asymDose = [
    0.0, 0.0, 1e-09, 1e-09, 3.1622776601683795e-09, 3.1622776601683795e-09, 1e-08,
    1e-08, 3.162277660168379e-08, 3.162277660168379e-08, 1e-07, 1e-07,
    3.162277660168379e-07, 3.162277660168379e-07, 1e-06, 1e-06, 3.162277660168379e-06,
    3.162277660168379e-06, 1e-05, 1e-05, 3.1622776601683795e-05, 3.1622776601683795e-05,
    0.0001, 0.0001, 0.00031622776601683794, 0.00031622776601683794, 0.001, 0.001,
  ]
  asymResponse = [
    5.0, 5.6, 8.34, 7.1, 10.65, 9.57, 16.18, 18.74, 22.64, 22.38, 37.01, 36.75, 54.87,
    52.8, 75.55, 77.0, 86.17, 87.94, 89.67, 90.89, 90.97, 94.18, 92.39, 95.46, 95.3,
    94.61, 89.96, 93.92,
  ]

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
    var doses, resp: seq[float]
    var curve: seq[int]
    for c in 0 .. 1:
      for i, d in ryegrassDose:
        doses.add d * float(c + 1)
        resp.add ryegrassRootLength[i]
        curve.add c
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
    let response = [
      0.037, 0.07, 0.059, 0.157, 0.207, 0.161, 0.498, 0.508, 1.048, 1.013, 1.632, 1.587,
      1.832, 1.768, 2.002, 1.956, 1.987, 2.026, 0.138, 0.129, 0.269, 0.263, 0.574,
      0.492, 1.175, 1.148, 1.645, 1.593, 1.868, 1.904, 2.021, 1.963, 2.006, 1.988,
      1.942, 1.981, 0.036, 0.05, 0.101, 0.088, 0.069, 0.066, 0.167, 0.162, 0.306, 0.305,
      0.799, 0.815, 1.375, 1.415, 1.789, 1.749, 1.944, 1.885,
    ]
    var doses: seq[float]
    var curve: seq[int]
    for c in 0 .. 2:
      for k in 0 .. 8:
        let dose = pow(10.0, -2.0 + 0.5 * float(k))
        doses.add [dose, dose]
        curve.add [c, c]
    let fit = fitLL4Multi(doses, response, curve, shared = {Bottom, Top, Hill})
    checkTest(fit.lackOfFit, 1.6517296515743194, 21, 27, 0.10926207189636181)

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
    var doses, resp: seq[float]
    var curve: seq[int]
    for c in 0 .. 1:
      for i, d in ryegrassDose:
        doses.add d * float(c + 1)
        resp.add ryegrassRootLength[i]
        curve.add c
    let joint = fitLL4Multi(doses, resp, curve, shared = {Bottom, Top, Hill})
    expect ValueError:
      discard joint.curves[0].lackOfFit
