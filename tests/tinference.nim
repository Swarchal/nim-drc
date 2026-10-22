import std/[math, strutils, unittest]
import drc
import ryegrass

# Reference values from R: nls() with summary(), vcov(), logLik(), AIC() and
# BIC(), and t-based Wald intervals est ± qt(0.975, df) · se. ECx intervals use
# the delta method on log10(ECx) with vcov() from the same fit.

proc close(a, b: float, rel = 1e-5): bool =
  abs(a - b) <= rel * max(1e-3, abs(b))

proc checkInterval(ci: ConfInt, lower, upper: float) =
  check close(ci.lower, lower)
  check close(ci.upper, upper)

suite "inference: ryegrass, all observations":
  let fit = fitLL4(ryegrassDose, ryegrassRootLength)

  test "standard errors":
    check close(fit.se[Bottom], 0.209506043286517868)
    check close(fit.se[Top], 0.189733260164579726)
    check close(fit.se[LogEC50], 0.026391645145418065)
    check close(fit.se[Hill], 0.458428813539406621)

  test "covariance is symmetric with SE² on the diagonal":
    for i in DrcParam:
      check close(fit.cov[i][i], fit.se[i] ^ 2, 1e-12)
      for j in DrcParam:
        check fit.cov[i][j] == fit.cov[j][i]

  test "residual error and goodness of fit":
    check fit.n == 24
    check fit.df == 20
    check close(fit.sigma, 0.51962556952875094, 1e-9)
    check close(fit.rSquared, 0.97825213902727193, 1e-9)
    check close(fit.logLik, -16.15514326645846, 1e-9)
    check close(fit.aic, 42.310286532916919, 1e-9)
    check close(fit.bic, 48.200555684656649, 1e-9)

  test "parameter confidence intervals":
    checkInterval(fit.confint(Bottom), 0.044388009507346338, 0.91843190606129899)
    checkInterval(fit.confint(Top), 7.397185682059396861, 8.18873897292719199)
    checkInterval(fit.confint(LogEC50), 0.430379106349147444, 0.54048312052225489)
    checkInterval(fit.confint(Hill), -3.938494423592238114, -2.02596292715891657)

  test "ECx estimates and intervals":
    # Point estimates also match drc's ED(LL.4 fit, c(10, 50, 90)).
    check close(fit.ecx(10.0), 1.4637092229506501)
    check close(fit.ecx(50.0), fit.ec50, 1e-12)
    check close(fit.ecx(90.0), 6.3886253071017132)
    checkInterval(fit.ecxCI(10.0), 1.1182326614376683, 1.915920329671231)
    checkInterval(fit.ec50CI, 2.69388533671466, 3.4712278420078828)
    checkInterval(fit.ecxCI(90.0), 4.8916512198533351, 8.343712885516025)

  test "wider intervals at higher confidence":
    let ci95 = fit.confint(LogEC50, 0.95)
    let ci99 = fit.confint(LogEC50, 0.99)
    check ci99.lower < ci95.lower
    check ci99.upper > ci95.upper

  test "no warnings for a good fit":
    check fit.warnings == {}

  test "invalid levels and percentages are rejected":
    expect ValueError:
      discard fit.confint(Hill, 1.0)
    expect ValueError:
      discard fit.confint(Hill, 0.0)
    expect ValueError:
      discard fit.ecx(0.0)
    expect ValueError:
      discard fit.ecxCI(100.0)

  test "summary lists parameters and EC50":
    let text = fit.summary
    for word in [
      "4PL fit to 24 observations", "Bottom", "Top", "LogEC50", "Hill", "EC50",
      "20 degrees of freedom",
    ]:
      check word in text
    check "Asymmetry" notin text

suite "inference: ryegrass, bottom fixed at zero":
  let fit = fitLL4(ryegrassDose, ryegrassRootLength, fixed = [(Bottom, 0.0)])

  test "fixed parameters have zero standard error":
    check fit.df == 21
    check fit.se[Bottom] == 0.0
    for param in DrcParam:
      check fit.cov[Bottom][param] == 0.0
    check fit.confint(Bottom) == (0.0, 0.0)
    check "(fixed)" in fit.summary

  test "matches R with three free parameters":
    check close(fit.se[Top], 0.207984311964259566)
    check close(fit.se[LogEC50], 0.026909116126382086)
    check close(fit.se[Hill], 0.318735984154862040)
    check close(fit.sigma, 0.56158024781796267, 1e-9)
    check close(fit.aic, 45.208267459473952, 1e-9)
    check close(fit.bic, 49.920482780865733, 1e-9)
    checkInterval(fit.confint(Top), 7.4229005553810561, 8.28795466467512654)
    checkInterval(fit.confint(LogEC50), 0.4577038387356418, 0.56962497962725989)
    checkInterval(fit.confint(Hill), -3.1331720610854665, -1.80747653016275933)

suite "inference: ryegrass, weighted":
  var weights: seq[float]
  for i in 0 ..< 8:
    weights.add [1.0, 2.0, 0.5]
  let fit = fitLL4(ryegrassDose, ryegrassRootLength, weights)

  test "matches weighted R nls":
    check close(fit.params[Bottom], 0.46740217012759200)
    check close(fit.params[Top], 7.66539482963887142)
    check close(fit.params[LogEC50], 0.48273738992792814)
    check close(fit.params[Hill], -3.25604935792272920)
    check close(fit.se[Bottom], 0.176650541719044263)
    check close(fit.se[Top], 0.164253335766425368)
    check close(fit.se[LogEC50], 0.022236071428808096)
    check close(fit.se[Hill], 0.451503874110482428)
    check close(fit.sigma, 0.49321627622910563, 1e-9)
    check close(fit.rSquared, 0.98303346034124284, 1e-9)
    check close(fit.logLik, -14.903285948322626, 1e-9)
    check close(fit.aic, 39.806571896645252, 1e-9)
    check close(fit.bic, 45.696841048384982, 1e-9)
    checkInterval(fit.confint(LogEC50), 0.436353757716641588, 0.52912102213921464)

  test "zero-weight observations do not count":
    var dose = @ryegrassDose
    var resp = @ryegrassRootLength
    var w = weights
    dose.add 3.0
    resp.add 1000.0
    w.add 0.0
    let withExtra = fitLL4(dose, resp, w)
    check withExtra.n == 24
    check withExtra.df == 20
    check close(withExtra.aic, fit.aic, 1e-9)
    check close(withExtra.se[Hill], fit.se[Hill])

suite "inference: warnings":
  test "EC50 outside the tested doses":
    # Only the upper part of a curve with EC50 = 1e-6 is observed.
    var dose, resp: seq[float]
    for k in 0 .. 5:
      let d = pow(10.0, -5.5 + 0.5 * float(k))
      dose.add d
      resp.add 100.0 / (1.0 + pow(10.0, (-6.0 - log10(d)) * 1.0)) + [0.3, -0.2][k mod 2]
    let fit = fitLL4(dose, resp, fixed = [(Bottom, 0.0)])
    check fwEc50OutsideRange in fit.warnings
    check fit.params[LogEC50] < fit.logDoseRange.lower

  test "no residual degrees of freedom":
    let fit = fitLL4([0.1, 1.0, 10.0, 100.0], [1.0, 3.0, 8.0, 9.0])
    check fit.df == 0
    check fwNoResidualDf in fit.warnings
    check fit.sigma.isNaN
    check fit.se[Hill].isNaN
    check fit.confint(Hill).lower.isNaN

  test "flat response is not identifiable":
    let fit = fitLL4(
      [0.1, 1.0, 10.0, 100.0, 0.1, 1.0, 10.0, 100.0],
      [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0],
    )
    check fwSingular in fit.warnings
    check fit.se[LogEC50].isNaN
