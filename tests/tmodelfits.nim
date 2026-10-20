import std/[math, strutils, unittest]
import drc
import ryegrass

# 5PL reference values from R nls() on synthetic asymmetric data (truth:
# Bottom 5, Top 95, LogEC50 −6, Hill 1.3, Asymmetry 0.35, noise sd 2), which
# agree with scipy.optimize.curve_fit. ECx intervals use the delta method with
# finite-difference gradients in R, independent of the analytic ones here.
const
  asymDose = [0.0, 0.0, 1e-09, 1e-09, 3.1622776601683795e-09, 3.1622776601683795e-09,
              1e-08, 1e-08, 3.162277660168379e-08, 3.162277660168379e-08, 1e-07, 1e-07,
              3.162277660168379e-07, 3.162277660168379e-07, 1e-06, 1e-06,
              3.162277660168379e-06, 3.162277660168379e-06, 1e-05, 1e-05,
              3.1622776601683795e-05, 3.1622776601683795e-05, 0.0001, 0.0001,
              0.00031622776601683794, 0.00031622776601683794, 0.001, 0.001]
  asymResponse = [5.0, 5.6, 8.34, 7.1, 10.65, 9.57, 16.18, 18.74, 22.64, 22.38, 37.01, 36.75,
                  54.87, 52.8, 75.55, 77.0, 86.17, 87.94, 89.67, 90.89, 90.97, 94.18, 92.39,
                  95.46, 95.3, 94.61, 89.96, 93.92]

proc close(a, b: float; rel = 1e-5): bool = abs(a - b) <= rel * abs(b)

proc checkInterval(ci: ConfInt; lower, upper: float) =
  check close(ci.lower, lower)
  check close(ci.upper, upper)

suite "5PL":
  let fit = fitLL5(asymDose, asymResponse)

  test "matches R nls":
    check fit.converged
    check fit.model == LL5
    check close(fit.params[Bottom], 4.88038610101201797)
    check close(fit.params[Top], 93.22614340859124127)
    check close(fit.params[LogEC50], -6.08146193945548585)
    check close(fit.params[Hill], 1.21810555178779012)
    check close(fit.params[Asymmetry], 0.39201413304850541)
    check close(fit.se[Bottom], 1.017247025067105293)
    check close(fit.se[Top], 0.580843140471135500)
    check close(fit.se[LogEC50], 0.090944639884384004)
    check close(fit.se[Hill], 0.144525766320753118)
    check close(fit.se[Asymmetry], 0.076409480367067448)
    check fit.df == 23
    check close(fit.sigma, 1.6304169158252952, 1e-9)
    check close(fit.rSquared, 0.99831419622159212, 1e-9)
    check close(fit.aic, 113.32747210739997, 1e-9)
    check close(fit.bic, 121.3206991684512, 1e-9)

  test "ECx is the dose of the given response, not 10^LogEC50":
    check close(fit.ecx(10.0), 6.6890636873642421e-09)
    check close(fit.ec50, 2.2637989538232376e-07)
    check close(fit.ecx(90.0), 2.1777365877654262e-06)
    check abs(log10(fit.ec50) - fit.params[LogEC50]) > 0.1
    let half = 0.5 * (fit.params[Bottom] + fit.params[Top])
    check close(fit.predict(fit.ec50), half, 1e-10)
    checkInterval(fit.ecxCI(10.0), 3.8859944314637203e-09, 1.151405999229946e-08)
    checkInterval(fit.ec50CI, 2.0174232507999178e-07, 2.5402630317159192e-07)
    checkInterval(fit.ecxCI(90.0), 1.6846299715056307e-06, 2.8151800252334224e-06)

  test "preferred over the 4PL for asymmetric data":
    check fit.aic < fitLL4(asymDose, asymResponse).aic

  test "summary includes Asymmetry":
    check "5PL fit to 28 observations" in fit.summary
    check "Asymmetry" in fit.summary

  test "Asymmetry fixed at 1 reproduces the 4PL":
    let fixed5 = fitLL5(ryegrassDose, ryegrassRootLength, fixed = [(Asymmetry, 1.0)])
    let four = fitLL4(ryegrassDose, ryegrassRootLength)
    for param in DrcParam:
      check close(fixed5.params[param], four.params[param], 1e-8)
      check close(fixed5.se[param] + 1.0, four.se[param] + 1.0, 1e-8)
    check close(fixed5.aic, four.aic, 1e-10)

  test "Asymmetry can only be fixed to a positive value, and only in a 5PL":
    expect ValueError: discard fitLL5(asymDose, asymResponse, fixed = [(Asymmetry, 0.0)])
    expect ValueError: discard fitLL4(asymDose, asymResponse, fixed = [(Asymmetry, 1.0)])

suite "3PL":
  test "is the 4PL with Bottom fixed":
    let three = fitLL3(ryegrassDose, ryegrassRootLength)
    let four = fitLL4(ryegrassDose, ryegrassRootLength, fixed = [(Bottom, 0.0)])
    check three.model == LL3
    check three.fixed == {Bottom, Asymmetry}
    check three.params == four.params
    check three.aic == four.aic
    # R nls with three free parameters (see tinference).
    check close(three.params[Top], 7.85542760957709, 1e-6)
    check "3PL fit" in three.summary

  test "Bottom can be fixed at another value":
    let fit = fitLL3(ryegrassDose, ryegrassRootLength, bottom = 0.5)
    check fit.params[Bottom] == 0.5
