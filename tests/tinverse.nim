import std/[math, unittest]
import drc
import ryegrass
import fixtures

# References from R: nls() fits in x = log10(dose), the dose x₀ where the
# fitted curve equals y₀ found by uniroot (and checked against the closed
# form), and the delta method on x₀ with finite-difference gradients:
# x₀ ± qt(0.5 + level/2, df)·sqrt(gᵀ·vcov·g) for doseAtCI, plus
# (∂x₀/∂y₀)²·σ²/(m·w) for calibrationCI, back-transformed to dose.

proc close(a, b: float, rel = 1e-5): bool =
  abs(a - b) <= rel * abs(b)

proc checkInterval(ci: ConfInt, lower, upper: float) =
  check close(ci.lower, lower)
  check close(ci.upper, upper)

suite "inverse prediction: ryegrass 4PL":
  let fit = fitLL4(ryegrassDose, ryegrassRootLength)

  test "dose at a response matches R":
    check close(fit.doseAt(2.0), 4.7907792055984002)
    check close(fit.doseAt(4.0), 3.1359262652104349)
    check close(fit.doseAt(7.0), 1.5088523801581555)

  test "the fitted curve passes through the inverse":
    for response in [0.6, 2.0, 4.0, 7.0, 7.7]:
      check close(fit.predict(fit.doseAt(response)), response, 1e-10)

  test "the relative ECx is the dose at the matching absolute response":
    let (b, t) = (fit.params[Bottom], fit.params[Top])
    for pct in [10.0, 50.0, 90.0]:
      # A falling curve: pct% of the way from Top down to Bottom.
      let response = t - pct / 100.0 * (t - b)
      check close(fit.doseAt(response), fit.ecx(pct), 1e-10)
      # The intervals differ: a fixed response level carries the plateaus'
      # uncertainty, a fixed fraction of the way between them does not.
      check fit.doseAtCI(response) != fit.ecxCI(pct)

  test "confidence intervals match R":
    checkInterval(fit.doseAtCI(2.0), 4.1051221325543876, 5.5909579924026671)
    checkInterval(fit.doseAtCI(4.0), 2.8250043305052461, 3.4810684835579959)
    checkInterval(fit.doseAtCI(7.0), 1.2061857008218901, 1.8874668333057179)
    checkInterval(fit.doseAtCI(4.0, 0.9), 2.8765519478541299, 3.4186879705658462)

  test "calibration intervals match R":
    checkInterval(fit.calibrationCI(2.0), 3.4123922117055625, 6.7259458974449231)
    checkInterval(fit.calibrationCI(4.0), 2.5044949781176449, 3.9265535074970823)
    checkInterval(fit.calibrationCI(7.0), 0.86122555796931233, 2.6434834452393869)
    checkInterval(fit.calibrationCI(4.0, 0.9), 2.6039362799117978, 3.7766029901352929)

  test "replicates and weight shrink only the measurement term":
    checkInterval(
      fit.calibrationCI(4.0, weight = 2.0, replicates = 3),
      2.7472397863390343,
      3.5796050966273607,
    )
    check fit.calibrationCI(4.0, weight = 3.0) == fit.calibrationCI(4.0, replicates = 3)

  test "responses outside the plateaus give NaN":
    for response in [fit.params[Bottom], fit.params[Top], -1.0, 0.4, 8.0, 100.0]:
      check fit.doseAt(response).isNaN
      check fit.doseAtCI(response).lower.isNaN
      check fit.calibrationCI(response).upper.isNaN

  test "invalid arguments are rejected":
    expect ValueError:
      discard fit.doseAt(NaN)
    expect ValueError:
      discard fit.doseAtCI(Inf)
    expect ValueError:
      discard fit.doseAtCI(4.0, level = 1.0)
    expect ValueError:
      discard fit.calibrationCI(4.0, weight = 0.0)
    expect ValueError:
      discard fit.calibrationCI(4.0, replicates = 0)

suite "inverse prediction: 5PL":
  let fit = fitLL5(asymDose, asymResponse)

  test "dose at a response and its intervals match R":
    check close(fit.doseAt(20.0), 2.0750018220892061e-08)
    check close(fit.doseAt(50.0), 2.3890520327457465e-07)
    check close(fit.doseAt(80.0), 1.4352609794322393e-06)
    checkInterval(fit.doseAtCI(20.0), 1.7140864111947981e-08, 2.5119110294283813e-08)
    checkInterval(fit.doseAtCI(50.0), 2.1591311627430226e-07, 2.6434566429561058e-07)
    checkInterval(fit.doseAtCI(80.0), 1.2371946258986813e-06, 1.6650363944027428e-06)
    checkInterval(
      fit.calibrationCI(20.0), 1.2465669750847073e-08, 3.4539921622590294e-08
    )
    checkInterval(
      fit.calibrationCI(50.0), 1.9247627371151677e-07, 2.9653367166288158e-07
    )
    checkInterval(
      fit.calibrationCI(80.0), 1.0476908650060367e-06, 1.9662041045560911e-06
    )

  test "the absolute and relative EC50 differ":
    # Plateaus at 4.9 and 93.2: a response of 50 is past the half-way point.
    check fit.doseAt(50.0) > fit.ec50
    let half = 0.5 * (fit.params[Bottom] + fit.params[Top])
    check close(fit.doseAt(half), fit.ec50, 1e-10)

  test "the fitted curve passes through the inverse":
    for response in [6.0, 20.0, 50.0, 80.0, 93.0]:
      check close(fit.predict(fit.doseAt(response)), response, 1e-10)

suite "inverse prediction: special cases":
  test "a flat curve gives NaN":
    let fit = fitLL4(
      [0.1, 1.0, 10.0, 100.0, 0.1, 1.0, 10.0, 100.0],
      [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0],
    )
    check fit.doseAt(5.0).isNaN

  test "fixed plateaus: absolute IC50 of normalised data":
    let doses = [1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4]
    let viability = [98.2, 95.1, 71.4, 22.9, 6.3, 2.0]
    let fit = fitLL4(doses, viability, fixed = [(Bottom, 0.0), (Top, 100.0)])
    # With the plateaus at 0 and 100, absolute and relative IC50 coincide, and
    # the fixed plateaus add nothing to the interval.
    check close(fit.doseAt(50.0), fit.ec50, 1e-10)
    check close(fit.doseAtCI(50.0).lower, fit.ec50CI.lower, 1e-10)
    check close(fit.doseAtCI(50.0).upper, fit.ec50CI.upper, 1e-10)

  test "curve views of a joint fit support inverse prediction":
    let (doses, resp, curve) = ryegrassCurves(4.0)
    let fit = fitLL4Multi(doses, resp, curve, shared = {Bottom, Top, Hill})
    check close(fit.curves[1].doseAt(4.0) / fit.curves[0].doseAt(4.0), 4.0, 1e-6)
