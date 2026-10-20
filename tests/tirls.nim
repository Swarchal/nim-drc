import std/[math, strutils, unittest]
import drc
import drc/linalg
import ryegrass

# Reference values from R's robustbase::nlrob (Huber, then bisquare started from
# the Huber fit) on a 4PL inhibition curve (truth: Bottom 2, Top 98, LogEC50
# −6.5, Hill −1.2, noise sd 3) with gross errors added at observations 4, 16
# and 25. Standard errors follow MASS::rlm's default (Huber) covariance,
# computed independently in R from finite-difference Jacobians at the nlrob fit.
const
  dose = [1e-09, 1e-09, 1e-09, 3.1622776601683795e-09, 3.1622776601683795e-09,
          3.1622776601683795e-09, 1e-08, 1e-08, 1e-08, 3.162277660168379e-08,
          3.162277660168379e-08, 3.162277660168379e-08, 1e-07, 1e-07, 1e-07,
          3.162277660168379e-07, 3.162277660168379e-07, 3.162277660168379e-07, 1e-06, 1e-06,
          1e-06, 3.162277660168379e-06, 3.162277660168379e-06, 3.162277660168379e-06, 1e-05,
          1e-05, 1e-05, 3.1622776601683795e-05, 3.1622776601683795e-05,
          3.1622776601683795e-05, 0.0001, 0.0001, 0.0001]
  response = [104.03, 90.24, 99.16, 95.92, 136.26, 96.97, 90.44, 95.81, 93.91, 102.27, 92.98,
              91.24, 77.88, 76.72, 75.56, 48.83, 16.45, 49.28, 24.15, 20.67, 21.35, 12.34, 9.33,
              6.18, 2.95, 35.12, 9.3, 1.57, 1.65, 5.39, -0.56, 1.22, 4.74]

proc close(a, b: float; rel: float): bool = abs(a - b) <= rel * max(1e-3, abs(b))

suite "loss functions":
  test "match robustbase's weights and ψ'":
    check close(psiWeight(Huber, 3.0, 1.345), 0.4483333333333333, 1e-12)
    check psiWeight(Huber, 0.5, 1.345) == 1.0
    check psiDerivative(Huber, 0.5, 1.345) == 1.0
    check psiDerivative(Huber, 3.0, 1.345) == 0.0
    check close(psiWeight(Tukey, 0.5, 4.685), 0.9773499, 1e-6)
    check close(psiWeight(Tukey, 3.0, 4.685), 0.3480560, 1e-6)
    check close(psiDerivative(Tukey, 0.5, 4.685), 0.9323091, 1e-6)
    check close(psiDerivative(Tukey, 3.0, 4.685), -0.6195708, 1e-6)
    check psiWeight(Tukey, 5.0, 4.685) == 0.0

  test "ψ is u times the weight":
    for loss in [Squared, Huber, Tukey]:
      for u in [-6.0, -1.0, 0.3, 2.0, 4.0]:
        check abs(psi(loss, u, defaultTuning(loss)) -
                  u * psiWeight(loss, u, defaultTuning(loss))) < 1e-14

  test "ψ' matches finite differences":
    for loss in [Huber, Tukey]:
      for u in [-3.7, -0.4, 0.9, 2.2]:
        let c = defaultTuning(loss)
        let numeric = (psi(loss, u + 1e-6, c) - psi(loss, u - 1e-6, c)) / 2e-6
        check abs(numeric - psiDerivative(loss, u, c)) < 1e-6

  test "MAD scale":
    check close(madScale([1.0, -2.0, 3.0]), 2.0 / 0.6745, 1e-14)
    check close(madScale([1.0, -2.0, 3.0, -4.0]), 2.5 / 0.6745, 1e-14)

suite "robust fits with outliers":
  let
    ls = fitLL4(dose, response)
    huber = fitLL4(dose, response, opts = FitOptions(loss: Huber))
    tukey = fitLL4(dose, response, opts = FitOptions(loss: Tukey))

  test "Huber matches nlrob":
    check huber.converged
    check fwRobustNotConverged notin huber.warnings
    check huber.tuning == 1.345
    check close(huber.params[Bottom], 3.09319936475949, 2e-5)
    check close(huber.params[Top], 98.355673985754407, 2e-5)
    check close(huber.params[LogEC50], -6.537473458202502, 2e-5)
    check close(huber.params[Hill], -1.1389656628274041, 2e-5)
    check close(huber.sigma, 2.7941781440171432, 2e-5)
    check close(huber.se[Bottom], 1.2976582257306499, 2e-5)
    check close(huber.se[Top], 1.3330832516095361, 2e-5)
    check close(huber.se[LogEC50], 0.031369300343063254, 2e-5)
    check close(huber.se[Hill], 0.084179940577245621, 2e-5)
    let ci = huber.ec50CI
    check close(ci.lower, 2.5024713776740475e-07, 2e-5)
    check close(ci.upper, 3.3626677960875993e-07, 2e-5)

  test "Tukey matches nlrob, rejecting the gross errors":
    check tukey.converged
    check close(tukey.params[Bottom], 2.4525384825582575, 1e-5)
    check close(tukey.params[Top], 97.13527747669437, 1e-5)
    check close(tukey.params[LogEC50], -6.5079575778279546, 1e-5)
    check close(tukey.params[Hill], -1.1397603765314477, 1e-5)
    check close(tukey.sigma, 2.280359795815404, 1e-5)
    check close(tukey.se[LogEC50], 0.025984811223164244, 1e-5)
    check close(tukey.se[Hill], 0.069838915856252015, 1e-5)
    let ci = tukey.ec50CI
    check close(ci.lower, 2.7472457867635697e-07, 1e-5)
    check close(ci.upper, 3.5090319969724354e-07, 1e-5)
    var rejected: seq[int]
    for i, rw in tukey.robustWeights:
      if rw == 0.0: rejected.add i
    check rejected == @[4, 9, 16, 25]  # the three planted errors and one large natural residual
    check close(tukey.robustWeights[0], 0.32131392401009551, 1e-5)
    check close(tukey.robustWeights[26], 0.60003789842593058, 1e-5)

  test "robust fits are closer to the truth than least squares":
    check abs(ls.params[LogEC50] - -6.5) > 0.15
    check abs(huber.params[LogEC50] - -6.5) < 0.05
    check abs(tukey.params[LogEC50] - -6.5) < 0.02

  test "no likelihood-based criteria for robust fits":
    check huber.loss == Huber
    check huber.aic.isNaN and huber.bic.isNaN and huber.logLik.isNaN
    check ls.robustWeights.len == 0
    check ls.aic.isFinite

  test "summary describes the robust fit":
    let text = tukey.summary
    check "robust (Tukey bisquare, tuning 4.685)" in text
    check "Robust residual scale (MAD)" in text
    check "0 (rejected): 4" in text
    check "AIC" notin text

  test "bands use the robust covariance and scale":
    let b = huber.confidenceBand(3e-7)
    check b.fit == huber.predict(3e-7)
    check b.upper - b.lower < ls.confidenceBand(3e-7).upper - ls.confidenceBand(3e-7).lower
    let p = huber.predictionBand(3e-7)
    check p.upper - p.lower > b.upper - b.lower

suite "robust fits: options and edge cases":
  test "a huge Huber tuning constant is least squares":
    let ls = fitLL4(ryegrassDose, ryegrassRootLength)
    let wide = fitLL4(ryegrassDose, ryegrassRootLength, opts = FitOptions(loss: Huber, tuning: 1e9))
    for param in Bottom .. Hill:
      check close(wide.params[param], ls.params[param], 1e-8)
    for rw in wide.robustWeights: check rw == 1.0

  test "on clean data a Huber fit stays close to least squares":
    let ls = fitLL4(ryegrassDose, ryegrassRootLength)
    let h = fitLL4(ryegrassDose, ryegrassRootLength, opts = FitOptions(loss: Huber))
    check h.converged
    check abs(h.params[LogEC50] - ls.params[LogEC50]) < 0.02

  test "user-excluded observations get weight 0 but are not counted as rejected":
    var d = @dose
    var r = @response
    var w = newSeq[float](d.len)
    for x in w.mitems: x = 1.0
    d.add 1e-6
    r.add 500.0
    w.add 0.0
    let fit = fitLL4(d, r, w, opts = FitOptions(loss: Tukey))
    check fit.robustWeights[^1] == 0.0
    check fit.n == dose.len
    check "0 (rejected): 4" in fit.summary
    check close(fit.params[LogEC50], -6.5079575778279546, 1e-5)

  test "5PL robust fits run":
    let fit = fitLL5(dose, response, opts = FitOptions(loss: Huber))
    check fit.converged
    check fit.se[Asymmetry].isFinite

  test "invalid tuning is rejected":
    expect ValueError: discard fitLL4(dose, response, opts = FitOptions(loss: Huber, tuning: 0.0))
    expect ValueError: discard fitLL4(dose, response, opts = FitOptions(loss: Tukey, tuning: -1.0))
