import std/[math, unittest]
import drc/[linalg, models]

type P = array[drcParamCount, float]

proc checkGradient(x: float; p: P) =
  ## Compares the analytic gradient against central finite differences.
  var grad, scratch: P
  discard ll5(x, p, grad)
  for j in 0 ..< drcParamCount:
    let step = 1e-6 * max(1.0, abs(p[j]))
    var hi = p
    var lo = p
    hi[j] += step
    lo[j] -= step
    let numeric = (ll5(x, hi, scratch) - ll5(x, lo, scratch)) / (2.0 * step)
    check abs(numeric - grad[j]) <= 1e-6 * max(1.0, abs(numeric))

suite "logistic model":
  test "4PL value at EC50 is the midpoint":
    var grad: P
    check abs(ll5(0.5, [1.0, 9.0, 0.5, 1.3, 1.0], grad) - 5.0) < 1e-12
    check abs(ll5(0.5, [1.0, 9.0, 0.5, -1.3, 1.0], grad) - 5.0) < 1e-12

  test "5PL value at LogEC50 is Bottom + (Top − Bottom) / 2^Asymmetry":
    var grad: P
    check abs(ll5(0.5, [1.0, 9.0, 0.5, 1.3, 0.4], grad) - (1.0 + 8.0 / pow(2.0, 0.4))) < 1e-12

  test "plateaus and direction":
    var grad: P
    # Small Asymmetry approaches Bottom slowly (as sᴬ), so test well out.
    for a in [1.0, 0.3, 3.0]:
      check abs(ll5(-50.0, [1.0, 9.0, 0.0, 1.0, a], grad) - 1.0) < 1e-8
      check abs(ll5(50.0, [1.0, 9.0, 0.0, 1.0, a], grad) - 9.0) < 1e-8
      check abs(ll5(-50.0, [1.0, 9.0, 0.0, -1.0, a], grad) - 9.0) < 1e-8
      check abs(ll5(50.0, [1.0, 9.0, 0.0, -1.0, a], grad) - 1.0) < 1e-8

  test "no overflow at extreme arguments":
    var grad: P
    for x in [-1e6, 1e6]:
      for h in [-50.0, 50.0]:
        for a in [0.2, 1.0, 5.0]:
          check ll5(x, [1.0, 9.0, 0.0, h, a], grad).isFinite
          for g in grad: check g.isFinite

  test "zero dose gives the low-dose plateau":
    var grad: P
    check ll5(NegInf, [1.0, 9.0, 0.3, 2.0, 0.5], grad) == 1.0
    check grad == [1.0, 0.0, 0.0, 0.0, 0.0]
    check ll5(NegInf, [1.0, 9.0, 0.3, -2.0, 0.5], grad) == 9.0
    check grad == [0.0, 1.0, 0.0, 0.0, 0.0]

  test "a flat curve (Hill = 0) is finite at zero dose":
    var grad: P
    let f = ll5(NegInf, [1.0, 9.0, 0.3, 0.0, 0.5], grad)
    check abs(f - (1.0 + 8.0 * pow(0.5, 0.5))) < 1e-12
    var flat: P
    check abs(f - ll5(2.0, [1.0, 9.0, 0.3, 0.0, 0.5], flat)) < 1e-12
    for g in grad: check g.isFinite

  test "non-positive asymmetry is rejected with NaN":
    var grad: P
    check ll5(0.0, [1.0, 9.0, 0.0, 1.0, 0.0], grad).isNaN
    check ll5(0.0, [1.0, 9.0, 0.0, 1.0, -1.0], grad).isNaN

  test "analytic gradient matches finite differences":
    for p in [[0.0, 100.0, -6.5, 1.0, 1.0], [0.48, 7.79, 0.485, -2.98, 1.0],
              [10.0, -3.0, 2.0, 0.4, 1.0], [5.0, 95.0, -6.0, 1.3, 0.35], [5.0, 95.0, -6.0, -0.8, 2.5]]:
      for x in [-9.0, -6.5, -6.0, -1.0, 0.0, 0.485, 1.5, 3.0]:
        checkGradient(x, p)

suite "initial guess":
  test "recovers the shape of clean data":
    var x, y: seq[float]
    var grad: P
    let truth: P = [5.0, 95.0, -6.0, -1.5, 1.0]
    for k in 0 .. 8:
      x.add -8.0 + 0.5 * float(k)
      y.add ll5(x[^1], truth, grad)
    let g = drcGuess(x, y)
    check abs(g[ord(LogEC50)] - truth[ord(LogEC50)]) < 0.2
    check g[ord(Hill)] < 0.0
    check abs(g[ord(Hill)] - truth[ord(Hill)]) < 0.5
    check g[ord(Asymmetry)] == 1.0

  test "fixed parameters are respected":
    let g = drcGuess([0.0, 1.0, 2.0, 3.0], [1.0, 2.0, 8.0, 9.0],
                     [0.0, 10.0, 0.0, 0.0, 0.5], {Bottom, Top, Asymmetry})
    check g[ord(Bottom)] == 0.0
    check g[ord(Top)] == 10.0
    check g[ord(Hill)] > 0.0
    check g[ord(Asymmetry)] == 0.5
