import std/[math, strutils, unittest]
import drc
import ryegrass

proc close(a, b: float; rel: float): bool = abs(a - b) <= rel * abs(b)

proc rosenbrock(p: array[3, float]): float =
  (1.0 - p[0]) ^ 2 + 100.0 * (p[1] - p[0] ^ 2) ^ 2 + p[2] ^ 2

suite "Nelder–Mead":
  test "minimises the Rosenbrock function":
    let r = nelderMead(rosenbrock, [-1.2, 1.0, 0.5], [true, true, true])
    check r.converged
    check abs(r.params[0] - 1.0) < 1e-6
    check abs(r.params[1] - 1.0) < 1e-6
    check abs(r.params[2]) < 1e-6

  test "fixed parameters do not move":
    let r = nelderMead(rosenbrock, [-1.2, 1.0, 0.5], [true, true, false])
    check r.params[2] == 0.5
    check abs(r.params[0] - 1.0) < 1e-6

  test "NaN counts as +Inf":
    proc f(p: array[1, float]): float = (if p[0] < 0.0: NaN else: (p[0] - 2.0) ^ 2)
    let r = nelderMead(f, [1.0], [true])
    check abs(r.params[0] - 2.0) < 1e-6

  test "returns the start when nothing is free":
    let r = nelderMead(rosenbrock, [-1.2, 1.0, 0.5], [false, false, false])
    check r.converged
    check r.params == [-1.2, 1.0, 0.5]
    check r.value == rosenbrock([-1.2, 1.0, 0.5])

  test "stops at maxEval":
    let r = nelderMead(rosenbrock, [-1.2, 1.0, 0.5], [true, true, true], NMOptions(maxEval: 20))
    check not r.converged
    check r.evaluations <= 25

suite "multi-start":
  # Five sparse, noisy points (found by a random search) where Levenberg–Marquardt
  # from the default start stops in a local minimum. 669.8355816 is also the
  # best SSR SciPy's curve_fit finds from 3000 random starts.
  const
    dose = [1e-9, 1.161749920676094e-8, 3.1629918890385245e-8, 0.00000021202993733263668,
            0.00002911371875035665]
    resp = [17.840260081080757, -16.695409693569324, -9.847221967902513,
            -4.111836129419825, 44.44746199771957]

  test "escapes a local minimum":
    let single = fitLL4(dose, resp, opts = FitOptions(multiStart: false, fallback: false))
    let multi = fitLL4(dose, resp)
    check single.ssr > 2000.0
    check close(multi.ssr, 669.8355816295317, 1e-8)

  test "gives the same answer as a single start on well-behaved data":
    let single = fitLL4(ryegrassDose, ryegrassRootLength,
                        opts = FitOptions(multiStart: false, fallback: false))
    let multi = fitLL4(ryegrassDose, ryegrassRootLength)
    check multi.params == single.params

suite "Nelder–Mead fallback":
  # Two iterations are not enough for Levenberg–Marquardt to converge.
  let starved = LMOptions(maxIter: 2)

  test "recovers when Levenberg–Marquardt does not converge":
    let fit = fitLL4(ryegrassDose, ryegrassRootLength, opts = FitOptions(lm: starved))
    check fit.usedFallback
    check fit.converged
    check fwNotConverged notin fit.warnings
    check close(fit.params[LogEC50], 0.48543111513116, 1e-6)
    check close(fit.params[Hill], -2.9822286375650706, 1e-6)
    check close(fit.ssr, 5.400214650161578, 1e-10)
    check "the Nelder–Mead fallback converged" in fit.summary

  test "without the fallback, non-convergence is reported":
    let fit = fitLL4(ryegrassDose, ryegrassRootLength,
                     opts = FitOptions(lm: starved, fallback: false))
    check not fit.usedFallback
    check not fit.converged
    check fwNotConverged in fit.warnings

  test "every parameter fixed, with Hill = 0 and zero-dose controls":
    # Previously the model gave NaN at zero dose and Nelder–Mead indexed out of bounds.
    let fit = fitLL4([0.0, 0.0, 1.0, 10.0, 100.0, 1000.0], [5.0, 6.0, 5.0, 6.0, 5.0, 6.0],
                     fixed = [(Bottom, 0.0), (Top, 10.0), (LogEC50, 1.0), (Hill, 0.0)])
    check fit.converged
    check abs(fit.ssr - 3.0) < 1e-12  # the flat curve is 5; residuals 0, 1, 0, 1, 0, 1
    check fit.predict(0.0) == 5.0

  test "is not used when Levenberg–Marquardt converges":
    check not fitLL4(ryegrassDose, ryegrassRootLength).usedFallback
