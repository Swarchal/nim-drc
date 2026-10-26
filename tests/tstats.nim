import std/[math, unittest]
import drc/stats

# Reference values from R: pbeta, pt and qt.

proc close(a, b: float, rel = 1e-10): bool =
  abs(a - b) <= rel * abs(b)

suite "incomplete beta":
  test "matches R pbeta":
    check close(incompleteBeta(2.0, 3.0, 0.3), 0.34829999999999989)
    check close(incompleteBeta(0.5, 0.5, 0.9), 0.79516723530086664)
    check close(incompleteBeta(10.0, 0.5, 0.01), 1.7700349655285852e-21)
    check close(incompleteBeta(20.0, 50.0, 0.01), 7.2367596828666965e-24)
    check close(incompleteBeta(0.3, 4.0, 0.2), 0.88935805582830962)
    check close(incompleteBeta(30.0, 1.5, 0.7), 8.0058741590440571e-05)
    check close(incompleteBeta(200.0, 200.0, 0.5), 0.5)

  test "edge cases":
    check incompleteBeta(2.0, 3.0, 0.0) == 0.0
    check incompleteBeta(2.0, 3.0, 1.0) == 1.0
    check incompleteBeta(0.0, 3.0, 0.5).isNaN
    check incompleteBeta(2.0, 3.0, NaN).isNaN

suite "Student's t":
  test "CDF matches R pt":
    let ts = [-3.0, -0.5, 0.0, 0.5, 2.2, 10.0]
    let expected = [
      0.019970984035859413, 0.32166498159093165, 0.5, 0.6783350184090684,
      0.95367366491018268, 0.99971899818864207,
    ]
    for i in 0 ..< ts.len:
      check close(studentTCdf(ts[i], 4.0), expected[i])
    check close(studentTCdf(1.5, 1.0), 0.81283295818900125)
    check close(studentTCdf(40.0, 3.0), 0.9999828096596054)

  test "quantile matches R qt":
    let ps = [0.975, 0.995, 0.9, 0.6]
    let expected = [
      (
        1.0,
        [
          12.706204736174694, 63.656741162871526, 3.0776835371752544,
          0.32491969623290634,
        ],
      ),
      (
        2.0,
        [
          4.3026527297494619, 9.9248432009182892, 1.8856180831641269,
          0.28867513459481314,
        ],
      ),
      (
        3.0,
        [
          3.1824463052837078, 5.8409093097333544, 1.6377443536962091,
          0.27667066233268983,
        ],
      ),
      (
        5.0,
        [
          2.5705818356363137, 4.0321429835552287, 1.4758840488244804,
          0.26718086570414518,
        ],
      ),
      (
        10.0,
        [2.2281388519862739, 3.16927267261695, 1.3721836411103356, 0.26018482949208016],
      ),
      (
        20.0,
        [
          2.0859634472658648, 2.8453397097861082, 1.3253407069850467,
          0.25674275385450185,
        ],
      ),
      (
        30.0,
        [
          2.0422724563012378, 2.7499956535672254, 1.3104150253913955,
          0.25560536495191277,
        ],
      ),
      (
        100.0,
        [
          1.9839715185235516, 2.6258905214380182, 1.2900747613465124,
          0.25402218245822766,
        ],
      ),
      (
        1000.0,
        [
          1.9623390808264076, 2.5807546980659501, 1.2823987214609247,
          0.25341451583949876,
        ],
      ),
    ]
    for (df, qs) in expected:
      for i in 0 ..< ps.len:
        check close(studentTQuantile(ps[i], df), qs[i])
    # Just above the median, where 1 − p carries only ~9 significant digits.
    check close(studentTQuantile(0.5000001, 1.0), 3.1415926533014113e-07, 1e-8)
    check close(studentTQuantile(0.025, 5.0), -2.570581835636315)

  test "quantile edge cases":
    check studentTQuantile(0.5, 7.0) == 0.0
    check studentTQuantile(0.0, 7.0) == NegInf
    check studentTQuantile(1.0, 7.0) == Inf
    check studentTQuantile(1.5, 7.0).isNaN
    check studentTQuantile(0.975, 0.0).isNaN

suite "normal quantile":
  test "matches R's qnorm":
    for (p, z) in [
      (1e-300, -37.047096299361201),
      (1e-20, -9.2623400897984052),
      (1e-8, -5.61200124417479),
      (0.001, -3.0902323061678132),
      (0.025, -1.9599639845400538),
      (0.1, -1.2815515655446008),
      (0.3, -0.52440051270804067),
      (0.7, 0.52440051270804067),
      (0.9, 1.2815515655446008),
      (0.975, 1.9599639845400536),
      (0.999, 3.0902323061678132),
      (1.0 - 1e-12, 7.0344869100478356),
    ]:
      check close(normalQuantile(p), z, 1e-14)

  test "edge cases":
    check normalQuantile(0.5) == 0.0
    check normalQuantile(0.0) == NegInf
    check normalQuantile(1.0) == Inf
    check normalQuantile(-0.1).isNaN
    check normalQuantile(NaN).isNaN

suite "log-likelihood":
  test "matches R logLik.nls":
    # Ryegrass LL4 fit: SSR 5.400214650161578 over 24 observations.
    check close(gaussianLogLik(5.400214650161578, 24), -16.15514326645846, 1e-12)
