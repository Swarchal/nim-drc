import std/[math, unittest]
import drc
import ryegrass

# References from R: delta-method variance gᵀ·vcov(nls)·g with finite-difference
# gradients, f ± qt(0.975, df)·sqrt(v) (confidence) and f ± qt·sqrt(v + σ²/w)
# (prediction). The 5PL data are those in tmodelfits.

proc close(a, b: float; rel = 1e-5): bool = abs(a - b) <= rel * max(1e-3, abs(b))

proc checkBand(b: Band; fit, lower, upper: float) =
  check close(b.fit, fit, 1e-7)
  check close(b.lower, lower)
  check close(b.upper, upper)

suite "bands: ryegrass 4PL":
  let fit = fitLL4(ryegrassDose, ryegrassRootLength)

  test "confidence band matches R":
    let b = fit.confidenceBand([0.0, 1.0, 3.0, 10.0, 50.0])
    checkBand(b[0], 7.7929623274932949, 7.3971856820392237, 8.1887389729473661)
    checkBand(b[1], 7.5411247294460635, 7.2023858070231102, 7.8798636518690168)
    checkBand(b[2], 4.2414615880532471, 3.6701670256148753, 4.8127561504916194)
    checkBand(b[3], 0.68887559118097264, 0.33418128688500359, 1.0435698954769417)
    checkBand(b[4], 0.48316729343851883, 0.048596532175530782, 0.91773805470150682)

  test "prediction band matches R":
    let b = fit.predictionBand([0.0, 1.0, 3.0, 10.0, 50.0])
    checkBand(b[0], 7.7929623274932949, 6.6390465598820541, 8.9468780951045357)
    checkBand(b[1], 7.5411247294460635, 6.4055074949787345, 8.6767419639133934)
    checkBand(b[2], 4.2414615880532471, 3.0162024633263136, 5.4667207127801802)
    checkBand(b[3], 0.68887559118097264, -0.45160259587696572, 1.8293537782389109)
    checkBand(b[4], 0.48316729343851883, -0.68462300974080337, 1.6509575966178409)

  test "a heavier prediction weight narrows only the noise term":
    checkBand(fit.predictionBand(3.0, weight = 4.0),
              4.2414615880532471, 3.4539985168648957, 5.028924659241599)

  test "at zero dose the confidence band is the plateau's interval":
    let b = fit.confidenceBand(0.0)
    let ci = fit.confint(Top)
    check close(b.lower, ci.lower, 1e-12)
    check close(b.upper, ci.upper, 1e-12)

  test "single-dose and sequence forms agree":
    check fit.confidenceBand(3.0) == fit.confidenceBand([3.0])[0]
    check fit.predictionBand(3.0) == fit.predictionBand([3.0])[0]

  test "invalid arguments are rejected":
    expect ValueError: discard fit.confidenceBand(-1.0)
    expect ValueError: discard fit.confidenceBand(Inf)
    expect ValueError: discard fit.predictionBand(1.0, weight = 0.0)
    expect ValueError: discard fit.confidenceBand(1.0, level = 1.0)

suite "bands: 5PL":
  const
    dose = [0.0, 0.0, 1e-09, 1e-09, 3.1622776601683795e-09, 3.1622776601683795e-09,
            1e-08, 1e-08, 3.162277660168379e-08, 3.162277660168379e-08, 1e-07, 1e-07,
            3.162277660168379e-07, 3.162277660168379e-07, 1e-06, 1e-06,
            3.162277660168379e-06, 3.162277660168379e-06, 1e-05, 1e-05,
            3.1622776601683795e-05, 3.1622776601683795e-05, 0.0001, 0.0001,
            0.00031622776601683794, 0.00031622776601683794, 0.001, 0.001]
    response = [5.0, 5.6, 8.34, 7.1, 10.65, 9.57, 16.18, 18.74, 22.64, 22.38, 37.01, 36.75,
                54.87, 52.8, 75.55, 77.0, 86.17, 87.94, 89.67, 90.89, 90.97, 94.18, 92.39,
                95.46, 95.3, 94.61, 89.96, 93.92]
  let fit = fitLL5(dose, response)

  test "confidence and prediction bands match R":
    let c = fit.confidenceBand([1e-9, 1e-7, 1e-6, 1e-4])
    let p = fit.predictionBand([1e-9, 1e-7, 1e-6, 1e-4])
    checkBand(c[0], 8.4489523995020441, 7.2367245722316538, 9.6611802267724336)
    checkBand(c[1], 36.147725068744329, 34.669213687630176, 37.626236449858482)
    checkBand(c[2], 75.109751551167946, 73.378064848713279, 76.841438253622613)
    checkBand(c[3], 93.125411452625599, 91.991153023947874, 94.259669881303324)
    checkBand(p[0], 8.4489523995020441, 4.8649458152511915, 12.032958983752897)
    checkBand(p[1], 36.147725068744329, 32.465116917606217, 39.830333219882441)
    checkBand(p[2], 75.109751551167946, 71.318400113197171, 78.901102989138721)
    checkBand(p[3], 93.125411452625599, 89.567020092225334, 96.683802813025864)
