import std/[math, unittest]
import drc
import ryegrass
import fixtures

# References from R: the exact ends of profile-likelihood intervals, found by
# uniroot on τ(β) = sqrt((S(β) − Ŝ)/σ̂²) − qt(0.5 + level/2, df), where S(β)
# is the SSR of nls() refitted with the parameter held at β. ECx intervals
# profile log10(ECx) in nls() with LogEC50 rewritten as
# log10(ECx) + log10(r^(−1/A) − 1)/Hill. R's own confint(nls fit), which
# interpolates a spline through a coarse profile, agrees to ~1e-4.

proc close(a, b: float, rel = 1e-5): bool =
  abs(a - b) <= rel * abs(b)

proc checkInterval(ci: ConfInt, lower, upper: float, rel = 1e-5) =
  check close(ci.lower, lower, rel)
  check close(ci.upper, upper, rel)

suite "profile likelihood: ryegrass 4PL":
  let fit = fitLL4(ryegrassDose, ryegrassRootLength)

  test "parameter intervals match R":
    checkInterval(fit.profileCI(Bottom), 0.010240384807400379, 0.90674598455888344)
    checkInterval(fit.profileCI(Top), 7.4053722219950346, 8.1924590594593916)
    checkInterval(fit.profileCI(LogEC50), 0.42992956339078325, 0.54100316108586066)
    checkInterval(fit.profileCI(Hill), -4.2328922557496984, -2.1820088143027823)
    checkInterval(fit.profileCI(Hill, 0.9), -3.9614478750706144, -2.2996989843439044)

  test "ECx intervals match R":
    checkInterval(fit.ecxProfileCI(10.0), 1.1023462816423895, 1.9030591389464846)
    checkInterval(fit.ec50ProfileCI, 2.6910983092900738, 3.4753869105315931)
    checkInterval(fit.ecxProfileCI(90.0), 5.0280824852950969, 8.799249697770561)

  test "for a 4PL, the EC50 interval is the LogEC50 one back-transformed":
    let ci = fit.profileCI(LogEC50)
    checkInterval(fit.ec50ProfileCI, pow(10.0, ci.lower), pow(10.0, ci.upper), 1e-8)

  test "the profile interval is asymmetric where the Wald one is not":
    # Hill's SSR rises more slowly towards steeper slopes.
    let wald = fit.confint(Hill)
    let prof = fit.profileCI(Hill)
    check prof.lower < wald.lower
    check fit.params[Hill] - prof.lower > prof.upper - fit.params[Hill]

  test "fixed parameters, robust fits and joint-fit curves":
    let unit = fitLL4(ryegrassDose, ryegrassRootLength, fixed = [(Hill, -3.0)])
    check unit.profileCI(Hill) == (-3.0, -3.0)
    check unit.profileCI(Asymmetry) == (1.0, 1.0)
    expect ValueError:
      discard
        fitLL4(ryegrassDose, ryegrassRootLength, fixed = [(LogEC50, 0.5)]).ec50ProfileCI
    expect ValueError:
      discard fitLL4(ryegrassDose, ryegrassRootLength, opts = FitOptions(loss: Huber))
        .profileCI(Hill)
    let (doses, resp, curve) = ryegrassCurves(2.0)
    let joint = fitLL4Multi(doses, resp, curve, shared = {Bottom, Top, Hill})
    expect ValueError:
      discard joint.curves[0].profileCI(Hill)

  test "invalid arguments are rejected":
    expect ValueError:
      discard fit.profileCI(Hill, 1.0)
    expect ValueError:
      discard fit.ecxProfileCI(0.0)
    expect ValueError:
      discard fitLL4(ryegrassDose, ryegrassRootLength, bounds = [(LogEC50, 0.0 .. 1.0)]).ec50ProfileCI

suite "profile likelihood: other fits":
  test "weighted fit matches R":
    var weights: seq[float]
    for i in 0 ..< 8:
      weights.add [1.0, 2.0, 0.5]
    let fit = fitLL4(ryegrassDose, ryegrassRootLength, weights)
    checkInterval(fit.profileCI(Hill), -4.5440144460412171, -2.4174446327785288)

  test "3PL matches R":
    let fit = fitLL3(ryegrassDose, ryegrassRootLength)
    checkInterval(fit.profileCI(Top), 7.4369064159013858, 8.2868091500773033)
    checkInterval(fit.profileCI(LogEC50), 0.45927329186872706, 0.56864741086487192)
    check fit.profileCI(Bottom) == (0.0, 0.0)

  test "5PL matches R":
    let fit = fitLL5(asymDose, asymResponse)
    checkInterval(fit.profileCI(Bottom), 2.7763340271183843, 6.917298531200581)
    checkInterval(fit.profileCI(LogEC50), -6.3048032858028424, -5.919346216116967)
    checkInterval(fit.profileCI(Hill), 0.96172600996518587, 1.6275520550338431)
    checkInterval(fit.profileCI(Asymmetry), 0.25776754890653869, 0.59771321728661286)
    checkInterval(
      fit.ecxProfileCI(10.0), 3.8571099990814115e-09, 1.1074662446952295e-08
    )
    checkInterval(fit.ec50ProfileCI, 2.0151654620179148e-07, 2.5378077017465028e-07)

  test "a bound the data do not exclude is an end of the interval":
    # R: refits by nls(algorithm = "port") within the same bounds.
    let doses = [1e-9, 1e-9, 1e-8, 1e-8, 1e-7, 1e-7, 1e-6, 1e-6, 1e-5, 1e-5]
    let viability = [101.3, 97.8, 96.1, 99.0, 78.2, 74.9, 31.5, 36.0, 4.1, 1.2]
    let fit =
      fitLL4(doses, viability, bounds = [(Bottom, 0.0 .. Inf), (Top, -Inf .. 100.0)])
    check fit.atBound == {Bottom, Top}
    let bottom = fit.profileCI(Bottom)
    check bottom.lower == 0.0
    check close(bottom.upper, 4.9726843127327065)
    let top = fit.profileCI(Top)
    check close(top.lower, 95.488117157774738)
    check top.upper == 100.0
    checkInterval(fit.profileCI(Hill), -1.0679594513682655, -0.75084447350589156)
    checkInterval(fit.profileCI(LogEC50), -6.4864503833554892, -6.2529438175975329)

suite "profile likelihood: a partial curve":
  # Only the upper half of the curve is observed, so the data do not bound
  # Bottom below or LogEC50 above. As LogEC50 → ∞ (and Bottom → −∞) the 4PL
  # tends to the power law Top − C·dose^k, whose fit in R has τ = 1.98, below
  # the cutoff qt(0.975, 6) = 2.45, so those ends are infinite. The finite
  # ends are from R as above, with optim() refits polished by nls().
  let doses = [1e-9, 1e-9, 1e-8, 1e-8, 1e-7, 1e-7, 1e-6, 1e-6, 1e-5, 1e-5]
  let resp = [100.4, 99.1, 98.7, 97.2, 93.0, 91.2, 80.1, 82.4, 68.0, 66.3]
  let fit = fitLL4(doses, resp)

  test "unbounded ends are infinite":
    let bottom = fit.profileCI(Bottom)
    check bottom.lower == -Inf
    check close(bottom.upper, 63.307179603036694)
    let logEC50 = fit.profileCI(LogEC50)
    check close(logEC50.lower, -6.1202157917817548)
    check logEC50.upper == Inf
    check fit.ec50ProfileCI.upper == Inf
    check close(fit.ec50ProfileCI.lower, pow(10.0, logEC50.lower), 1e-8)

  test "finite ends match R":
    check close(fit.profileCI(Top).lower, 98.39667350103646)
    check close(fit.profileCI(Hill).lower, -0.77776085647574922)

  test "ends approached at infinity agree with a full refit":
    # Top's and Hill's upper ends lie where the profile's minimum is at
    # Bottom → −∞, which R's refits stop short of. A multi-start fit with the
    # parameter fixed at the end must give τ = t quantile.
    let tq = studentTQuantile(0.975, float(fit.df))
    for param in [Top, Hill]:
      let refit = fitLL4(doses, resp, fixed = [(param, fit.profileCI(param).upper)])
      check close(sqrt(refit.ssr - fit.ssr) / fit.sigma, tq, 1e-4)

suite "profile likelihood: relative potency":
  # R: indexed nls() fits of the README's potency data, reparametrised with
  # log10 ECx(b) and d = log10(ECx(a)/ECx(b)) in place of LogEC50 for curves
  # a and b (LogEC50 = log10 ECx + log10(1/r − 1)/Hill), with d profiled as
  # above.
  let (doses, curve) = potencyDesign()
  const response = potencyResponse
  let parallel = fitLL4Multi(doses, response, curve, shared = {Bottom, Top, Hill})
  let separate = fitLL4Multi(doses, response, curve)

  test "parallel curves match R":
    checkInterval(
      parallel.relativePotencyProfileCI(0, 1), 3.4861125628264502, 4.1554784812591503
    )
    checkInterval(
      parallel.relativePotencyProfileCI(0, 2), 0.18109733960896374, 0.21578108265642343
    )
    checkInterval(
      parallel.relativePotencyProfileCI(2, 1, 10.0),
      17.636929648581379,
      21.019096402654903,
    )

  test "separate curves match R":
    checkInterval(
      separate.relativePotencyProfileCI(0, 1), 3.1506373939275805, 4.1381293525018625
    )
    # At EC10 the ratio depends on each curve's own Hill.
    checkInterval(
      separate.relativePotencyProfileCI(2, 1, 10.0),
      13.486145663710765,
      24.959736158078126,
    )

  test "swapping the curves inverts the interval":
    let ab = separate.relativePotencyProfileCI(0, 1)
    let ba = separate.relativePotencyProfileCI(1, 0)
    check close(ab.lower, 1.0 / ba.upper, 1e-8)
    check close(ab.upper, 1.0 / ba.lower, 1e-8)

  test "unsuitable fits are rejected":
    expect ValueError:
      discard parallel.relativePotencyProfileCI(1, 1)
    expect ValueError:
      discard parallel.relativePotencyProfileCI(0, 3)
    let sharedL = fitLL4Multi(doses, response, curve, shared = {LogEC50})
    expect ValueError:
      discard sharedL.relativePotencyProfileCI(0, 1, 10.0)
