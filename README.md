# drc

Dose-response curve fitting for Nim, with no dependencies beyond the standard library.

- Three-, four- and five-parameter logistic fits (3PL, 4PL / Hill, 5PL)
- Levenberg–Marquardt with multi-start, and a Nelder–Mead fallback
- Zero-dose controls, fixed parameters and weighted fits
- Robust fits (Huber, Tukey bisquare) that downweight or reject outliers
- Standard errors, t-based confidence intervals, R², AIC/BIC
- EC50 and any ECx (EC10, EC90, …) with confidence intervals
- Confidence and prediction bands for the fitted curve
- Results validated against R (`nls`, `robustbase::nlrob`, and the published `drc` estimates for its `ryegrass` dataset)

## Installation

Requires Nim ≥ 2.0. From a clone of this repository:

```sh
nimble install    # or: nimble develop
nimble test       # run the test suite
```

## The models

Responses are fitted against x = log10(dose) with the five-parameter logistic:

```
f(x) = Bottom + (Top − Bottom) / (1 + 10^((LogEC50 − x) · Hill))^Asymmetry
```

| Function | Model | Free parameters |
|---|---|---|
| `fitLL4` | 4PL (Hill equation) | Bottom, Top, LogEC50, Hill; Asymmetry is fixed at 1 |
| `fitLL3` | 3PL | as 4PL, with Bottom fixed (at 0 by default) |
| `fitLL5` | 5PL | all five, for curves that are not symmetric about their midpoint |

- **Hill > 0** gives a rising curve and **Hill < 0** a falling one (inhibition).
  Fits are reported with Bottom < Top.
- **A dose of 0** is treated as the low-dose plateau, so untreated controls can
  stay in the data.
- **`fit.params`** is indexed by `Bottom`, `Top`, `LogEC50`, `Hill` and `Asymmetry`.
- **In a 5PL, LogEC50 is a location parameter**, not the log of the half-way
  dose. Use `fit.ec50` and `fit.ecx(p)`, which give the actual doses for
  every model.

## Examples

### Fit a curve and print a summary

```nim
import drc

# Ryegrass root length vs ferulic acid concentration (mM), from the R drc package
let doses = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.94, 0.94, 0.94, 1.88, 1.88, 1.88,
             3.75, 3.75, 3.75, 7.5, 7.5, 7.5, 15.0, 15.0, 15.0, 30.0, 30.0, 30.0]
let rootLength = [7.58, 8.0, 8.3285714, 7.25, 7.375, 7.9625,
                  8.3555556, 6.9142857, 7.75, 6.8714286, 6.45, 5.9222222,
                  1.925, 2.8857143, 4.2333333, 1.1875, 0.8571429, 1.0571429,
                  0.6875, 0.525, 0.825, 0.25, 0.22, 0.44]

let fit = fitLL4(doses, rootLength)
echo fit.summary
```

```
4PL fit to 24 observations

              Estimate   Std.Error     Lower 95%     Upper 95%
Bottom        0.481410    0.209506     0.0443880      0.918432
Top            7.79296    0.189733       7.39719       8.18874
LogEC50       0.485431   0.0263916      0.430379      0.540483
Hill          -2.98223    0.458429      -3.93849      -2.02596
EC50           3.05796                   2.69389       3.47123

Residual standard error: 0.519626 on 20 degrees of freedom
R²: 0.978252   AIC: 42.3103   BIC: 48.2006
Stopped after 7 iterations: relative change in SSR below ftol
```

### Normalised data with fixed plateaus

For responses normalised to controls, fix the plateaus at 0 and 100:

```nim
import std/strformat
import drc

let doses = [1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4]   # molar
let viability = [98.2, 95.1, 71.4, 22.9, 6.3, 2.0] # % of untreated control

let fit = fitLL4(doses, viability, fixed = [(Bottom, 0.0), (Top, 100.0)])
let ci = fit.ec50CI
echo &"IC50 = {fit.ec50 * 1e9:.1f} nM (95% CI {ci.lower * 1e9:.1f}–{ci.upper * 1e9:.1f} nM)"
echo &"Hill = {fit.params[Hill]:.2f} ± {fit.se[Hill]:.2f}"
```

```
IC50 = 272.8 nM (95% CI 227.9–326.5 nM)
Hill = -0.89 ± 0.04
```

Fixed parameters have a standard error of 0 and do not count towards the
degrees of freedom. A 3PL with only Bottom fixed is `fitLL3(doses, viability)`.

### Asymmetric curves (5PL)

```nim
import std/[math, strformat]
import drc

var doses = @[0.0, 0.0]                 # two untreated controls
for k in 0 .. 12:                       # 1 nM to 1 mM in half-log steps, duplicates
  let dose = pow(10.0, -9.0 + 0.5 * float(k))
  doses.add [dose, dose]
let responses = [5.0, 5.6, 8.34, 7.1, 10.65, 9.57, 16.18, 18.74, 22.64, 22.38, 37.01, 36.75,
                 54.87, 52.8, 75.55, 77.0, 86.17, 87.94, 89.67, 90.89, 90.97, 94.18, 92.39,
                 95.46, 95.3, 94.61, 89.96, 93.92]

let four = fitLL4(doses, responses)
let five = fitLL5(doses, responses)
echo &"AIC: 4PL {four.aic:.1f}, 5PL {five.aic:.1f}"
echo five.summary
```

```
AIC: 4PL 127.5, 5PL 113.3
5PL fit to 28 observations

              Estimate   Std.Error     Lower 95%     Upper 95%
Bottom         4.88039     1.01725       2.77605       6.98472
Top            93.2261    0.580843       92.0246       94.4277
LogEC50       -6.08146   0.0909446      -6.26960      -5.89333
Hill           1.21811    0.144526      0.919131       1.51708
Asymmetry     0.392014   0.0764095      0.233949      0.550079
EC50       2.26380e-07               2.01742e-07   2.54026e-07

Residual standard error: 1.63042 on 23 degrees of freedom
R²: 0.998314   AIC: 113.327   BIC: 121.321
Stopped after 12 iterations: relative change in SSR below ftol
```

The lower AIC says the asymmetry is real. Note that EC50 (226 nM) is not
10^LogEC50 (829 nM) here. 5PL parameters are strongly correlated and need
data across both plateaus; prefer the 4PL unless the asymmetry is clear.

### ECx values, predictions and model comparison

```nim
import std/strformat
import drc

let doses = [0.0, 0.0, 0.94, 0.94, 1.88, 1.88, 3.75, 3.75, 7.5, 7.5, 15.0, 15.0, 30.0, 30.0]
let resp = [7.58, 8.0, 8.36, 6.91, 6.87, 6.45, 1.93, 2.89, 1.19, 0.86, 0.69, 0.53, 0.25, 0.22]
let fit = fitLL4(doses, resp)

for pct in [10.0, 50.0, 90.0]:
  let ci = fit.ecxCI(pct)
  echo &"EC{pct:g} = {fit.ecx(pct):.3f}  [{ci.lower:.3f}, {ci.upper:.3f}]"

echo &"predicted response at dose 5: {fit.predict(5.0):.3f}"

# Does a free Hill slope fit better than a standard slope of -1?
let unitHill = fitLL4(doses, resp, fixed = [(Hill, -1.0)])
echo &"AIC free Hill: {fit.aic:.2f}, Hill = -1: {unitHill.aic:.2f}"
```

```
EC10 = 1.618  [1.216, 2.154]
EC50 = 2.885  [2.511, 3.315]
EC90 = 5.144  [4.057, 6.521]
predicted response at dose 5: 1.335
AIC free Hill: 23.53, Hill = -1: 47.86
```

`ecx(fit, p)` is the dose giving p% of the change from the low-dose plateau to
the high-dose plateau. ECx intervals are computed on the log-dose scale and
back-transformed, so they are asymmetric in dose.

### Confidence and prediction bands

Continuing the example above:

```nim
let grid = [0.5, 1.0, 2.0, 4.0, 8.0, 16.0]
let conf = fit.confidenceBand(grid)
let pred = fit.predictionBand(grid)
echo "dose   fitted   95% confidence    95% prediction"
for i, dose in grid:
  echo &"{dose:>4g} {conf[i].fit:>8.3f}   [{conf[i].lower:.3f}, {conf[i].upper:.3f}]   [{pred[i].lower:6.3f}, {pred[i].upper:.3f}]"
```

```
dose   fitted   95% confidence    95% prediction
 0.5    7.774   [7.225, 8.324]   [ 6.603, 8.946]
   1    7.657   [7.190, 8.123]   [ 6.522, 8.791]
   2    6.342   [5.575, 7.109]   [ 5.054, 7.629]
   4    2.161   [1.472, 2.850]   [ 0.919, 3.404]
   8    0.684   [0.281, 1.088]   [-0.426, 1.794]
  16    0.548   [0.099, 0.996]   [-0.579, 1.675]
```

- **Confidence band:** where the true mean response lies at each dose. It uses
  the delta method with the parameter covariance.
- **Prediction band:** where a new observation is expected to fall. It adds
  the residual variance σ² (or σ²/`weight`, for a new observation with that
  relative weight).
- Both take a single dose or a sequence (as here, for plotting), and return
  `(fit, lower, upper)` tuples.

### Weights

Weights are relative, as in R's `nls`: observation i is taken to have error
variance σ²/wᵢ, and σ is estimated from the fit. A weight of 0 excludes an
observation.

```nim
let fit = fitLL4(doses, resp, weights = [1.0, 1.0, 2.0, 2.0, 1.0, 1.0, 1.0,
                                         1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0])
```

### Robust fitting (outliers)

Least squares lets a few bad wells drag the curve. A robust fit downweights
large residuals, either gently (`Huber`) or rejecting them outright beyond a
cut-off (`Tukey` bisquare):

```nim
import std/strformat
import drc

# % activity of an inhibitor, triplicates, with three wells gone wrong
let doses = [1e-9, 1e-9, 1e-9, 3.16e-9, 3.16e-9, 3.16e-9, 1e-8, 1e-8, 1e-8,
             3.16e-8, 3.16e-8, 3.16e-8, 1e-7, 1e-7, 1e-7, 3.16e-7, 3.16e-7, 3.16e-7,
             1e-6, 1e-6, 1e-6, 3.16e-6, 3.16e-6, 3.16e-6, 1e-5, 1e-5, 1e-5,
             3.16e-5, 3.16e-5, 3.16e-5, 1e-4, 1e-4, 1e-4]
let activity = [104.03, 90.24, 99.16, 95.92, 136.26, 96.97, 90.44, 95.81, 93.91,
                102.27, 92.98, 91.24, 77.88, 76.72, 75.56, 48.83, 16.45, 49.28,
                24.15, 20.67, 21.35, 12.34, 9.33, 6.18, 2.95, 35.12, 9.3,
                1.57, 1.65, 5.39, -0.56, 1.22, 4.74]

for loss in [Squared, Huber, Tukey]:
  let fit = fitLL4(doses, activity, opts = FitOptions(loss: loss))
  let ci = fit.ec50CI
  echo &"{$loss:<15} IC50 {fit.ec50 * 1e9:5.1f} nM  (95% CI {ci.lower * 1e9:5.1f}–{ci.upper * 1e9:5.1f})"

let tukey = fitLL4(doses, activity, opts = FitOptions(loss: Tukey))
for i, w in tukey.robustWeights:
  if w < 0.5: echo &"well {i}: dose {doses[i]:g}, activity {activity[i]}, weight {w:.2f}"
```

```
least squares   IC50 206.8 nM  (95% CI 137.9–310.2)
Huber           IC50 290.0 nM  (95% CI 250.1–336.2)
Tukey bisquare  IC50 310.4 nM  (95% CI 274.6–350.8)
well 0: dose 1e-09, activity 104.03, weight 0.32
well 1: dose 1e-09, activity 90.24, weight 0.35
well 4: dose 3.16e-09, activity 136.26, weight 0.00
well 9: dose 3.16e-08, activity 102.27, weight 0.00
well 16: dose 3.16e-07, activity 16.45, weight 0.00
well 25: dose 1e-05, activity 35.12, weight 0.00
```

The data were simulated with a true IC50 of 316 nM, and wells 4, 16 and 25 had
gross errors added. Both robust fits land much closer than least squares, with
narrower intervals. Tukey also rejected well 9, a genuine but unusually high
point. Rejection is automatic, so check `robustWeights` before trusting it.

How a robust fit works:
- **The method** is M-estimation by iteratively reweighted least squares, as
  R's `robustbase::nlrob`. Residuals are scaled by their robust spread
  (median |r| / 0.6745), and each round reweights and refits. Tukey starts
  from the Huber fit, because its loss has local minima.
- **`tuning`** in `FitOptions` sets the cut-off in units of that spread. The
  defaults (1.345 Huber, 4.685 Tukey) lose only 5% efficiency when there are
  no outliers.
- **Standard errors** use Huber's covariance correction (as `MASS::rlm`), and
  `fit.sigma` is the robust scale. `logLik`, `aic` and `bic` are NaN, since a
  robust fit is not a likelihood fit.
- **`fit.robustWeights`** holds each observation's final weight: 1 means not
  downweighted, 0 means rejected.

### Checking fit quality

Problems are reported in `fit.warnings` rather than raised:

```nim
let fit = fitLL4([0.01, 0.1, 1.0, 10.0], [99.0, 97.0, 90.0, 82.0], fixed = [(Bottom, 0.0)])
if not fit.converged or fit.warnings != {}:
  for w in fit.warnings: echo "warning: ", w
```

```
warning: EC50 lies outside the tested dose range
```

| Warning | Meaning |
|---|---|
| `fwNotConverged` | neither Levenberg–Marquardt nor the Nelder–Mead fallback converged |
| `fwEc50OutsideRange` | the EC50 is extrapolated beyond the tested doses |
| `fwSingular` | the parameters are not identifiable (e.g. a flat response); standard errors are NaN |
| `fwNoResidualDf` | as many free parameters as observations; standard errors are NaN |
| `fwRobustNotConverged` | the robust reweighting hit `irlsMaxIter` without settling |

Invalid input raises `ValueError`. This covers negative or non-finite doses,
non-finite responses, mismatched lengths, too few observations, a non-positive
Asymmetry, and confidence levels or percentages out of range.

### How fits are found

1. **Starting values** come from the data: the plateaus from the extreme mean
   responses, LogEC50 from where the means cross half-way, and Hill from a
   linearised regression.
2. **Multi-start:** Levenberg–Marquardt also runs from starts with LogEC50
   shifted by ±1 and Hill scaled by 0.4 and 2.5, and keeps the best fit. This
   catches local minima on sparse or noisy data. A 5PL is fitted symmetric
   first, then with Asymmetry released from 1, 0.5 and 2.
3. **Fallback:** if Levenberg–Marquardt still has not converged, Nelder–Mead
   (which needs no derivatives) minimises the SSR, and Levenberg–Marquardt
   polishes its result. `fit.usedFallback` records this.

Each 4PL fit takes about 0.2 ms with multi-start. For very large screens,
`FitOptions(multiStart: false)` is about 13× faster:

```nim
let fit = fitLL4(doses, viability, opts = FitOptions(multiStart: false))
```

`FitOptions` also holds the optimiser settings (`lm: LMOptions(...)` and
`nelderMead: NMOptions(...)`) and turns the fallback on or off (`fallback`).

### Fitting your own model

The Levenberg–Marquardt optimiser works with any model that returns its value
and writes its gradient with respect to the parameters:

```nim
import std/[math, strformat]
import drc

proc expDecay(x: float; p: array[2, float]; grad: var array[2, float]): float =
  let e = exp(-p[1] * x)
  grad = [e, -p[0] * x * e]    # ∂f/∂A, ∂f/∂k
  p[0] * e

let x = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
let y = [10.1, 6.0, 3.7, 2.2, 1.4, 0.8]
let w = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
let r = levenbergMarquardt(expDecay, x, y, w, p0 = [5.0, 0.1], free = [true, true])
echo &"A = {r.params[0]:.3f}, k = {r.params[1]:.3f}, converged: {r.status.converged}"
```

```
A = 10.070, k = 0.505, converged: true
```

For the parameter covariance, build JᵀWJ at the solution with
`normalEquations`, pass it to `unscaledCovariance`, and multiply by
σ² = SSR / (n − k).

## Layout

| File | Contents |
|---|---|
| `src/drc.nim` | public API: `fitLL3`/`fitLL4`/`fitLL5`, `DrcFit`, `FitOptions`, `summary`, `confint`, `ecx`, `ecxCI`, `predict`, `confidenceBand`, `predictionBand` |
| `src/drc/models.nim` | the 5PL model with its analytic gradient, and initial guesses |
| `src/drc/lm.nim` | Levenberg–Marquardt for weighted least squares |
| `src/drc/neldermead.nim` | Nelder–Mead simplex minimisation |
| `src/drc/robust.nim` | Huber and Tukey loss functions (ψ, weights, ψ′) and the MAD scale |
| `src/drc/stats.nim` | incomplete beta, Student's t CDF and quantile, covariance, log-likelihood |
| `src/drc/linalg.nim` | Cholesky solve and inverse for small fixed-size matrices |

## Validation

The tests compare against reference values from R:
- **Parameters, SEs and intervals:** `nls` fits of the ryegrass data (unweighted,
  weighted, and with Bottom fixed), compared on parameters, standard errors,
  σ, R², log-likelihood, AIC, BIC and confidence intervals; and a 5PL fit of
  synthetic asymmetric data, including ECx intervals.
- **Bands:** delta-method confidence and prediction bands, computed in R with
  finite-difference gradients, for 4PL and 5PL fits.
- **Robust fits:** `robustbase::nlrob` coefficients, scale and weights for
  Huber and Tukey fits of data with outliers. Standard errors are checked
  against `MASS::rlm`'s covariance formula, computed independently in R.
- **Distributions:** `qt`, `pt` and `pbeta`.
- **Cross-checks:** the parameter estimates also agree with SciPy's `curve_fit`
  and with the LL.4 and ED10/50/90 estimates published for `drc`.
- **Robustness:** a sparse dataset where a single start gets stuck in a local
  minimum (multi-start reaches the best SSR SciPy finds from 3000 random
  starts), and fits where Levenberg–Marquardt is starved of iterations.

Standard errors match `nls` and differ slightly from those `drc` reports
(Bottom 0.2095 vs 0.2122), likely because `drc` uses a numerical Hessian.
