# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

Dose-response curve fitting (3PL/4PL/5PL logistic) in Nim 2.x. **No dependencies beyond the Nim standard library** — this is a design goal; don't add nimble packages. C's `<math.h>` via `importc` is acceptable (see `log1p`/`expm1` in `models.nim`).

## Commands

```sh
nimble test                                        # build and run every tests/t*.nim
nim c -r --hints:off --path:src tests/tfit.nim     # one test file
nim c -r --hints:off --path:src tests/tfit.nim "fitLL4: ryegrass reference::zero doses dropped"   # one test ("suite::test", globs allowed)
nim c -r -d:release --path:src <file>.nim          # benchmarks / timing
nimble figures                                     # re-render docs/images (README figures)
nimble vignettes                                   # re-render docs/vignettes; `nim c -r vignettes/knit.nim intervals` for one
```

`figures` and `vignettes` need nim-plotter (`$PLOTTER_PATH`, default `../nim-plotter/src`) and `vl-convert` on PATH. The `dataframes` vignette also needs racoon (`$RACOON_PATH`, default `../Racoon/src`, else the installed package; the installed 0.1.0 prints tables without the dtype line, so render with the current checkout). plotter and racoon are not drc dependencies: only `examples/` and `vignettes/` import them. `vignettes/data/plate_screen.csv` is simulated data.

`nimble test` stops at the first failing file, so a failure can hide later files. Compiled test binaries land next to their sources in `tests/`.

## Architecture

The layers, from bottom to top:

- `src/drc/linalg.nim`: Cholesky solve and inverse, plus `pinFixed`, which decouples fixed parameters in a linear system. **Container-generic:** every routine takes fixed-size `array[N, …]` (single-curve fits) or runtime-sized `seq` (joint fits); `ensureLen`/`ensureShape`/`setZero` are no-ops for arrays, so the array code is what it was before.
- `src/drc/lm.nim` and `src/drc/neldermead.nim`: generic optimisers. `F` is any callable model `fn(x, p, var grad): float` (or `f(p): float` for Nelder–Mead), so closures work. LM's body is the private `lmCore`, with two public overloads: `[K: static int]` arrays → `LMResult[K]`, and `seq` → `LMDynResult`. `normalEquations` returns the true JᵀWJ. Pinning is applied only to the per-solve copy (in LM) or the per-covariance copy (`stats.unscaledCovariance`). Nelder–Mead is array-only.
- `src/drc/robust.nim`: Huber/Tukey ψ, weights ψ(u)/u, ψ′, and the MAD scale (median|r|/0.6745, as `nlrob`/`rlm`).
- `src/drc/stats.nim`: incomplete beta, Student t CDF/quantile (Newton inside a bisection bracket), F upper tail (`fTail`), covariance, Gaussian log-likelihood matching R's `logLik.nls`.
- `src/drc/models.nim`: **one model, `ll5`** (five-parameter logistic in x = log10 dose), plus `drcGuess` for starting values.
- `src/drc.nim`: public API. `fitCurve(model, …)` returns a `DrcFit`, and `fitLL3`/`fitLL4`/`fitLL5` are one-line wrappers over it with identical arguments. For LL3, `resolveFixed` adds Bottom = 0 unless `fixed` gives Bottom.
- `src/drc/multi.nim`: joint fits of several curves (`fitMulti(model, …)` and the `fitLL4Multi` etc. wrappers → `MultiFit`; `param`/`cov` accessors; `relativePotency`, `fTest`). It is `include`d at the end of `drc.nim` (not imported), so it can use that module's private helpers (`fitCurve`, `logEcxAndGradient`, `tQuantile`, `resolveFixed`, …).
- `src/drc/bootstrap.nim`: `bootstrap` for `DrcFit` and `MultiFit` → generic `Bootstrap[T]` (one resampling core, `resampleFits`), and `bootCI` (percentile, as `boot.ci(type = "perc")`). Refits start from the original estimates with one LM run (the private `fitCurveFrom`/`fitMultiFrom` with `warm`), falling back to a full fit if that does not converge. It is `include`d after `multi.nim`, for the same reason.

Design decisions that span files:

- **4PL is 5PL with `Asymmetry` fixed at 1, and 3PL is 4PL with `Bottom` fixed.** There is one parameter enum, `DrcParam`. Separate per-model enums can't work: Nim can't resolve duplicate enum member names inside `fixed = [(Bottom, 0.0)]`. For 3PL/4PL fits, `fit.fixed` therefore always contains `Asymmetry`.
- **Parameter vectors have two representations.** Optimisers and models use positional `array[drcParamCount, float]` indexed with `ord(Param)`. `DrcFit` exposes `array[DrcParam, float]`, and `toVec`/`toParams` convert at the boundary. Avoid generic type aliases with `static int` params (e.g. `Vec[N] = array[N, float]`): they broke type matching for `var` params.
- **Aliasing:** Nim passes arrays larger than 24 bytes by reference and builds a call's result in place, so in `r.fit = levenbergMarquardt(..., r.fit.params, ...)` the start and the result share memory. The array `levenbergMarquardt` is `{.noinit.}`, and `lmCore` resets its result only after copying `p0`, so this is safe (tested in `tests/tfit.nim`). A new proc that takes a large array and returns a plain object (no `seq`/`string` fields) needs the same treatment. Results with such fields go through a temporary, so `irls` and the `seq` LM overload are unaffected (tested in `tests/tmulti.nim`); never put `{.noinit.}` on those, since resetting an uninitialised seq frees garbage.
- **Zero doses** become x = -Inf. `ll5` returns the low-dose plateau limit there, so controls stay in the data. Hill = 0 is handled as a flat curve.
- **In 5PL, `LogEC50` is a location parameter, not log10(EC50).** `ec50`/`ecx` solve for the actual dose; `logEcxAndGradient` gives the delta-method gradient that includes Asymmetry.
- **Fit pipeline** (`fitCurve` → `searchLM`):
  1. Run LM multi-start: LogEC50 ±1 and Hill ×0.4/×2.5. Candidates are ranked by `better`: lower SSR first, convergence only breaks near-ties.
  2. For 5PL, fit symmetric first, then release Asymmetry.
  3. If nothing converged, run Nelder–Mead and polish its result with LM.
  4. For `FitOptions.loss` Huber/Tukey, run IRLS (`irls`, matching `robustbase::nlrob`): each round sets s = MAD scale of the √w-standardised residuals (floored at 1e-6 × the responses' weighted RMS spread), weights wᵢ·ψ(uᵢ)/uᵢ, and refits with LM from the current parameters. It converges only when a converged refit leaves the residuals unchanged, and stops early if a refit stalls. Tukey starts from a Huber IRLS fit run to nlrob's tolerance (1e-6).
  5. Compute covariance with `normalEquations` (user weights) at the final parameters, independent of which optimiser produced them. Least squares uses σ²(JᵀWJ)⁻¹. Robust fits use Huber's correction as in `MASS::rlm`, κ²·S/mean(ψ′)². For robust fits, `sigma` is the MAD scale and `logLik`/`aic`/`bic` are NaN.
- **Profile intervals** (`profileInterval`/`profileEnd`) take a τ(β) callback (`ProfileTau`). `pinnedTau` holds one parameter (or the reparametrised log10 ECx / log potency) at β. `constrainedTau` profiles any function g(p) by refitting subject to g(p) = β: an augmented Lagrangian, with the constraint as one extra observation (x = NaN) of weight μ, so the bounded LM handles bounds as usual. It is used where LogEC50 is fixed, bounded or shared, and returns τ = Inf where no fit within the bounds has g = β. If g depends on no free parameter, the interval is the point estimate (`constrainedInterval`). Profiles need least-squares fits; robust fits have no likelihood, so `logLik`/`aic`/`bic`, profiles, `fTest` and `lackOfFit` are unavailable for them by design.
- Confidence/prediction bands (`bandAt`) are gᵀ·cov·g from `ll5`'s gradient, plus σ²/weight for prediction.
- **Convergence:** `LMStatus.converged` includes `lmStationary` (no downhill step and gradient ≈ 0, or a near-perfect fit) but not `lmStalled`. After the fallback, `DrcFit.converged` can be true while `status` (the final LM polish) is not; `fallbackConverged` records the fallback's own outcome. `DrcFit.converged` is false exactly when `fwNotConverged` or `fwRobustNotConverged` is set.
- **Bounds** (`bounds = [(Hill, -10.0 .. -0.1)]`, a `Box` of ±Inf-padded `lower`/`upper` vectors internally): projected LM in `lmCore`. Each trial point is clamped into the box, and the active set (free parameters on a bound whose gradient g = JᵀWr points outwards) is pinned for that solve, as fixed ones are. `stationary` checks only the non-active parameters. Nelder–Mead clamps each vertex. Starts are clamped in `fitCurve`/`fitMulti` before `searchLM`, because the symmetric 5PL stage holds Asymmetry at its start value, which must lie within the box. A parameter can't be both fixed and bounded. Parameters that end exactly on a bound go in `atBound` (and raise `fwAtBound`), but they still count as free in df/AIC/covariance, matching `nls(algorithm = "port")` and SciPy.
- **Fixed parameters:** standard error 0, covariance rows 0, and they don't count toward `df`/AIC. Standard errors are NaN when unidentifiable or when df = 0, and problems are reported in `fit.warnings`, not raised. Invalid input raises `ValueError`.
- **Weights** are relative (variance σ²/wᵢ, as in R's `nls`). Zero-weight points are excluded from `n`, `df`, AIC and the starting guess.
- **Joint fits** (`fitMulti`):
  - The joint vector has one slot per shared parameter and one per curve for each unshared one, ordered by `DrcParam` then curve; the private `slots[c][param]` is the slot in the private `joint`/`jointFree`/`jointCov`. The public view is `fit.param(c, p)` and `fit.cov(a, pa, b, pb)`. Fixed parameters are global, count as shared, and are pinned with the usual machinery.
  - LM sees the observation index as x, and a closure maps the joint vector to that observation's curve's 5-vector, calls `ll5` and scatters the gradient.
  - Starts come from per-curve `fitCurve` fits (falling back to `drcGuess` for curves too small to fit alone), with shared parameters at the median over curves, plus a start from the pooled fit. The multi-start grid, the symmetric-then-asymmetric 5PL stage, the Nelder–Mead fallback (the `seq` overload of `nelderMead`) and the robust IRLS stage mirror `fitCurve`. The per-curve start fits are least squares even for a robust joint fit.
  - `fit.curves[c]` is a `DrcFit` view: curve c's 5×5 covariance block, the joint df and pooled σ (the robust scale for a robust fit), its own n/ssr/R²/dose range and robust weights. `logLik`/`aic`/`bic` are NaN on views, so single-curve inference works on them unchanged. A view has no link to its `MultiFit`, so profiles and the bootstrap of a curve go through the `MultiFit` (`profileCI(fit, curve, param)`, `bootCI` with `f.curves[c]`).

## Vignettes

- `vignettes/<name>.nim` are literate programs: column-0 `## ` lines are Markdown prose, everything else is code. `vignettes/knit.nim` runs each one and writes `docs/vignettes/<name>.md` with every chunk's code and real output, plus `docs/vignettes/README.md` (an index built from each title and first paragraph). Never edit the generated `.md` files. The order is the `vignettes` list in `knit.nim`.
- Code from a column-0 `#| hide` line to the next prose or `#| hide` line runs but isn't shown. Prose and `#| hide` lines must not split a statement.
- `vignettes/drcplot.nim` holds the plotting helpers, shared with `examples/figures.nim`. `figure(chart, name, alt)` writes `docs/vignettes/figures/<vignette>/<name>.png` and prints a marker that knit turns into an image link.
- Output must be deterministic (fixed seeds, no timings), so re-rendering only changes what actually changed.

## Testing conventions

- Numeric results are validated against **external references**, not against this code's own output. The reference values are hard-coded in the tests:
  - R `nls` (joint fits use indexed parameters, e.g. `l[c]`; bounded fits use `algorithm = "port"`), `anova`, `qt`, `pt`, `pbeta` and `pf`.
  - `robustbase::nlrob` for robust fits.
  - SciPy `curve_fit`.
  - `drc`'s published ryegrass estimates (the data are in `tests/ryegrass.nim`).
- R's `drc` package is not installed. Base R (`Rscript`, with `MASS` and `robustbase`) and Python with SciPy are.
- When adding a numeric feature, compute the reference independently: R/SciPy, or finite differences for gradients. `tests/tmodels.nim` checks every analytic gradient against central differences.
- Synthetic-data tests write out the curve formula on purpose, as an independent check on `ll5`.
- R's `nls` stops at its relative-offset tolerance (1e-5, and it fails to converge if tightened on these data), so it can be ~1e-7 short of the optimum. Compare estimates against it in units of their standard error (`near` in `tests/tmulti.nim`, 2e-5·SE), not relative to the value; a Bottom near 0 fails relative checks.
- **Keep README examples truthful.** Their output blocks are copied from real runs, so re-run the snippets after behaviour changes.


## Formatting conventions

- Code should be formatted with `nph` or if `nph` is not available, manually formatted to a similar style.

## Version control
- Do not make any commits.
