# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

Dose-response curve fitting (3PL/4PL/5PL logistic) in Nim 2.x. **No dependencies beyond the Nim standard library** — this is a design goal; don't add nimble packages. C's `<math.h>` via `importc` is acceptable (see `log1p`/`expm1` in `models.nim`).

## Commands

```sh
nimble test                                        # build and run every tests/t*.nim
nim c -r --hints:off --path:src tests/tfit.nim     # one test file
nim c -r --hints:off --path:src tests/tfit.nim "fitLL4: ryegrass reference::zero doses dropped"   # one test ("suite::test", globs allowed)
nim c -r -d:release --path:src <file>.nim          # benchmarks / timing
```

`nimble test` stops at the first failing file, so a failure can hide later files. Compiled test binaries land next to their sources in `tests/`.

## Architecture

The layers, from bottom to top:

- `src/drc/linalg.nim`: Cholesky solve and inverse for fixed-size `array[K, array[K, float]]`, plus `pinFixed`, which decouples fixed parameters in a linear system.
- `src/drc/lm.nim` and `src/drc/neldermead.nim`: generic optimisers over `[K: static int; F]`. `F` is any callable model `fn(x, p, var grad): float` (or `f(p): float` for Nelder–Mead), so closures work. `normalEquations` returns the true JᵀWJ. Pinning is applied only to the per-solve copy (in LM) or the per-covariance copy (`stats.unscaledCovariance`).
- `src/drc/robust.nim`: Huber/Tukey ψ, weights ψ(u)/u, ψ′, and the MAD scale (median|r|/0.6745, as `nlrob`/`rlm`).
- `src/drc/stats.nim`: incomplete beta, Student t CDF/quantile (Newton inside a bisection bracket), covariance, Gaussian log-likelihood matching R's `logLik.nls`.
- `src/drc/models.nim`: **one model, `ll5`** (five-parameter logistic in x = log10 dose), plus `drcGuess` for starting values.
- `src/drc.nim`: public API. `fitLL3`/`fitLL4`/`fitLL5` are thin wrappers over a single private `fitCurve`. It returns a `DrcFit`.

Design decisions that span files:

- **4PL is 5PL with `Asymmetry` fixed at 1, and 3PL is 4PL with `Bottom` fixed.** There is one parameter enum, `DrcParam`. Separate per-model enums can't work: Nim can't resolve duplicate enum member names inside `fixed = [(Bottom, 0.0)]`. For 3PL/4PL fits, `fit.fixed` therefore always contains `Asymmetry`.
- **Parameter vectors have two representations.** Optimisers and models use positional `array[drcParamCount, float]` indexed with `ord(Param)`. `DrcFit` exposes `array[DrcParam, float]`, and `toVec`/`toParams` convert at the boundary. Avoid generic type aliases with `static int` params (e.g. `Vec[N] = array[N, float]`): they broke type matching for `var` params.
- **Aliasing:** Nim passes arrays larger than 24 bytes by reference and builds a call's result in place, so in `r.fit = levenbergMarquardt(..., r.fit.params, ...)` the start and the result share memory. `levenbergMarquardt` is `{.noinit.}` and resets `result` only after copying `p0`, so this is safe (tested in `tests/tfit.nim`). A new proc that takes a large array and returns a plain object (no `seq`/`string` fields) needs the same treatment; results with such fields go through a temporary, so `irls` is unaffected.
- **Zero doses** become x = -Inf. `ll5` returns the low-dose plateau limit there, so controls stay in the data. Hill = 0 is handled as a flat curve.
- **In 5PL, `LogEC50` is a location parameter, not log10(EC50).** `ec50`/`ecx` solve for the actual dose; `logEcxAndGradient` gives the delta-method gradient that includes Asymmetry.
- **Fit pipeline** (`fitCurve` → `searchLM`):
  1. Run LM multi-start: LogEC50 ±1 and Hill ×0.4/×2.5. Candidates are ranked by `better`: lower SSR first, convergence only breaks near-ties.
  2. For 5PL, fit symmetric first, then release Asymmetry.
  3. If nothing converged, run Nelder–Mead and polish its result with LM.
  4. For `FitOptions.loss` Huber/Tukey, run IRLS (`irls`, matching `robustbase::nlrob`): each round sets s = MAD scale of the √w-standardised residuals (floored at 1e-6 × the responses' weighted RMS spread), weights wᵢ·ψ(uᵢ)/uᵢ, and refits with LM from the current parameters. It converges only when a converged refit leaves the residuals unchanged, and stops early if a refit stalls. Tukey starts from a Huber IRLS fit run to nlrob's tolerance (1e-6).
  5. Compute covariance with `normalEquations` (user weights) at the final parameters, independent of which optimiser produced them. Least squares uses σ²(JᵀWJ)⁻¹. Robust fits use Huber's correction as in `MASS::rlm`, κ²·S/mean(ψ′)². For robust fits, `sigma` is the MAD scale and `logLik`/`aic`/`bic` are NaN.
- Confidence/prediction bands (`bandAt`) are gᵀ·cov·g from `ll5`'s gradient, plus σ²/weight for prediction.
- **Convergence:** `LMStatus.converged` includes `lmStationary` (no downhill step and gradient ≈ 0, or a near-perfect fit) but not `lmStalled`. After the fallback, `DrcFit.converged` can be true while `status` (the final LM polish) is not; `fallbackConverged` records the fallback's own outcome. `DrcFit.converged` is false exactly when `fwNotConverged` or `fwRobustNotConverged` is set.
- **Fixed parameters:** standard error 0, covariance rows 0, and they don't count toward `df`/AIC. Standard errors are NaN when unidentifiable or when df = 0, and problems are reported in `fit.warnings`, not raised. Invalid input raises `ValueError`.
- **Weights** are relative (variance σ²/wᵢ, as in R's `nls`). Zero-weight points are excluded from `n`, `df`, AIC and the starting guess.

## Testing conventions

- Numeric results are validated against **external references**, not against this code's own output. The reference values are hard-coded in the tests:
  - R `nls`, `qt`, `pt` and `pbeta`.
  - `robustbase::nlrob` for robust fits.
  - SciPy `curve_fit`.
  - `drc`'s published ryegrass estimates (the data are in `tests/ryegrass.nim`).
- R's `drc` package is not installed. Base R (`Rscript`, with `MASS` and `robustbase`) and Python with SciPy are.
- When adding a numeric feature, compute the reference independently: R/SciPy, or finite differences for gradients. `tests/tmodels.nim` checks every analytic gradient against central differences.
- Synthetic-data tests write out the curve formula on purpose, as an independent check on `ll5`.
- **Keep README examples truthful.** Their output blocks are copied from real runs, so re-run the snippets after behaviour changes.
