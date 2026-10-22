## Loss functions for robust M-estimation, fitted by iteratively reweighted
## least squares (IRLS). Residuals are measured in units of a robust scale,
## u = r / s, and a tuning constant c sets where the loss stops being quadratic.

import std/algorithm

type Loss* = enum
  Squared = "least squares"
  Huber = "Huber" ## quadratic within ±c, linear beyond: outliers are downweighted
  Tukey = "Tukey bisquare" ## redescending: residuals beyond ±c get weight 0

func defaultTuning*(loss: Loss): float =
  ## Tuning constants giving 95% efficiency when the errors are normal.
  case loss
  of Squared: Inf
  of Huber: 1.345
  of Tukey: 4.685

func psiWeight*(loss: Loss, u, c: float): float =
  ## The IRLS weight ψ(u)/u (1 at u = 0).
  case loss
  of Squared:
    1.0
  of Huber:
    (if abs(u) <= c: 1.0 else: c / abs(u))
  of Tukey:
    if abs(u) < c:
      (1.0 - (u / c) * (u / c)) * (1.0 - (u / c) * (u / c))
    else:
      0.0

func psi*(loss: Loss, u, c: float): float =
  ## ψ(u), the derivative of the loss. Written out rather than as u·ψ(u)/u,
  ## which would give NaN at u = ±Inf.
  case loss
  of Squared:
    u
  of Huber:
    clamp(u, -c, c)
  of Tukey:
    (if abs(u) < c: u * psiWeight(Tukey, u, c) else: 0.0)

func psiDerivative*(loss: Loss, u, c: float): float =
  ## ψ'(u).
  case loss
  of Squared:
    1.0
  of Huber:
    (if abs(u) <= c: 1.0 else: 0.0)
  of Tukey:
    let v = (u / c) * (u / c)
    if abs(u) < c:
      (1.0 - v) * (1.0 - 5.0 * v)
    else:
      0.0

func madScale*(residuals: openArray[float]): float =
  ## Robust scale median(|r|)/0.6745, as R's nlrob and MASS::rlm (residuals
  ## are centred on zero by the fit, so no median is subtracted).
  var a = newSeq[float](residuals.len)
  for i, r in residuals:
    a[i] = abs(r)
  a.sort()
  let m = a.len div 2
  let median =
    if a.len mod 2 == 1:
      a[m]
    else:
      0.5 * (a[m - 1] + a[m])
  median / 0.6745
