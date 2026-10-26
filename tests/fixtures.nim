# Data sets shared by several tests.
import std/math
import ryegrass

const
  # Synthetic asymmetric data (truth: Bottom 5, Top 95, LogEC50 −6, Hill 1.3,
  # Asymmetry 0.35, noise sd 2), for 5PL fits.
  asymDose* = [
    0.0, 0.0, 1e-09, 1e-09, 3.1622776601683795e-09, 3.1622776601683795e-09, 1e-08,
    1e-08, 3.162277660168379e-08, 3.162277660168379e-08, 1e-07, 1e-07,
    3.162277660168379e-07, 3.162277660168379e-07, 1e-06, 1e-06, 3.162277660168379e-06,
    3.162277660168379e-06, 1e-05, 1e-05, 3.1622776601683795e-05, 3.1622776601683795e-05,
    0.0001, 0.0001, 0.00031622776601683794, 0.00031622776601683794, 0.001, 0.001,
  ]
  asymResponse* = [
    5.0, 5.6, 8.34, 7.1, 10.65, 9.57, 16.18, 18.74, 22.64, 22.38, 37.01, 36.75, 54.87,
    52.8, 75.55, 77.0, 86.17, 87.94, 89.67, 90.89, 90.97, 94.18, 92.39, 95.46, 95.3,
    94.61, 89.96, 93.92,
  ]
  # The README's simulated potency assay: a reference and two samples (4× and
  # 1/5× as potent), nine half-log doses in duplicate; see `potencyDesign`.
  potencyResponse* = [
    0.037, 0.07, 0.059, 0.157, 0.207, 0.161, 0.498, 0.508, 1.048, 1.013, 1.632, 1.587,
    1.832, 1.768, 2.002, 1.956, 1.987, 2.026, 0.138, 0.129, 0.269, 0.263, 0.574, 0.492,
    1.175, 1.148, 1.645, 1.593, 1.868, 1.904, 2.021, 1.963, 2.006, 1.988, 1.942, 1.981,
    0.036, 0.05, 0.101, 0.088, 0.069, 0.066, 0.167, 0.162, 0.306, 0.305, 0.799, 0.815,
    1.375, 1.415, 1.789, 1.749, 1.944, 1.885,
  ]

proc potencyDesign*(): tuple[doses: seq[float], curve: seq[int]] =
  ## The doses and curve ids for `potencyResponse`.
  for c in 0 .. 2:
    for k in 0 .. 8:
      let dose = pow(10.0, -2.0 + 0.5 * float(k))
      result.doses.add [dose, dose]
      result.curve.add [c, c]

proc ryegrassCurves*(
    scale: float
): tuple[doses, responses: seq[float], curve: seq[int]] =
  ## The ryegrass data twice: curve 0 as it is, and curve 1 with its doses
  ## multiplied by `scale` (so `scale` times less potent).
  for c in 0 .. 1:
    for i, d in ryegrassDose:
      result.doses.add d * (if c == 0: 1.0 else: scale)
      result.responses.add ryegrassRootLength[i]
      result.curve.add c
