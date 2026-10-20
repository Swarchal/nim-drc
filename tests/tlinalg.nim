import std/[math, unittest]
import drc/linalg

suite "linalg":
  test "solves an SPD system":
    let a: array[3, array[3, float]] = [[4.0, 12.0, -16.0], [12.0, 37.0, -43.0], [-16.0, -43.0, 98.0]]
    let b: array[3, float] = [1.0, 2.0, 3.0]
    var x: array[3, float]
    check solveSpd(a, b, x)
    for i in 0 ..< 3:
      var s = 0.0
      for j in 0 ..< 3: s += a[i][j] * x[j]
      check abs(s - b[i]) < 1e-10

  test "cholesky factor matches known result":
    let a: array[3, array[3, float]] = [[4.0, 12.0, -16.0], [12.0, 37.0, -43.0], [-16.0, -43.0, 98.0]]
    var l: array[3, array[3, float]]
    check cholesky(a, l)
    let expected: array[3, array[3, float]] = [[2.0, 0.0, 0.0], [6.0, 1.0, 0.0], [-8.0, 5.0, 3.0]]
    for i in 0 ..< 3:
      for j in 0 ..< 3:
        check abs(l[i][j] - expected[i][j]) < 1e-12

  test "rejects non positive definite matrices":
    var x: array[2, float]
    check not solveSpd([[1.0, 2.0], [2.0, 1.0]], [1.0, 1.0], x)
    check not solveSpd([[0.0, 0.0], [0.0, 1.0]], [1.0, 1.0], x)
    check not solveSpd([[NaN, 0.0], [0.0, 1.0]], [1.0, 1.0], x)

  test "SPD inverse":
    let a: array[3, array[3, float]] = [[4.0, 12.0, -16.0], [12.0, 37.0, -43.0], [-16.0, -43.0, 98.0]]
    var inv: array[3, array[3, float]]
    check invertSpd(a, inv)
    for i in 0 ..< 3:
      for j in 0 ..< 3:
        check inv[i][j] == inv[j][i]
        var s = 0.0
        for k in 0 ..< 3: s += a[i][k] * inv[k][j]
        check abs(s - (if i == j: 1.0 else: 0.0)) < 1e-10
    var bad: array[2, array[2, float]]
    check not invertSpd([[1.0, 2.0], [2.0, 1.0]], bad)
