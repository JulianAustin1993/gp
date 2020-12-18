package utils

import breeze.linalg.{DenseMatrix, DenseVector, isClose, lowerTriangular, sum, upperTriangular}
import breeze.numerics.abs
import org.scalatest.FunSuite


class utilsTest extends FunSuite {
  test("Cholesky Decomp") {
    val X = DenseMatrix.rand[Double](5, 5)
    val C = X.t * X
    val L = utils.regularisedCholeskyDecomposition(0.0)(C)
    assert(sum(abs((L * L.t) - C)) < 1E-4)
  }

  test("forward solve single rhs") {
    val A = lowerTriangular(DenseMatrix.rand[Double](5, 5))
    val x = DenseVector.rand[Double](5)
    val b = A * x
    assert(isClose(utils.forwardSolve(A, b), x))
  }

  test("forward solve multiple rhs") {
    val A = lowerTriangular(DenseMatrix.rand[Double](5, 5))
    val x = DenseMatrix.rand[Double](5, 3)
    val b = A * x
    assert(sum(abs(utils.forwardSolve(A, b) - x)) < 1E-4)
  }

  test("back solve single rhs") {
    val A = upperTriangular(DenseMatrix.rand[Double](5, 5))
    val x = DenseVector.rand[Double](5)
    val b = A * x
    assert(isClose(utils.backSolve(A, b), x))
  }

  test("back solve multiple rhs") {
    val A = upperTriangular(DenseMatrix.rand[Double](5, 5))
    val x = DenseMatrix.rand[Double](5, 3)
    val b = A * x
    assert(sum(abs(utils.backSolve(A, b) - x)) < 1E-4)
  }

}
