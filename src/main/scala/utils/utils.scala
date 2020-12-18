package utils

import breeze.linalg.{DenseMatrix, DenseVector, NotConvergedException, lowerTriangular}
import com.github.fommil.netlib.BLAS.{getInstance => blas}
import com.github.fommil.netlib.LAPACK.{getInstance => lapack}
import org.netlib.util.intW
import spire.implicits.cfor

object utils {

  /**
   * Regularised Cholesky decomposition of symmetric matrix. Note we don't check for symmetry or non-empty.
   * Regularisation can be skipped by setting to 0.0
   * Modified from:
   * https://github.com/scalanlp/breeze/blob/master/math/src/main/scala/breeze/linalg/functions/cholesky.scala
   *
   * @param regWeight Regularisation weight to add to diagonal
   * @param X         Matrix to decompose using a regularised cholesky decomposition.
   * @return Lower triangular decomposition if try is successful.
   */
  def regularisedCholeskyDecomposition(regWeight: Double = 0.0)(X: DenseMatrix[Double]): DenseMatrix[Double] = {
    val A: DenseMatrix[Double] = lowerTriangular(X)
    if (regWeight != 0)
      cfor(0)(_ < A.cols, _ + 1) {
        i => A(i, i) += regWeight
      }
    val N = X.rows
    val info = new intW(0)
    lapack.dpotrf(
      "L" /* lower triangular */ ,
      N /* number of rows */ ,
      A.data,
      scala.math.max(1, N) /* LDA */ ,
      info
    )
    // A value of info.`val` < 0 would tell us that the i-th argument
    // of the call to dpotrf was erroneous (where i == |info.`val`|).
    assert(info.`val` >= 0)

    if (info.`val` > 0)
      throw new NotConvergedException(NotConvergedException.Iterations)

    A
  }

  /**
   * Forward solve a lower-triangular linear system
   * with a single RHS.
   * Credit:
   * https://github.com/darrenjw/scala-glm/blob/master/src/main/scala/scalaglm/Utils.scala
   *
   * @param A A lower-triangular matrix
   * @param y A single vector RHS
   * @return The solution, x, of the linear system A x = y
   */
  def forwardSolve(A: DenseMatrix[Double],
                   y: DenseVector[Double]): DenseVector[Double] = {
    val yc = y.copy
    blas.dtrsv("L", "N", "N", A.cols, A.toArray, A.rows, yc.data, 1)
    yc
  }

  /**
   * Forward solve an lower-triangular linear system
   * with a single RHS
   *
   * @param A An lower-triangular matrix
   * @param Y A matrix with columns corresponding to RHSs
   * @return Matrix of solutions, X, to the linear system A X = Y
   */
  def forwardSolve(A: DenseMatrix[Double], Y: DenseMatrix[Double]): DenseMatrix[Double] = {
    val yc = Y.copy
    blas.dtrsm("L", "L", "N", "N", yc.rows, yc.cols, 1.0, A.toArray, A.rows, yc.data, yc.rows)
    yc
  }

  /**
   * Back solve an upper-triangular linear system
   * with a single RHS
   * credit:
   * https://github.com/darrenjw/scala-glm/blob/master/src/main/scala/scalaglm/Utils.scala
   *
   * @param A An upper-triangular matrix
   * @param y A single vector RHS
   * @return The solution, x, of the linear system A x = y
   */
  def backSolve(A: DenseMatrix[Double],
                y: DenseVector[Double]): DenseVector[Double] = {
    val yc = y.copy
    blas.dtrsv("U", "N", "N", A.cols, A.toArray,
      A.rows, yc.data, 1)
    yc
  }

  /**
   * Back solve an upper-triangular linear system
   * with multiple RHSs
   * credit:
   * https://github.com/darrenjw/scala-glm/blob/master/src/main/scala/scalaglm/Utils.scala
   *
   * @param A An upper-triangular matrix
   * @param Y A matrix with columns corresponding to RHSs
   * @return Matrix of solutions, X, to the linear system A X = Y
   */
  def backSolve(A: DenseMatrix[Double], Y: DenseMatrix[Double]): DenseMatrix[Double] = {
    val yc = Y.copy
    blas.dtrsm("L", "U", "N", "N", yc.rows, yc.cols, 1.0, A.toArray, A.rows, yc.data, yc.rows)
    yc
  }
}
