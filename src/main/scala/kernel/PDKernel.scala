/**
 * Positive Definite kernel implementation.
 * modified from:
 * https://github.com/unibas-gravis/scalismo/blob/master/src/main/scala/scalismo/kernels/Kernel.scala
 */
package kernel

import breeze.linalg.{DenseMatrix, diag}
import geometry.{Domain, EuclideanVector, RealSpace}
import spire.syntax.cfor._

/**
 * Abstract implementation of a positive definite kernel.
 *
 * @tparam D : Dimensionality of the space the kernel is defined over.
 */
abstract class PDKernel[D] {
  self =>

  def domain: Domain[D]

  def apply(x: EuclideanVector[D], y: EuclideanVector[D]): Double = if (this.domain.isDefinedAt(x) && this.domain.isDefinedAt(y))
    k(x, y)
  else if (!this.domain.isDefinedAt(x)) throw new IllegalArgumentException(s"$x is outside of the domain") else {
    throw new IllegalArgumentException(s"$y is outside of the domain")
  }

  def +(that: PDKernel[D]): PDKernel[D] = new PDKernel[D] {
    override protected def k(x: EuclideanVector[D], y: EuclideanVector[D]): Double = self.k(x, y) + that.k(x, y)

    override def domain: Domain[D] = Domain.intersection(self.domain, that.domain)
  }

  def *(that: PDKernel[D]): PDKernel[D] = new PDKernel[D] {
    override protected def k(x: EuclideanVector[D], y: EuclideanVector[D]): Double = self.k(x, y) * that.k(x, y)

    override def domain: Domain[D] = Domain.intersection(self.domain, that.domain)
  }

  def *(s: Double): PDKernel[D] = new PDKernel[D] {
    override protected def k(x: EuclideanVector[D], y: EuclideanVector[D]): Double = self.k(x, y) * s

    override def domain: Domain[D] = self.domain
  }

  def compose(phi: EuclideanVector[D] => EuclideanVector[D]): PDKernel[D] = new PDKernel[D] {
    override protected def k(x: EuclideanVector[D], y: EuclideanVector[D]): Double = self.k(phi(x), phi(y))

    override def domain: Domain[D] = self.domain
  }

  protected def k(x: EuclideanVector[D], y: EuclideanVector[D]): Double
}

/**
 * Zero kernel for utility.
 *
 * @tparam D : Dimensionality of the space the kernel is defined over.
 */
class ZeroKernel[D]() extends PDKernel[D] {
  override def domain: Domain[D] = RealSpace[D]

  override protected def k(x: EuclideanVector[D], y: EuclideanVector[D]): Double = 0.0
}

/**
 * Companion objects for the Kernel.
 */
object Kernel {
  def zeroKernel[D] = new ZeroKernel[D]()

  /**
   * Compute a kernel matrix from a sequence of points.
   *
   * @param xs     Sequence of euclidean vectors to calculate kernel for.
   * @param kernel The positive definite kernel to calculate.
   * @tparam D space of the kernel and vectors.
   * @return DenseMatrix of the kernel evalauted as all points xs(i), xs(j).
   */
  def computeKernelMatrix[D](xs: Seq[EuclideanVector[D]], kernel: PDKernel[D]): DenseMatrix[Double] = {
    val K = DenseMatrix.zeros[Double](xs.size, xs.size)
    cfor(0)(_ < xs.size, _ + 1) {
      i => {
        cfor(0)(_ <= i, _ + 1) {
          j => K(i, j) = kernel(xs(i), xs(j))
        }
      }
    }
    K + K.t - diag(diag(K))
  }

  /**
   *
   * Compute a kernel matrix from a sequence of points.
   *
   * @param xs     First sequence of euclidean vectors to calculate kernel for.
   * @param ys     Second sequence of euclidean vectors to calculate kernel for.
   * @param kernel The positive definite kernel to calculate.
   * @tparam D space of the kernel and vectors.
   * @return DenseMatrix of the kernel evalauted as all points xs(i), ys(j).
   */
  def computeKernelMatrix[D](xs: Seq[EuclideanVector[D]], ys: Seq[EuclideanVector[D]], kernel: PDKernel[D]): DenseMatrix[Double] = {
    val K = DenseMatrix.zeros[Double](xs.size, ys.size)
    cfor(0)(_ < xs.size, _ + 1) {
      i => {
        cfor(0)(_ < ys.size, _ + 1) {
          j => K(i, j) = kernel(xs(i), ys(j))
        }
      }
    }
    K
  }

  /**
   * Calculate a kernel matrix for single point x agains a vector of points.
   *
   * @param x      Euclidean vectors to calculate kernel for.
   * @param xs     Second sequence of euclidean vectors to calculate kernel for.
   * @param kernel The positive definite kernel to calculate.
   * @tparam D space of the kernel and vectors.
   * @return DenseMatrix of the kernel evalauted as all points x, xs(i).
   */
  def computeKernelVector[D](x: EuclideanVector[D], xs: Seq[EuclideanVector[D]], kernel: PDKernel[D]): DenseMatrix[Double] = {
    val K = DenseMatrix.zeros[Double](1, xs.size)
    cfor(0)(_ < xs.size, _ + 1) {
      i => K(0, i) = kernel(x, xs(i))
    }
    K
  }
}
