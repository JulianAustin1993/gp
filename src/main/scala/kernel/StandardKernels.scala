/**
 * Implementation of some standard kernels for ease of use.
 * Modified from:
 * https://github.com/unibas-gravis/scalismo/blob/master/src/main/scala/scalismo/kernels/StandardKernels.scala
 */
package kernel

import geometry.{Domain, EuclideanVector, RealSpace}

object StandardKernels {

  def gaussianKernel[D](sigma: Double) = new GaussianKernel[D](sigma = sigma)

  def exponentiatedQuadraticKernel[D](sigma: Double, lengthscale: Double) = new ExponentiatedQuadraticKernel[D](sigma = sigma, lengthscale = lengthscale)

}

/**
 * Simple Gaussian kernel over N-dimensional space.
 *
 * @param sigma lengthscale of kernel.
 * @tparam D : Dimensionality of the space the kernel is defined over.
 */
class GaussianKernel[D](sigma: Double) extends PDKernel[D] {
  private val sigma2 = sigma * sigma

  override def domain: Domain[D] = RealSpace[D]

  override protected def k(x: EuclideanVector[D], y: EuclideanVector[D]): Double = {
    val r: EuclideanVector[D] = x - y
    scala.math.exp(-r.normSquared / sigma2)
  }
}

/**
 * Isotropic Exponentiated Quadratic (Squared exponential) kernel over N-dimensional space.
 *
 * @param sigma       variance parameter.
 * @param lengthscale lengthscale parameter
 * @tparam D : Dimensionality of the space the kernel is defined over.
 */
class ExponentiatedQuadraticKernel[D](sigma: Double, lengthscale: Double) extends PDKernel[D] {
  private val sigma2 = sigma * sigma
  private val l2 = 2 * lengthscale * lengthscale

  override def domain: Domain[D] = RealSpace[D]

  override protected def k(x: EuclideanVector[D], y: EuclideanVector[D]): Double = {
    val r: EuclideanVector[D] = x - y
    sigma2 * scala.math.exp(-r.normSquared / l2)
  }
}