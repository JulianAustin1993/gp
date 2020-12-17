package kernel

import geometry.{Domain, EuclideanVector, RealSpace}

object StandardKernels {

  def gaussianKernel[D](sigma: Double) = new GaussianKernel[D](sigma = sigma)

  def exponentiatedQuadraticKernel[D](sigma: Double, lengthscale: Double) = new ExponentiatedQuadraticKernel[D](sigma = sigma, lengthscale = lengthscale)

}

class GaussianKernel[D](sigma: Double) extends PDKernel[D] {
  private val sigma2 = sigma * sigma

  override def domain: Domain[D] = RealSpace[D]

  override protected def k(x: EuclideanVector[D], y: EuclideanVector[D]): Double = {
    val r: EuclideanVector[D] = x - y
    scala.math.exp(-r.normSquared / sigma2)
  }
}

class ExponentiatedQuadraticKernel[D](sigma: Double, lengthscale: Double) extends PDKernel[D] {
  private val sigma2 = sigma * sigma
  private val l2 = 2 * lengthscale * lengthscale

  override def domain: Domain[D] = RealSpace[D]

  override protected def k(x: EuclideanVector[D], y: EuclideanVector[D]): Double = {
    val r: EuclideanVector[D] = x - y
    sigma2 * scala.math.exp(-r.normSquared / l2)
  }
}