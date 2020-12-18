/**
 * Implementation of Scalar Gaussian process
 * modified from:
 * https://github.com/unibas-gravis/scalismo/blob/master/src/main/scala/scalismo/statisticalmodel/GaussianProcess.scala
 */
package statsmodels

import breeze.linalg.{DenseMatrix, DenseVector}
import breeze.stats.distributions.RandBasis
import geometry._
import kernel.Kernel.{computeKernelMatrix, computeKernelVector}
import kernel.PDKernel
import spire.implicits.cforRange
import utils.utils.{backSolve, forwardSolve, regularisedCholeskyDecomposition}

import scala.util.Try

/**
 * Implementation of a scalar Gaussian process defined over N-dimensional space.
 *
 * @param mean        Mean field of the Gaussian process.
 * @param cov         Positive definite covariance kernel of the process.
 * @param NDSpace$D$0 Space of the process.
 * @tparam D Dimension of the process.
 */
class GaussianProcess[D: NDSpace](val mean: ScalarField[D], val cov: PDKernel[D]) {
  /** Sample the GP at point x. * */
  def sampleAtPoint(x: EuclideanVector[D])(implicit rand: RandBasis): DenseVector[Double] = {
    require(domain.isDefinedAt(x))
    val marginalDist = marginal(x)
    marginalDist.sample()(rand)
  }

  /** Calculate the marginal distribution at point x * */
  def marginal(x: EuclideanVector[D]): MultivariateNormalDistribution = {
    require(domain.isDefinedAt(x))
    val meanVec: DenseVector[Double] = DenseVector[Double](mean.f(x))
    val covMat: DenseMatrix[Double] = DenseMatrix(cov.apply(x, x))
    MultivariateNormalDistribution(meanVec, covMat)
  }

  /** Sample the GP at points xs. * */
  def sampleAtPoints(xs: Seq[EuclideanVector[D]])(implicit rand: RandBasis): DenseVector[Double] = {
    require(xs.forall(domain.isDefinedAt))
    val marginalDist = marginal(xs)
    marginalDist.sample()(rand)
  }

  /** Calculate the marginal distribution of points xs. * */
  def marginal(xs: Seq[EuclideanVector[D]]): MultivariateNormalDistribution = {
    require(xs.forall(domain.isDefinedAt))
    val meanVec: DenseVector[Double] = DenseVector[Double](xs.map(mean.f).toArray)
    val covMat: DenseMatrix[Double] = computeKernelMatrix(xs, cov)
    MultivariateNormalDistribution(meanVec, covMat)
  }

  def domain: Domain[D] = Domain.intersection(mean.domain, cov.domain)
}

object GaussianProcess {

  /** Gaussian process with mean and covariance */
  def apply[D: NDSpace](mean: ScalarField[D], cov: PDKernel[D]): GaussianProcess[D] = {
    new GaussianProcess[D](mean, cov)
  }

  /** Gaussian process with zero-mean and covariance */
  def apply[D: NDSpace](cov: PDKernel[D]): GaussianProcess[D] = {
    val zeroField = ScalarField[D](RealSpace[D], _ => 0.0)
    new GaussianProcess[D](zeroField, cov)
  }

  /**
   * Gaussian process regression with error .
   *
   * @param gp           Prior gaussian process with mean and covariance.
   * @param trainingData Sequence of triplet of point, value and error distribution.
   * @tparam D Dimension of space.
   * @return Tuple of posterior Gaussian process and log marginal likelihood of data.
   */
  def regression[D: NDSpace](gp: GaussianProcess[D], trainingData: IndexedSeq[(EuclideanVector[D], Double, MultivariateNormalDistribution)]): (GaussianProcess[D], Double) = {
    val (xs, ys, errorDists) = trainingData.unzip3

    val meanVec = DenseVector(xs.map(gp.mean).toArray)
    val yVec = DenseVector(ys.toArray)
    val fVec = yVec - meanVec

    val Ky = computeKernelMatrix(xs, gp.cov)
    cforRange(0 until Ky.rows) {
      i => Ky(i, i) += errorDists(i).cov(0, 0)
    }
    val L: DenseMatrix[Double] = {
      Iterator.iterate(1E-10)(w => w * 2)
        .map(w => Try {
          regularisedCholeskyDecomposition(w)(Ky)
        })
        .dropWhile(_.isFailure)
        .next()
        .get
    }

    val alpha = backSolve(L.t, forwardSolve(L, fVec))

    def Kxstar(x: EuclideanVector[D]): DenseMatrix[Double] = computeKernelVector(x, xs, gp.cov)

    def posteriorMean(x: EuclideanVector[D]): Double = {
      val m = Kxstar(x).t * alpha
      m(0)
    }

    val posteriorKernel = new PDKernel[D] {
      override def domain: Domain[D] = gp.domain

      override def k(x: EuclideanVector[D], y: EuclideanVector[D]): Double = {
        val v_x = forwardSolve(L, Kxstar(x))
        val v_y = forwardSolve(L, Kxstar(y))
        val K = gp.cov(x, y) - (v_x.t * v_y)
        K(0, 0)
      }
    }
    val logDetKy: Double = {
      var acc = 0.0
      cforRange(trainingData.indices) {
        i => acc += math.log(L(i, i))
      }
      2 * acc
    }
    val const = trainingData.length * math.log(2.0 * math.Pi)
    val posteriorGp = new GaussianProcess[D](ScalarField(gp.domain, posteriorMean), posteriorKernel)
    val logMarginalLikelihood = -0.5 * ((alpha.t * alpha) + logDetKy + const)
    (posteriorGp, logMarginalLikelihood)
  }

  /**
   *
   * Gaussian process liklihood of observed data.
   *
   * @param gp           Prior gaussian process with mean and covariance.
   * @param trainingData Sequence of triplet of point, value and error distribution.
   * @tparam D Dimension of space.
   * @return log marginal liklihood of data under gp.
   */
  def logMarginalLiklihood[D: NDSpace](gp: GaussianProcess[D], trainingData: IndexedSeq[(EuclideanVector[D], Double, MultivariateNormalDistribution)]): Double = {
    val (xs, ys, errorDists) = trainingData.unzip3

    val meanVec = DenseVector(xs.map(gp.mean).toArray)
    val yVec = DenseVector(ys.toArray)
    val fVec = yVec - meanVec

    val Ky = computeKernelMatrix(xs, gp.cov)
    cforRange(0 until Ky.rows) {
      i => Ky(i, i) += errorDists(i).cov(0, 0)
    }
    val L: DenseMatrix[Double] = {
      Iterator.iterate(1E-10)(w => w * 2)
        .map(w => Try {
          regularisedCholeskyDecomposition(w)(Ky)
        })
        .dropWhile(_.isFailure)
        .next()
        .get
    }

    val alpha = backSolve(L.t, forwardSolve(L, fVec))
    val logDetKy: Double = {
      var acc = 0.0
      cforRange(trainingData.indices) {
        i => acc += math.log(L(i, i))
      }
      2 * acc
    }
    val const = trainingData.length * math.log(2.0 * math.Pi)
    -0.5 * ((alpha.t * alpha) + logDetKy + const)
  }

}