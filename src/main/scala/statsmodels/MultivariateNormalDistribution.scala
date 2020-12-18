/**
 * Implementation of multivariate normal distribution.
 * Modified from:
 * https://github.com/unibas-gravis/scalismo/blob/master/src/main/scala/scalismo/statisticalmodel/MultivariateNormalDistribution.scala
 */
package statsmodels

import breeze.linalg.{DenseMatrix, DenseVector, diag}
import breeze.stats.distributions.{Gaussian, RandBasis}
import spire.implicits.cforRange
import utils.utils.{forwardSolve, regularisedCholeskyDecomposition}

import scala.util.Try

/**
 * Multivariate normal implementation
 *
 * @param mean mean vector
 * @param cov  covariance matrix
 */
case class MultivariateNormalDistribution(mean: DenseVector[Double], cov: DenseMatrix[Double]) {
  require(cov.rows == cov.cols)
  require(mean.size == cov.rows)

  /* Calculate cholesky of cov with regularisedCholeskyDecomposition */
  lazy val root: DenseMatrix[Double] = {
    Iterator.iterate(1E-10)(w => w * 2)
      .map(w => Try {
        regularisedCholeskyDecomposition(w)(cov)
      })
      .dropWhile(_.isFailure)
      .next()
      .get
  }
  /* Forward solve root for inverse calculation */
  lazy val rootInv: DenseMatrix[Double] = forwardSolve(root, DenseMatrix.eye[Double](dim))
  lazy val covInverse: DenseMatrix[Double] = {
    rootInv.t * rootInv
  }
  /* Log determinant calculated using cholesky decomposition */
  lazy val logCovDet: Double = {
    var acc = 0.0
    cforRange(0 until dim) {
      i => acc += math.log(root(i, i))
    }
    2 * acc
  }
  lazy val logNormFactor: Double = -0.5 * (logCovDet + dim * math.log(2.0 * math.Pi))
  val dim: Int = mean.size

  def pdf(x: DenseVector[Double]): Double = {
    math.exp(logpdf(x))
  }

  def logpdf(x: DenseVector[Double]): Double = {
    if (x.size != dim) throw new IllegalArgumentException(s"Invalid vector dimensionality (provided ${x.size} should be $dim)")
    val exponent = -0.5 * mahalanobisDistance2(x)
    logNormFactor + exponent
  }

  /* Squared mahalonobis distance */
  def mahalanobisDistance2(x: DenseVector[Double]): Double = {
    val alpha = rootInv * (x - mean)
    alpha dot alpha
  }

  /* sample form multivariate normal */
  def sample()(implicit rand: RandBasis): DenseVector[Double] = {
    val standardNormal = Gaussian(0, 1)(rand)
    val normalSamples = standardNormal.sample(dim)
    val u = DenseVector[Double](normalSamples.toArray)
    mean + (root * u)
  }

  /* Marginal distributino */
  def marginal(subspace: IndexedSeq[Int]): MultivariateNormalDistribution = {
    val subMean = mean(subspace).toDenseVector
    val subCov = cov(subspace, subspace).toDenseMatrix
    MultivariateNormalDistribution(subMean, subCov)
  }

  /* Conditional distribution */
  def conditional(observations: IndexedSeq[(Int, Double)]): MultivariateNormalDistribution = {
    val (obsIdx: IndexedSeq[Int], obsVals: IndexedSeq[Double]) = observations.unzip
    val unknownIdx = (0 until dim).filter(e => !obsIdx.contains(e))
    val meanUn = mean(unknownIdx).toDenseVector
    val meanObs = mean(obsIdx).toDenseVector

    val covUnUn = cov(unknownIdx, unknownIdx).toDenseMatrix
    val covUnObs = cov(unknownIdx, obsIdx).toDenseMatrix
    val covObsUn = covUnObs.t
    val covObsObs = cov(obsIdx, obsIdx).toDenseMatrix

    val diff = DenseVector(obsVals.toArray) - meanObs
    val LObsObs = {
      Iterator.iterate(1E-10)(w => w * 2)
        .map(w => Try {
          regularisedCholeskyDecomposition(w)(covObsObs)
        })
        .dropWhile(_.isFailure)
        .next()
        .get
    }
    val LObsObsInverse = forwardSolve(LObsObs, DenseMatrix.eye[Double](diff.size))
    val newMean = meanUn + covUnObs * (LObsObsInverse * diff)
    val newCov = covUnUn - covUnObs * LObsObsInverse * covObsUn
    MultivariateNormalDistribution(newMean, newCov)
  }
}

object MultivariateNormalDistribution {
  /**
   * Create multivariate normal from variance along specified vectors.
   *
   * @param mean     Mean vector
   * @param mainAxis Sequence of axis directions and variances.
   * @return Multivariate normal distribution.
   */
  def apply(mean: DenseVector[Double], mainAxis: Seq[(DenseVector[Double], Double)]): MultivariateNormalDistribution = {
    val dim = mean.size
    require(mainAxis.lengthCompare(dim) == 0)

    val cov: DenseMatrix[Double] = {
      val Phi = DenseMatrix.zeros[Double](dim, mainAxis.size)
      val sigma2 = DenseVector.zeros[Double](mainAxis.size)

      cforRange(mainAxis.indices) {
        i => {
          val (phi, sig2) = mainAxis(i)
          Phi(::, i) := phi
          sigma2(i) = sig2
        }
      }
      Phi * diag(sigma2) * Phi
    }
    MultivariateNormalDistribution(mean, cov)
  }

  /**
   * Estimate a multivariate normal distribution from observed data.
   *
   * @param samples observed samples of the data.
   * @return MultivariateNormalDistribution with mean and covariance estimated form the data.
   */
  def estimateFromData(samples: Seq[DenseVector[Double]]): MultivariateNormalDistribution = {
    val numSamples = samples.length
    require(numSamples > 0)
    val sampleDim = samples.head.length
    require(samples.forall(s => s.length == sampleDim))

    val zeroVec = DenseVector.zeros[Double](sampleDim)
    val sampleMean: DenseVector[Double] = {
      samples.foldLeft((zeroVec, 1))((acc, s) => (acc._1 + (s - acc._1) * (1.0 / acc._2), acc._2 + 1))._1
    }

    val zeroMatrix: DenseMatrix[Double] = DenseMatrix.zeros[Double](sampleDim, sampleDim)

    def outer(v1: DenseVector[Double], v2: DenseVector[Double]): DenseMatrix[Double] = v1.toDenseMatrix.t * v2.toDenseMatrix

    val ddof: Int = if (numSamples == 1) 1 else numSamples - 1
    val sampleCov = samples.foldLeft(zeroMatrix)((acc, s) => {
      val s0 = s - sampleMean
      acc + outer(s0, s0)
    }) * (1.0 / ddof)
    new MultivariateNormalDistribution(sampleMean, sampleCov)
  }
}
