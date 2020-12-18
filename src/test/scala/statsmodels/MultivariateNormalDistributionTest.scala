package statsmodels

import breeze.linalg.{DenseMatrix, DenseVector, sum}
import breeze.numerics.abs
import breeze.stats.sampling.standardBasis
import org.scalatest.FunSuite
import statsmodels.MultivariateNormalDistribution.estimateFromData
import utils.utils.regularisedCholeskyDecomposition

class MultivariateNormalDistributionTest extends FunSuite {
  test("Case class construction") {
    val X = DenseMatrix.rand[Double](5, 5)
    val C = X.t * X
    val mean = DenseVector.zeros[Double](5)
    val L = regularisedCholeskyDecomposition(0.0)(C)
    val mvn = MultivariateNormalDistribution(mean, C)
    assert(sum(abs(L - mvn.root)) < 1E-4)
  }
  test("Sample and estimate") {
    val X = DenseMatrix.rand[Double](5, 5)
    val C = X.t * X
    val mean = DenseVector.rand[Double](5)
    val L = regularisedCholeskyDecomposition(0.0)(C)
    val mvn = MultivariateNormalDistribution(mean, C)
    val nSamples = 100000
    val samples = for (i <- 1 to nSamples) yield mvn.sample()
    val estimated = estimateFromData(samples)
    assert(sum(abs(mean - estimated.mean)) < 1)
    assert(sum(abs(C - estimated.cov)) < 1)
  }

  test("pdf") {
    val C = DenseMatrix(1.0)
    val mean = DenseVector(0.0)
    val mvn = MultivariateNormalDistribution(mean, C)
    assert(math.abs(mvn.pdf(DenseVector(0.0)) - 0.3989422) < 1E-4)
    assert(math.abs(mvn.pdf(DenseVector(1.0)) - 0.24197) < 1E-4)
  }

  test("conditional") {
    val C = DenseMatrix((1.0, 0.5), (0.5, 1.0))
    val mean = DenseVector(1.0, 2.0)
    val mvn = MultivariateNormalDistribution(mean, C)
    val conditional = mvn.conditional(Vector((1, 2.0)))
    assert(math.abs(conditional.mean(0) - mean(0)) < 1E-4)
    assert(math.abs(conditional.cov(0, 0) - 0.75) < 1E-4)
  }

  test("marginal") {
    val C = DenseMatrix((1.0, 0.5), (0.5, 1.0))
    val mean = DenseVector(1.0, 2.0)
    val mvn = MultivariateNormalDistribution(mean, C)
    val marginal = mvn.marginal(Vector(1))
    assert(math.abs(marginal.mean(0) - mean(1)) < 1E-4)
    assert(math.abs(marginal.cov(0, 0) - C(1, 1)) < 1E-4)
  }
}
