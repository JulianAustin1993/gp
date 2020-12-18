package statsmodels

import breeze.linalg.DenseVector
import breeze.numerics.abs
import breeze.stats.sampling.standardBasis
import geometry._
import kernel.Kernel.computeKernelMatrix
import kernel.StandardKernels.exponentiatedQuadraticKernel
import org.scalatest.FunSuite

class GaussianProcessTest extends FunSuite {

  test("GP zero mean construction and marginal") {
    val kernel = exponentiatedQuadraticKernel[_1D](1.0, math.sqrt(0.5))
    val gp = GaussianProcess(kernel)
    val xs = List(EuclideanVector(2.0), EuclideanVector(4.0), EuclideanVector(10.0))
    val marginalGp = gp.marginal(xs)
    val C = computeKernelMatrix(xs, kernel)
    assert(marginalGp == MultivariateNormalDistribution(DenseVector.zeros[Double](3), C))
  }

  test("GP non-zero mean construction and marginal") {
    val kernel = exponentiatedQuadraticKernel[_2D](5.0, 5.0)
    val meanField = ScalarField[_2D](BoxDomain2D(EuclideanVector(0.0, 0.0), EuclideanVector(5.0, 5.0)), (v: EuclideanVector[_2D]) => v(0) * 2 + v(1) * 2)
    val gp = GaussianProcess(meanField, kernel)
    val xs = List(EuclideanVector(2.0, 2.0), EuclideanVector(4.0, 3.0), EuclideanVector(4.9, 1.5))
    val marginalGp = gp.marginal(xs)
    val C = computeKernelMatrix(xs, kernel)
    assert(marginalGp == MultivariateNormalDistribution(DenseVector(xs.map(v => v(0) * 2 + v(1) * 2).toArray), C))
  }
  test("GP points outside domain") {
    val kernel = exponentiatedQuadraticKernel[_2D](5.0, 5.0)
    val meanField = ScalarField[_2D](BoxDomain2D(EuclideanVector(0.0, 0.0), EuclideanVector(5.0, 5.0)), (v: EuclideanVector[_2D]) => v(0) * 2 + v(1) * 2)
    val gp = GaussianProcess(meanField, kernel)
    val xs = List(EuclideanVector(2.0, 2.0), EuclideanVector(4.0, 3.0), EuclideanVector(5.01, 1.5))
    assertThrows[IllegalArgumentException] {
      gp.marginal(xs)
    }
  }

  test("test variance for GP") {
    val kernel = exponentiatedQuadraticKernel[_1D](1.0, 0.1)
    val gp = GaussianProcess(kernel)
    val numPoints = 3
    val domain = BoxDomain1D(EuclideanVector(0.0), EuclideanVector(5.0))
    val sampler = UniformBoxDomainSampler(domain, numPoints)
    val (pts, _) = sampler.sample().unzip
    val numSamples = 500

    def sampleValueForIthPoint(i: Int) = for (_ <- 0 until numSamples) yield {
      val ptData = gp.sampleAtPoints(pts)
      ptData(i)
    }

    def testAtIthPoint(i: Int) {
      val sampleValuesAtPt = sampleValueForIthPoint(i)
      val meanAtPt = sampleValuesAtPt.sum / numSamples
      val varAtPt = sampleValuesAtPt.foldLeft(0.0)((acc, e) => acc + (e - meanAtPt) * (e - meanAtPt)) / numSamples
      assert(abs(meanAtPt - 0.0) < 3e-1)
      assert(abs(varAtPt - 1.0) < 3e-1)
    }

    for (i <- 0 until numPoints) testAtIthPoint(i)
  }
}
