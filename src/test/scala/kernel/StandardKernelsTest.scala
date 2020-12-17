package kernel

import breeze.linalg.{DenseMatrix, sum}
import breeze.numerics.abs
import geometry.{EuclideanVector, _1D}
import kernel.Kernel.computeKernelMatrix
import kernel.StandardKernels.{exponentiatedQuadraticKernel, gaussianKernel}
import org.scalatest.FunSuite

class StandardKernelsTest extends FunSuite {
  test("Validation of 1D gaussian kernel") {
    val kernel = gaussianKernel[_1D](1.0)
    val xs = List(EuclideanVector(1.0), EuclideanVector(2.0), EuclideanVector(5.0))
    val truth = DenseMatrix(
      (1.0, math.exp(-1.0), math.exp(-16.0)),
      (math.exp(-1.0), 1.0, math.exp(-9.0)),
      (math.exp(-16.0), math.exp(-9.0), 1.0))
    assert(sum(abs(computeKernelMatrix(xs, kernel) - truth)) < 1E-4)
  }

  test("Validation of 1D exponentiated quadratic kernel") {
    val kernel = exponentiatedQuadraticKernel[_1D](2.0, math.sqrt(0.5))
    val validKern = gaussianKernel[_1D](1.0)
    val xs = List(EuclideanVector(1.0), EuclideanVector(2.0), EuclideanVector(5.0))
    assert(sum(abs(computeKernelMatrix(xs, kernel) - 4.0 * computeKernelMatrix(xs, validKern))) < 1E-4)
  }

  test("Kernel arithmetic") {
    val kernel1 = exponentiatedQuadraticKernel[_1D](2.0, math.sqrt(0.5))
    val kernel2 = gaussianKernel[_1D](1.0)
    val xs = List(EuclideanVector(1.0), EuclideanVector(2.0), EuclideanVector(5.0))
    val K1 = computeKernelMatrix(xs, kernel1)
    val K2 = computeKernelMatrix(xs, kernel2) * 3.0
    val kern = kernel1 + (kernel2 * 3.0)
    assert(sum(abs(computeKernelMatrix(xs, kern) - (K1 + K2))) < 1E-4)
  }

  test("Kernel function arithmetic") {
    val kernel1 = exponentiatedQuadraticKernel[_1D](2.0, math.sqrt(0.5))

    def phi(x: EuclideanVector[_1D]): EuclideanVector[_1D] = x * 0.5

    val xs = List(EuclideanVector(2.0), EuclideanVector(4.0), EuclideanVector(10.0))
    val validK = 4.0 * computeKernelMatrix(xs.map(phi), gaussianKernel[_1D](1.0))
    val K = computeKernelMatrix(xs, kernel1.compose(phi))
    assert(sum(abs(K - validK)) < 1E-4)

  }

}
