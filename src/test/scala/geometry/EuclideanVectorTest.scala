package geometry

import breeze.linalg.{DenseVector, isClose}
import org.scalatest.FunSuite

class EuclideanVectorTest extends FunSuite {
  test("EuclideanVector1D indexing") {
    val vec = EuclideanVector1D(1.4)
    assert(vec(0) == 1.4)
    assertThrows[IndexOutOfBoundsException] {
      vec(1)
    }
  }
  test("EuclideanVector1D unit, zeros and ones") {
    assert(EuclideanVector1D.ones(0) == 1)
    assert(EuclideanVector1D.zero(0) == 0)
    assert(EuclideanVector1D.unitX(0) == 1)
  }

  test("EuclideanVector2D indexing") {
    val vec = EuclideanVector2D(1.4, -1.0)
    assert(vec(0) == 1.4)
    assert(vec(1) == -1.0)
    assertThrows[IndexOutOfBoundsException] {
      vec(2)
    }
  }
  test("EuclideanVector2D unit, zeros and ones") {
    assert(EuclideanVector2D.ones(0) == 1)
    assert(EuclideanVector2D.ones(1) == 1)
    assert(EuclideanVector2D.zero(0) == 0)
    assert(EuclideanVector2D.zero(1) == 0)
    assert(EuclideanVector2D.unitX(0) == 1)
    assert(EuclideanVector2D.unitX(1) == 0)
    assert(EuclideanVector2D.unitY(0) == 0)
    assert(EuclideanVector2D.unitY(1) == 1)
  }

  test("EuclideanVector") {
    assert(EuclideanVector(0.0)(0) == 0)
    assert(EuclideanVector(1.0, 2.0)(1) == 2)
    val vec = EuclideanVector[_2D](Array(2.0, 4.0))
    assert(vec(1) == 4.0)
    assert(vec(0) == 2.0)
    val breeze = DenseVector[Double](2.0, 3.0)
    val vec2 = EuclideanVector.fromBreezeVector[_2D](breeze)
    assert(vec2(0) == breeze(0) && vec2(1) == breeze(1))
  }

  test("EuclideanVector[_1D] scala arithmetic") {
    val x = EuclideanVector[_1D](Array(2.0))
    assert((x * 4.0) (0) == 8.0)
    assert((x / 4.0) (0) == 0.5)
    assert((-x) (0) == -2.0)
  }
  test("EuclideanVector[_2D] scalar arithmetic") {
    val y = EuclideanVector[_2D](Array(2.0, 4.0))
    assert((y * 4.0) (0) == 8.0 && (y * 4.0) (1) == 16.0)
    assert((y / 4.0) (0) == 0.5 && (y / 4.0) (1) == 1.0)
    assert((-y) (0) == -2.0 && (-y) (1) == -4.0)
  }

  test("EuclideanVector[_1D] vector arithmetic") {
    val x = EuclideanVector[_1D](Array(2.0))
    val y = EuclideanVector[_1D](Array(3.0))
    assert((x + y) (0) == 5.0)
    assert((x - y) (0) == -1.0)
    assert((x :* y) (0) == 6.0)
    assert(x.normSquared == 4.0 && x.norm == 2.0)
    assert((x dot y) == 6.0)
    assert(x.normalize(0) == 1.0)
  }

  test("EuclideanVector[_2D] vector arithmetic") {
    val x = EuclideanVector[_2D](Array(2.0, 4.0))
    val y = EuclideanVector[_2D](Array(2.0, 2.0))
    assert((x + y) (0) == 4.0 && (x + y) (1) == 6.0)
    assert((x - y) (0) == 0.0 && (x - y) (1) == 2.0)
    assert((x :* y) (0) == 4.0 && (x :* y) (1) == 8.0)
    assert(x.normSquared == 20.0 && x.norm == math.sqrt(20))
    assert((x dot y) == 12.0)
    assert(y.normalize(0) == 2.0 / math.sqrt(8.0))
  }

  test("EuclideanVector conversions") {
    val x = EuclideanVector2D(2.0, 4.0)
    assert(x.toArray(0) == 2.0, x.toArray(1) == 4.0)
    assert(isClose(x.toBreezeVector, DenseVector[Double](2.0, 4.0)))
  }

  test("EuclideanVector mapping") {
    val x = EuclideanVector2D(1.0, 5.0)
    val y = EuclideanVector1D(6.0)
    val fx = x.map(v => math.pow(v, 3.0))
    val gy = y.map(v => math.pow(v, 2.0) + v - 3.0)
    assert(fx(0) == 1.0 && fx(1) == 125.0)
    assert(gy(0) == 39.0)
  }

  test("spireVectorSpace") {
    val x = EuclideanVector2D(1.0, 5.0)
    val y = EuclideanVector2D(6.0, 0.0)
    val h = 4.0 *: x + y
    assert(h(0) == 10.0 && h(1) == 20.0)
    assert(h == EuclideanVector2D(10.0, 20.0))
  }
}
