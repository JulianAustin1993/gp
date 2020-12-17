package geometry

import org.scalatest.FunSuite

class DomainTest extends FunSuite {
  test("Reals Domain construction") {
    val domain = RealSpace[_1D]
    val vec = EuclideanVector(1.0)
    assert(domain.isDefinedAt(vec))
  }

  test("BoxDomain1D positive construction") {
    val domain = BoxDomain1D(EuclideanVector(0.0), EuclideanVector(5.0))
    assert(domain.isDefinedAt(EuclideanVector(2.0)))
    assert(!domain.isDefinedAt(EuclideanVector(-1.0)))
  }


  test("BoxDomain2D positive construction") {
    val domain = BoxDomain2D(EuclideanVector(0.0, 0.0), EuclideanVector(5.0, 4.0))
    assert(domain.isDefinedAt(EuclideanVector(2.0, 1.0)))
    assert(!domain.isDefinedAt(EuclideanVector(-1.0, 1.0)))
  }

  test("BoxDomain2D general construction") {
    val domain = BoxDomain(EuclideanVector(0.0, 0.0), EuclideanVector(5.0, 4.0))
    assert(domain.isDefinedAt(EuclideanVector(2.0, 1.0)))
    assert(!domain.isDefinedAt(EuclideanVector(-1.0, 1.0)))
  }

  test("BoxDomain2D negative construction") {
    assertThrows[IllegalArgumentException] {
      val domain = BoxDomain[_2D](EuclideanVector(0.0, 0.0), EuclideanVector(-5.0, -4.0))
    }
  }
}
