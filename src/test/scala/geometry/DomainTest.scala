package geometry

import geometry.Domain.{intersection, union}
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

  test("Domain union") {
    val d1 = BoxDomain1D(EuclideanVector(0.0), EuclideanVector(5.0))
    val d2 = BoxDomain1D(EuclideanVector(6.0), EuclideanVector(7.0))
    val d3 = union(d1, d2)
    assert(d3.isDefinedAt(EuclideanVector(4.5)) && d3.isDefinedAt(EuclideanVector(6.5)))
  }

  test("Domain intersection") {
    val d1 = BoxDomain1D(EuclideanVector(0.0), EuclideanVector(5.0))
    val d2 = BoxDomain1D(EuclideanVector(4.0), EuclideanVector(7.0))
    val d3 = intersection(d1, d2)
    assert(d3.isDefinedAt(EuclideanVector(4.5)))
    assert(!d3.isDefinedAt(EuclideanVector(6.0)))
  }
}
