package geometry

import org.scalatest.FunSuite

class FieldTest extends FunSuite {

  test("Scalar field construction") {
    def f(v: EuclideanVector[_1D]): Double = (v * 4).norm

    val dom = BoxDomain1D(EuclideanVector(0.0), EuclideanVector(5.0))
    val sf = ScalarField(dom, f)
    assert(sf.f(EuclideanVector(2.0)) == (EuclideanVector(2.0) * 4).norm)
  }

  test("Scalar field arithmetic") {
    def f(v: EuclideanVector[_2D]): Double = (v * 4).norm

    def g(v: EuclideanVector[_2D]): Double = (v * 2).norm

    val dom = BoxDomain2D(EuclideanVector(0.0, 0.0), EuclideanVector(5.0, 5.0))
    val sf = ScalarField(dom, f)
    val sg = ScalarField(dom, g)
    val vec = EuclideanVector(1.0, 1.0)
    assert((sf + sg).f(vec) == (sf.f(vec) + sg.f(vec)))
    assert((sf - sg).f(vec) == (sf.f(vec) - sg.f(vec)))
    assert((sf :* sg).f(vec) == (sf.f(vec) * sg.f(vec)))
  }

  test("Vector field construction") {
    def f(v: EuclideanVector[_1D]): EuclideanVector[_1D] = (v * 4) + EuclideanVector1D(3.0)

    val dom = BoxDomain1D(EuclideanVector(0.0), EuclideanVector(5.0))
    val vf = VectorField(dom, f)
    assert(vf.f(EuclideanVector(2.0)) == (EuclideanVector(2.0) * 4) + EuclideanVector1D(3.0))
  }

  test("Vector field arithmetic") {
    def f(v: EuclideanVector[_2D]): EuclideanVector[_2D] = v * 4

    def g(v: EuclideanVector[_2D]): EuclideanVector[_2D] = v * 3

    val dom = BoxDomain2D(EuclideanVector(0.0, 0.0), EuclideanVector(5.0, 5.0))
    val vf = VectorField(dom, f)
    val vg = VectorField(dom, g)
    val vec = EuclideanVector(1.0, 1.0)
    assert((vf + vg).f(vec) == (vf.f(vec) + vg.f(vec)))
    assert((vf - vg).f(vec) == (vf.f(vec) - vg.f(vec)))
    assert((vf * 4.0).f(vec) == (vf.f(vec) * 4.0))
  }
}
