package geometry

import geometry.Dim.{OneDSpace, TwoDSpace}

class DimTest extends org.scalatest.FunSuite {
  test("OneDSpace test") {
    val d = OneDSpace
    assert(d.dimensionality == 1)
  }

  test("TwoDSpace test") {
    val d = TwoDSpace
    assert(d.dimensionality == 2)
  }

  test("NDSpace test") {
    val d = NDSpace(OneDSpace)
    assert(d.dimensionality == 1)
  }

}
