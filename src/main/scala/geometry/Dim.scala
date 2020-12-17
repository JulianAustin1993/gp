package geometry

sealed trait Dim

trait _1D extends Dim

trait _2D extends Dim


trait NDSpace[D] extends EuclideanVector.Create[D] {
  def dimensionality: Int
}

object NDSpace {
  def apply[D](implicit ndSpace: NDSpace[D]): NDSpace[D] = ndSpace
}

object Dim {

  implicit object OneDSpace extends NDSpace[_1D] with EuclideanVector.Create1D {
    override def dimensionality: Int = 1
  }

  implicit object TwoDSpace extends NDSpace[_2D] with EuclideanVector.Create2D {
    override def dimensionality: Int = 2
  }

}