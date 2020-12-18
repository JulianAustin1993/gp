/**
 * Dimension traits and objects to use as the basis for defining N-dimensional spaces.
 * Modified from:
 * https://github.com/unibas-gravis/scalismo/blob/master/src/main/scala/scalismo/geometry/Dim.scala
 */
package geometry

/**
 * Marker trait to distinguis dimensions.
 */
sealed trait Dim

/**
 * Marker for 1-dimensional space
 */
trait _1D extends Dim

/**
 * Marker for 2-dimensional space
 */
trait _2D extends Dim

/**
 * N-dimensional space of dimension D, with a vector creater
 *
 * @tparam D : Dimension of space.
 */
trait NDSpace[D] extends EuclideanVector.Create[D] {
  def dimensionality: Int
}

/**
 * Companion object for NDSpace trait.
 */
object NDSpace {
  def apply[D](implicit ndSpace: NDSpace[D]): NDSpace[D] = ndSpace
}

/**
 * Companion object for Dim trait
 */
object Dim {

  /**
   * Implicit 1-dimensional space.
   */
  implicit object OneDSpace extends NDSpace[_1D] with EuclideanVector.Create1D {
    override def dimensionality: Int = 1
  }

  /**
   * Implicit 2-dimensional space.
   */
  implicit object TwoDSpace extends NDSpace[_2D] with EuclideanVector.Create2D {
    override def dimensionality: Int = 2
  }

}