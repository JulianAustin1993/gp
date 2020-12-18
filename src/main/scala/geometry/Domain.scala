/**
 * Domain classes and objects defining a domain implementation for models.
 * Modified from:
 * https://github.com/unibas-gravis/scalismo/blob/master/src/main/scala/scalismo/common/Domain.scala
 */
package geometry

/**
 * Domain trait, giving a domain with a check if a vector is in the domain.
 *
 * @tparam D Dimension of the space holding the domain.
 */
trait Domain[D] {
  def isDefinedAt(v: EuclideanVector[D]): Boolean
}

/**
 * Domain objects
 */
object Domain {
  /**
   * Create a domain from a predicate.
   *
   * @param chi predicate function defining which vectors are in the domain.
   * @tparam D Dimension of the domain
   * @return Domain where vectors are defined at chi(v)==True
   */
  def fromPredicate[D](chi: EuclideanVector[D] => Boolean): Domain[D] = (v: EuclideanVector[D]) => chi(v)

  /**
   * Intersection of domains
   *
   * @param thisDomain First domain
   * @param thatDomain Second domain
   * @tparam D dimension of the space
   * @return Domain which is the intersection between the two domains.
   */
  def intersection[D](thisDomain: Domain[D], thatDomain: Domain[D]): Domain[D] = (v: EuclideanVector[D]) => thisDomain.isDefinedAt(v) && thatDomain.isDefinedAt(v)

  /**
   * Union of domains.
   *
   * @param thisDomain First domain
   * @param thatDomain Second domain
   * @tparam D dimension of the space
   * @return Domain which is the union of the two domains.
   */
  def union[D](thisDomain: Domain[D], thatDomain: Domain[D]): Domain[D] = (v: EuclideanVector[D]) => thisDomain.isDefinedAt(v) || thatDomain.isDefinedAt(v)
}

/**
 * Real space domain implementation. Vectors are defined everywhere.
 *
 * @tparam D Dimension of the space holding the domain.
 */
class RealSpace[D] extends Domain[D] {
  override def isDefinedAt(v: EuclideanVector[D]): Boolean = true
}

/**
 * Companion object for RealSpace class.
 */
object RealSpace {
  def apply[D] = new RealSpace[D]
}

/**
 * Define a domain which is a box over the space.
 *
 * @tparam D Dimension of the space holding the domain.
 */
trait BoxDomain[D] extends Domain[D] {
  val origin: EuclideanVector[D]
  val oppositeCorner: EuclideanVector[D]
  val extent: EuclideanVector[D] = oppositeCorner - origin
  val volume: Double = (0 until origin.dimensionality).foldLeft(1.0)((prod, i) => prod * (oppositeCorner(i) - origin(i)))

  override def isDefinedAt(v: EuclideanVector[D]): Boolean = {
    def isInsideAxis(i: Int) = v(i) >= origin(i) && v(i) <= oppositeCorner(i)

    (0 until v.dimensionality).forall(i => isInsideAxis(i))
  }
}

/**
 * Case class of box domain in 1-dimension
 *
 * @param origin         lower left hand corner of domain.
 * @param oppositeCorner upper right hand corner of domain.
 */
case class BoxDomain1D(origin: EuclideanVector1D, oppositeCorner: EuclideanVector1D) extends BoxDomain[_1D] {
  require(origin.x <= oppositeCorner.x, "Origin must be smaller than oppositeCorner")

  override def isDefinedAt(v: EuclideanVector[_1D]): Boolean = {
    val vec: EuclideanVector1D = v
    vec.x >= origin.x && vec.x <= oppositeCorner.x
  }
}

/**
 * Case class of box domain in 2-dimensions.
 *
 * @param origin         lower left hand corner of domain.
 * @param oppositeCorner upper right hand corner of doomain.
 */
case class BoxDomain2D(origin: EuclideanVector2D, oppositeCorner: EuclideanVector2D) extends BoxDomain[_2D] {
  require(origin.x <= oppositeCorner.x && origin.y <= oppositeCorner.y, "Origin must be the lower left corner")

  override def isDefinedAt(v: EuclideanVector[_2D]): Boolean = {
    val vec: EuclideanVector2D = v
    vec.x >= origin.x && vec.x <= oppositeCorner.x && vec.y >= origin.y && vec.y <= oppositeCorner.y
  }
}

/**
 * Box domain apply methods.
 */
object BoxDomain {
  def apply(origin: EuclideanVector1D, oppositeCorner: EuclideanVector1D): BoxDomain1D = BoxDomain1D(origin, oppositeCorner)

  def apply(origin: EuclideanVector2D, oppositeCorner: EuclideanVector2D): BoxDomain2D = BoxDomain2D(origin, oppositeCorner)

  def apply[D: NDSpace](orig: EuclideanVector[D], oppCorner: EuclideanVector[D]): BoxDomain[D] = new BoxDomain[D] {
    require((0 until origin.dimensionality).forall(i => orig(i) <= oppCorner(i)), " Origin must be the lower left corner")
    override lazy val origin: EuclideanVector[D] = orig
    override lazy val oppositeCorner: EuclideanVector[D] = oppCorner
  }
}