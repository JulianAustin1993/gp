package geometry

trait Domain[D] {
  def isDefinedAt(v: EuclideanVector[D]): Boolean
}

object Domain {
  def fromPredicate[D](chi: EuclideanVector[D] => Boolean): Domain[D] = (v: EuclideanVector[D]) => chi(v)

  def intersection[D](thisDomain: Domain[D], thatDomain: Domain[D]): Domain[D] = (v: EuclideanVector[D]) => thisDomain.isDefinedAt(v) && thatDomain.isDefinedAt(v)

  def union[D](thisDomain: Domain[D], thatDomain: Domain[D]): Domain[D] = (v: EuclideanVector[D]) => thisDomain.isDefinedAt(v) || thatDomain.isDefinedAt(v)
}

class RealSpace[D] extends Domain[D] {
  override def isDefinedAt(v: EuclideanVector[D]): Boolean = true
}

object RealSpace {
  def apply[D] = new RealSpace[D]
}

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

case class BoxDomain1D(origin: EuclideanVector1D, oppositeCorner: EuclideanVector1D) extends BoxDomain[_1D] {
  require(origin.x <= oppositeCorner.x, "Origin must be smaller than oppositeCorner")

  override def isDefinedAt(v: EuclideanVector[_1D]): Boolean = {
    val vec: EuclideanVector1D = v
    vec.x >= origin.x && vec.x <= oppositeCorner.x
  }
}

case class BoxDomain2D(origin: EuclideanVector2D, oppositeCorner: EuclideanVector2D) extends BoxDomain[_2D] {
  require(origin.x <= oppositeCorner.x && origin.y <= oppositeCorner.y, "Origin must be the lower left corner")

  override def isDefinedAt(v: EuclideanVector[_2D]): Boolean = {
    val vec: EuclideanVector2D = v
    vec.x >= origin.x && vec.x <= oppositeCorner.x && vec.y >= origin.y && vec.y <= oppositeCorner.y
  }
}

object BoxDomain {
  def apply(origin: EuclideanVector1D, oppositeCorner: EuclideanVector1D): BoxDomain1D = BoxDomain1D(origin, oppositeCorner)

  def apply(origin: EuclideanVector2D, oppositeCorner: EuclideanVector2D): BoxDomain2D = BoxDomain2D(origin, oppositeCorner)

  def apply[D: NDSpace](orig: EuclideanVector[D], oppCorner: EuclideanVector[D]): BoxDomain[D] = new BoxDomain[D] {
    require((0 until origin.dimensionality).forall(i => orig(i) <= oppCorner(i)), " Origin must be the lower left corner")
    override lazy val origin: EuclideanVector[D] = orig
    override lazy val oppositeCorner: EuclideanVector[D] = oppCorner
  }
}