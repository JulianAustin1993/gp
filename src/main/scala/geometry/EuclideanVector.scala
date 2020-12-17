package geometry

import breeze.linalg.DenseVector
import spire.algebra.{Field, VectorSpace}
import spire.std.any.DoubleAlgebra

import scala.language.implicitConversions


abstract class EuclideanVector[D: NDSpace] {
  def apply(i: Int): Double

  def dimensionality: Int = implicitly[NDSpace[D]].dimensionality

  /* Scalar arithmetic */
  def *(s: Double): EuclideanVector[D]

  def *:(s: Double): EuclideanVector[D] = this * s

  def unary_- : EuclideanVector[D] = this * (-1.0)

  /* Vector arithmetic */
  def +(that: EuclideanVector[D]): EuclideanVector[D]

  def -(that: EuclideanVector[D]): EuclideanVector[D]

  def :*(that: EuclideanVector[D]): EuclideanVector[D]

  def dot(that: EuclideanVector[D]): Double

  def normSquared: Double

  def normalize: EuclideanVector[D] = this / norm

  def /(s: Double): EuclideanVector[D] = this * (1.0 / s)

  def norm: Double = math.sqrt(normSquared)

  /* Conversions */
  def toArray: Array[Double]

  def toBreezeVector: DenseVector[Double] = DenseVector(toArray)

  /* functions */
  def mapWithIndex(f: (Double, Int) => Double): EuclideanVector[D]

  def map(f: Double => Double): EuclideanVector[D] = mapWithIndex((v, _) => f(v))

}

case class EuclideanVector1D(x: Double) extends EuclideanVector[_1D] {
  override def apply(i: Int): Double = i match {
    case 0 => x
    case _ => throw new IndexOutOfBoundsException("EuclideanVector1D has only 1 element.")
  }

  override def *(s: Double): EuclideanVector1D = EuclideanVector1D(x * s)

  override def +(that: EuclideanVector[_1D]): EuclideanVector1D = EuclideanVector1D(x + that.x)

  override def -(that: EuclideanVector[_1D]): EuclideanVector1D = EuclideanVector1D(x - that.x)

  override def :*(that: EuclideanVector[_1D]): EuclideanVector1D = EuclideanVector1D(x * that.x)

  override def dot(that: EuclideanVector[_1D]): Double = x * that.x

  override def normSquared: Double = x * x

  override def toArray: Array[Double] = Array(x)

  override def mapWithIndex(f: (Double, Int) => Double): EuclideanVector1D = EuclideanVector1D(f(x, 0))
}

case class EuclideanVector2D(x: Double, y: Double) extends EuclideanVector[_2D] {
  override def apply(i: Int): Double = i match {
    case 0 => x
    case 1 => y
    case _ => throw new IndexOutOfBoundsException("EuclideanVector2D has only 2 elements.")
  }

  override def *(s: Double): EuclideanVector2D = EuclideanVector2D(x * s, y * s)

  override def +(that: EuclideanVector[_2D]): EuclideanVector2D = EuclideanVector2D(x + that.x, y + that.y)

  override def -(that: EuclideanVector[_2D]): EuclideanVector2D = EuclideanVector2D(x - that.x, y - that.y)

  override def :*(that: EuclideanVector[_2D]): EuclideanVector2D = EuclideanVector2D(x * that.x, y * that.y)

  override def dot(that: EuclideanVector[_2D]): Double = x * that.x + y * that.y

  override def normSquared: Double = x * x + y * y

  override def toArray: Array[Double] = Array(x, y)

  override def mapWithIndex(f: (Double, Int) => Double): EuclideanVector2D = EuclideanVector2D(f(x, 0), f(y, 1))
}

object EuclideanVector1D {
  val zero: EuclideanVector1D = EuclideanVector1D(0.0)
  val ones: EuclideanVector1D = EuclideanVector1D(1.0)
  val unitX: EuclideanVector1D = EuclideanVector1D(1.0)
}

object EuclideanVector2D {
  val unitX: EuclideanVector2D = EuclideanVector2D(1.0, 0.0)
  val unitY: EuclideanVector2D = EuclideanVector2D(0.0, 1.0)

  val zero: EuclideanVector2D = EuclideanVector2D(0.0, 0.0)
  val ones: EuclideanVector2D = EuclideanVector2D(1.0, 1.0)
}

object EuclideanVector {

  def apply(x: Double): EuclideanVector[_1D] = EuclideanVector1D(x)

  def apply(x: Double, y: Double): EuclideanVector[_2D] = EuclideanVector2D(x, y)

  def fromBreezeVector[D: NDSpace](breeze: DenseVector[Double]): EuclideanVector[D] = {
    val dim = implicitly[NDSpace[D]].dimensionality
    require(breeze.size == dim, s"Invalid size of breeze vector (${breeze.size} != $dim)")
    EuclideanVector.apply[D](breeze.data)
  }

  def apply[D: NDSpace](d: Array[Double])(implicit builder: Create[D]): EuclideanVector[D] = builder.createVector(d)

  trait Create[D] {
    val zero: EuclideanVector[D]

    def createVector(data: Array[Double]): EuclideanVector[D]
  }

  trait Create1D extends Create[_1D] {
    override val zero: EuclideanVector[_1D] = EuclideanVector1D.zero

    override def createVector(data: Array[Double]): EuclideanVector[_1D] = {
      require(data.length == 1, "Creation of Vector failed: provided Array has invalid length")
      EuclideanVector1D(data(0))
    }
  }

  trait Create2D extends Create[_2D] {
    override val zero: EuclideanVector[_2D] = EuclideanVector2D.zero

    override def createVector(data: Array[Double]): EuclideanVector[_2D] = {
      require(data.length == 2, "Creation of Vector failed: provided Array has invalid length")
      EuclideanVector2D(data(0), data(1))
    }
  }

  def zeros[D: NDSpace
  ](
     implicit builder: Create[D]
   ):
  EuclideanVector[D] = builder.zero


  /** spire VectorSpace implementation for Vector */
  implicit def spireVectorSpace[D: NDSpace]: VectorSpace[EuclideanVector[D], Double] = new VectorSpace[EuclideanVector[D], Double] {
    implicit override def scalar: Field[Double] = Field[Double]

    override def timesl(r: Double, v: EuclideanVector[D]): EuclideanVector[D] = v.map(f => f * r)

    override def negate(x: EuclideanVector[D]): EuclideanVector[D] = x.map(f => -f)

    override def zero: EuclideanVector[D] = zeros[D]

    override def plus(x: EuclideanVector[D], y: EuclideanVector[D]): EuclideanVector[D] =
      x.mapWithIndex((f, i) => f + y(i))
  }

  implicit def parametricToConcrete1D(p: EuclideanVector[_1D]): EuclideanVector1D = p.asInstanceOf[EuclideanVector1D]

  implicit def parametricToConcrete2D(p: EuclideanVector[_2D]): EuclideanVector2D = p.asInstanceOf[EuclideanVector2D]

}

