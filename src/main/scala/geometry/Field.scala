/**
 * Implementation of a field defined over a domain with function.
 * Modified from:
 * https://github.com/unibas-gravis/scalismo/blob/master/src/main/scala/scalismo/common/Field.scala
 */
package geometry

/**
 * Field trait giving an method to apply a function over a vector in the domain.
 *
 * @tparam D Dimension of the space of the domain
 * @tparam A Type giving the image that function over the domain defines.
 */
trait Field[D, A] extends (EuclideanVector[D] => A) {
  self =>
  val f: EuclideanVector[D] => A

  def domain: Domain[D]

  override def apply(v: EuclideanVector[D]): A = {
    if (!isDefinedAt(v)) throw new IllegalArgumentException(s"Vector $v is outside domain")
    f(v)
  }

  def isDefinedAt(v: EuclideanVector[D]): Boolean = domain.isDefinedAt(v)

  def liftValues: (EuclideanVector[D] => Option[A]) = new Field[D, Option[A]] {
    override val f = { (v: EuclideanVector[D]) => if (self.isDefinedAt(v)) Some(self.f(v)) else None }

    override def domain: RealSpace[D] = RealSpace[D]
  }
}

/**
 * Companion objects for Field
 */
object Field {
  def apply[D, A](dom: Domain[D], fun: EuclideanVector[D] => A) = new Field[D, A] {
    override def domain: Domain[D] = dom

    override val f: (EuclideanVector[D] => A) = fun
  }

  def lift[D, A](fl: A => A): Field[D, A] => Field[D, A] = {
    img: Field[D, A] =>
      new Field[D, A] {
        override def apply(v: EuclideanVector[D]) = fl(img.apply(v))

        override val f = img.f

        def domain: Domain[D] = img.domain
      }
  }
}

/** Differntiable field extension. includes a method to differentiate a field. */
trait DifferentiableField[D, A, dA] extends Field[D, A] {
  self =>
  val df: EuclideanVector[D] => dA

  def differentiate: Field[D, dA] = {
    Field(domain, df)
  }
}

/**
 * Companion object for DifferentiableField trait.
 */
object DifferentiableField {
  def apply[D, A, dA](domain: Domain[D], f: EuclideanVector[D] => A, df: EuclideanVector[D] => dA): DifferentiableField[D, A, dA] = {
    val outerdf = df
    val outerf = f
    val outerdomain = domain

    new DifferentiableField[D, A, dA] {
      override val df: EuclideanVector[D] => dA = outerdf
      override val f: EuclideanVector[D] => A = outerf

      override def domain: Domain[D] = outerdomain
    }
  }
}

/**
 * Scalar field where the function is from Euclidean n-dimensional space to the reals.
 *
 * @param domain Domain over which function f is defined.
 * @param f      Function which is implemented as part of the field.
 * @tparam D Dimension of the space of the domain
 */
case class ScalarField[D](domain: Domain[D], f: EuclideanVector[D] => Double) extends Field[D, Double] {

  def +(that: ScalarField[D]): ScalarField[D] = {
    def f(v: EuclideanVector[D]): Double = this.f(v) + that.f(v)

    new ScalarField[D](Domain.intersection[D](domain, that.domain), f)
  }

  def -(that: ScalarField[D]): ScalarField[D] = {
    def f(v: EuclideanVector[D]): Double = this.f(v) - that.f(v)

    new ScalarField[D](Domain.intersection[D](domain, that.domain), f)
  }

  def :*(that: ScalarField[D]): ScalarField[D] = {
    def f(v: EuclideanVector[D]): Double = this.f(v) * that.f(v)

    new ScalarField[D](Domain.intersection[D](domain, that.domain), f)
  }

  def *(s: Double): ScalarField[D] = {
    def f(v: EuclideanVector[D]): Double = this.f(v) * s

    new ScalarField[D](domain, f)
  }

  def compose(t: EuclideanVector[D] => EuclideanVector[D]): ScalarField[D] = {
    val f = this.f.compose(t)
    val newDomain = Domain.fromPredicate[D](v => isDefinedAt(t(v)))
    new ScalarField[D](newDomain, f)
  }
}

/**
 * Vector field where function maps from N-dimensional euclidean space to M-dimensional euclidean space.
 *
 * @param domain Domain over which function f is defined.
 * @param f      Function which is implemented as part of the field.
 * @tparam D Dimension of the space of the domain
 */
case class VectorField[D](domain: Domain[D], f: EuclideanVector[D] => EuclideanVector[D]) extends Field[D, EuclideanVector[D]] {

  def +(that: VectorField[D]): VectorField[D] = {
    def f(v: EuclideanVector[D]): EuclideanVector[D] = this.f(v) + that.f(v)

    new VectorField[D](Domain.intersection[D](domain, that.domain), f)
  }

  def -(that: VectorField[D]): VectorField[D] = {
    def f(v: EuclideanVector[D]): EuclideanVector[D] = this.f(v) - that.f(v)

    new VectorField[D](Domain.intersection[D](domain, that.domain), f)
  }

  def *(s: Double): VectorField[D] = {
    def f(v: EuclideanVector[D]): EuclideanVector[D] = this.f(v) * s

    new VectorField[D](domain, f)
  }

  def compose(t: EuclideanVector[D] => EuclideanVector[D]): VectorField[D] = {
    val f = this.f.compose(t)
    val newDomain = Domain.fromPredicate[D](v => isDefinedAt(t(v)))
    new VectorField[D](newDomain, f)
  }
}

