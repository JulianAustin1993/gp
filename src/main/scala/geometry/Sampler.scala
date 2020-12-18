/**
 * Implementing a sampler over a domain.
 * Modified from:
 * https://github.com/unibas-gravis/scalismo/blob/master/src/main/scala/scalismo/numerics/Sampler.scala
 */
package geometry

import breeze.stats.distributions.{RandBasis, Uniform}

/**
 * General Sampler trait
 *
 * @tparam D Dimension of the space.
 */
trait Sampler[D] {

  val numberOfPoints: Int

  val volumeOfSampleRegion: Double

  def sample(): IndexedSeq[(EuclideanVector[D], Double)]

}

/**
 * Sampler for a box domain.
 *
 * @param domain      BoxDomain to be sampled over.
 * @param n           number of samples to take
 * @param NDSpace$D$0 Space over which the domain is defined.
 * @param rand        Implicit RandBasis
 * @tparam D Dimension of the space.
 */
case class UniformBoxDomainSampler[D: NDSpace](domain: BoxDomain[D], n: Int)(implicit rand: RandBasis) extends Sampler[D] {
  override val numberOfPoints: Int = n

  override val volumeOfSampleRegion: Double = domain.volume
  val p: Double = 1.0 / domain.volume

  override def sample(): IndexedSeq[(EuclideanVector[D], Double)] = {
    val ndSpace = implicitly[NDSpace[D]]
    val randGens = for (i <- 0 until ndSpace.dimensionality) yield {
      Uniform(domain.origin(i), domain.oppositeCorner(i))(rand)
    }
    for (_ <- 0 until numberOfPoints) yield (EuclideanVector.apply[D](randGens.map(r => r.draw()).toArray), p)
  }
}

/**
 * Fixed Sampler for a box domain so that multiple calls give same sampling points.
 *
 * @param domain      BoxDomain to be sampled over.
 * @param n           number of samples to take
 * @param NDSpace$D$0 Space over which the domain is defined.
 * @param rand        Implicit RandBasis
 * @tparam D Dimension of the space.
 */
case class FixedPointsUniformBoxDomainSampler[D: NDSpace](domain: BoxDomain[D], n: Int)(implicit rand: RandBasis) extends Sampler[D] {
  override val numberOfPoints: Int = n
  override val volumeOfSampleRegion: Double = domain.volume
  private val samplePoints = UniformBoxDomainSampler(domain, n).sample()

  override def sample(): IndexedSeq[(EuclideanVector[D], Double)] = samplePoints
}
