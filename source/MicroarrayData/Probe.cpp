/**
 * @file Probe.cpp Defines probe functions.
 * @author $Author$ 
 * @author Mario Fasold
 * @date $Date$
 */
#include <string>
#include <boost/lexical_cast.hpp>

#include "Probe.hpp"

using namespace std;
using namespace larrpack;

/**
 * Constructor
 *
 * @todo describe parameters
 */
Probe::Probe(const std::string& chromosome, const char strand, const std::string& sequence, 
             const int location, const IntensityType pm,
             const std::pair<int,int>& pmPosition,
             const std::string& probesetId) : 
  mChromosome(chromosome), mStrand(strand), mSequence(sequence), mLocation(location), 
  mPm(pm), mPmPosition(pmPosition), mProbesetId(probesetId)
{
}

/**
 * Returns a string representation of the probe.
 *
 * @return string representation of the probe
 */
std::string Probe::toString() 
{
  return mProbesetId + " " + boost::lexical_cast<std::string>(mLocation)
    + " " + mSequence + " " + mStrand + 
    " " + boost::lexical_cast<std::string>(mPmPosition.first) +
    "/" + boost::lexical_cast<std::string>(mPmPosition.second)
    + "\n";
}

/**
 * Returns the location of the probe on the chromsome.
 *
 * @return Location of the probe on the chromsome.
 */
int Probe::getLocation() const 
{
  return mLocation;
}

/**
 * Returns the perfect match intensity.
 *
 * @return Perfect match intensity
 */
IntensityType Probe::getPm() const
{
  return mPm;
}

/**
 * Returns position of the PM spot on the array.
 *
 * @return Position of the PM spot on the array.
 */
pair<int,int> Probe::getPositionPm() const
{
  return mPmPosition;
}

/**
 * Returns the chromosome (or generally any sequence id) the probe is 
 * located on.
 *
 * @return The chromsome (or generally any sequence id).
 */
std::string Probe::getChromosome() const 
{
  return mChromosome;
}

/**
 * Returns the sequence of the PM-probe.
 *
 * @return The sequence of the PM-probe.
 */
std::string Probe::getSequence() const 
{
  return mSequence;
}

/**
 * Returns the strand of the PM-probe.
 *
 * @return The strand of the PM-probe.
 */
char Probe::getStrand() const
{
  return mStrand;
}

/**
 * Returns the id of the corresponding probeset.
 *
 * @return The id of the corresponding probeset.
 *
 */
std::string Probe::getProbesetId() const 
{
  return mProbesetId;
}

/**
 * @debug we use this function to mark specific probes.
 * 
 */
void Probe::setProbesetId(string id) 
{
  mProbesetId = id;
}
