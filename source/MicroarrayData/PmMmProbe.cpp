/**
 * @file PmMmProbe.cpp
 * @author $Author$ 
 * @date $Date$
 */
#include <string>
#include <cmath>
#include <boost/lexical_cast.hpp>

#include "PmMmProbe.hpp"
#include "SequenceUtil.hpp"

using namespace std;
using namespace larrpack;

/**
 * Constructor
 *
 * @todo describe parameters
 */
PmMmProbe::PmMmProbe(const std::string& chromosome, const char strand, const std::string& sequence, 
                     const int location, const IntensityType pm, const IntensityType mm,
                     const std::pair<int,int>& pmPosition, const std::pair<int,int>& mmPosition,
                     const std::string& probesetId) 
  : Probe(chromosome, strand, sequence, location, pm, pmPosition, probesetId),
    mMm(mm), mMmPosition(mmPosition)
{
}

/**
 * Returns the mismatch intensity.
 *
 * @return Perfect match intensity
 */
IntensityType PmMmProbe::getMm() const 
{
  return mMm;
}

/**
 * Returns log(PM) - log(MM). Logarithms to the base 10 are used.
 *
 * @return log(PM) - log(MM)
 */
IntensityType PmMmProbe::getDeltaLogI() const
{
  return log10(getPm()) - log10(getMm());
}

/**
 * Returns log(PM) + log(MM). Logarithms to the base 10 are used.
 *
 * @return log(PM) + log(MM)
 */
IntensityType PmMmProbe::getSumLogI() const
{
  return log10(getPm()) + log10(getMm());
}


/**
 * Returns position of the MM spot on the array.
 *
 * @return Position of the MM spot on the array.
 */
pair<int,int> PmMmProbe::getPositionMm() const
{
  return mMmPosition;
}

/**
 * Returns the sequence of the MM-probe.
 *
 * @return The sequence of the MM-probe.
 */
std::string PmMmProbe::getSequenceMm() const
{
  return sequenceutil::complementMiddlebase(getSequence());
}


/**
* Returns the sensitivity of the PM probe sequence.
*
* @param sensitivityArray the array holding the calculated sensitivities.
*
* @return Y^{PM}_p The sensitivity of the PM sequence.)
*/
double PmMmProbe::getSensitivityPm(std::valarray<double>& sensitivityArray) const
{
  ///todo: remove sequenceLength from here and anywhere else (e.g AnalyzedMicroarray::calculateSequenceSensitivity)
   //      to just one location. (as static constant?)
  size_t sequenceLength = 25;
  string sequence = this->getSequence();
  double sensitivity = 0.0;
  for (size_t posIndex = 0; posIndex < sequence.length(); ++posIndex){
    sensitivity += sensitivityArray[sequenceutil::baseToIndex(sequence[posIndex])*sequenceLength + posIndex];
  }
  return sensitivity;
}


/**
* Returns the sensitivity of the MM probe sequence.
*
* @param sensitivityArray the array holding the calculated sensitivities.
*
* @return Y^{MM}_p The sensitivity of the MM sequence.)
*/
double PmMmProbe::getSensitivityMm(std::valarray<double>& sensitivityArray) const
{
  /// @todo: remove sequenceLength from here and anywhere else (e.g AnalyzedMicroarray::calculateSequenceSensitivity)
   //      to just one location. (as static constant?)
  size_t sequenceLength = 25;
  string sequence = this->getSequence(); //ask mario if this is really copied and not a reference (damn c++)
  sequence[(int) sequenceLength/2] = sequenceutil::complement(sequence[(int) sequenceLength/2]);
  double sensitivity = 0.0;
  for (size_t posIndex = 0; posIndex < sequence.length(); ++posIndex){
    sensitivity += sensitivityArray[sequenceutil::baseToIndex(sequence[posIndex])*sequenceLength + posIndex];
  }
  return sensitivity;
}


/**
 * Extracts all PmMmProbes from a list of probes.
 *
 * @param probes The probes.
 *
 * @return A list of all PmMmProbes within the probes.
 *
 * @note Copying big vectors may be time-consuming...
 */
PmMmProbePtrVector larrpack::extractPmMmProbesFromProbelist(const ProbePtrVector& probes)
{
  PmMmProbePtrVector pmMmProbes;
  for (ProbePtrVector::const_iterator pi = probes.begin(); pi != probes.end(); ++pi) {
    // For later: Use dynamic cast in case of virtual classes
    // PmMmProbePtr pmMmProbe = boost::shared_dynamic_cast<PmMmProbe>(*pi);
    PmMmProbePtr pmMmProbe = boost::shared_static_cast<PmMmProbe>(*pi);
    if (pmMmProbe) {
        pmMmProbes.push_back(pmMmProbe);
    }
  }
  return pmMmProbes;
}
