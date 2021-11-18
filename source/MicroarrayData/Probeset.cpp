/**
 * @file Probeset.cpp Defines a functions for aggregated probes (PmMmProbes).
 * @author $Author$
 * @author Mario Fasold
 * @date $Date$
 */ 
#include <vector>

#include "Probe.hpp"
#include "Probeset.hpp"
#include "MathUtil.hpp"
#include "SequenceUtil.hpp"

#include "AveragingProcedure.hpp"


using namespace std;
using namespace larrpack;
using namespace sequenceutil;

/// Set intensity function to default value
IntensitiesFunction larrpack::computeAverageIntensity = Mean<IntensityType>();

/**
 * Constructor which initializes the probeset with 
 * some probes. Also initializes the corrected probe
 * intensities with the original PM/MM intensities.
 *
 * @param probes List of probes
 *
 */
Probeset::Probeset(const PmMmProbePtrVector& probes) 
  :  mProbes(probes), mProbesetId("none")
{
}

/**
 * Constructor which initializes the probeset with a slice
 * of a probe-vector (iterators).Also initializes the corrected probe 
 * intensities with the original PM/MM intensities.
 *
 * @param probesetBegin Iterator pointing to the begin of the slice of probes,
 *                      which make up this probeset
 * @param probesetEnd   Iterator pointing to element after the slice of probes
 */
Probeset::Probeset(PmMmProbePtrVector::const_iterator probesetBegin, 
                   PmMmProbePtrVector::const_iterator probesetEnd,
                   std::string probesetId) 
  : 
    mProbesetId(probesetId)
{
  copy(probesetBegin, probesetEnd, back_inserter(mProbes));
//   size_t probeIndex = 0;
//   for (PmMmProbePtrVector::const_iterator probeIt = probesetBegin; probeIt != probesetEnd; ++probeIt, ++probeIndex) {
//     // Save all probe pointers in local member variable
//     mProbes.push_back(*probeIt);
//   }
  //cout << mProbesetId << "  " << mCurrentCorrectedPmLogIs.size() << endl;
}

/**
 * Returns the number of elements within the set.
 *
 * @return The number of elements within the set.
 */
size_t Probeset::getSize() const
{
  return mProbes.size();
}


/**
 * Returns the i-th probe contained in this set.
 *
 * @note Using this function is neither good for speed (probe copy) 
 *       nor for encapsulation...
 * @note Perhaps return a pointer for speed?
 *
 * @return The probe i within the set.
 */
PmMmProbe Probeset::getProbe(const size_t probeIndex) const
{
//	cout << mProbes.size() << endl;
//	cout << (*mProbes[probeIndex]).getLocation() <<endl; 
  return *mProbes[probeIndex];
}

/**
 * Returns a pointer to the i-th probe contained in this set.
 *
 * @note Using this function not good for encapsulation...
 * @todo Remove this function.
 *
 * @return The probe i within the set.
 */
PmMmProbePtr Probeset::getProbePtr(const size_t probeIndex) const
{
  return mProbes[probeIndex];
}

/**
 * Returns the probesetId of the probeset
 * 
 * if there is no id (like in tiling arrays)
 * the id is set to "no_id"
 * 
 * @return The probeset id
 */
std::string Probeset::getProbesetId() const
{
  return mProbesetId;
}


/**
 * @debug
 * 
 */
void Probeset::setProbesetId(string id)
{
	mProbesetId = id;
}

/**
 * Returns an array with PM or MM intensity values of each probe of the probeset.
 *
 * @param pmOrMm pm or mm probe type
 *
 * @return An array with PM or MM intensity values of the probeset.
 */
IntensityArray Probeset::getProbeIntensities(const IntensityMode pmOrMm) const {
  IntensityArray intensities(getSize());
  if (pmOrMm == kPmI) {
    transform(beginProbe(), endProbe(), &intensities[0], boost::mem_fn(&PmMmProbe::getPm));
  } 
  else if (pmOrMm == kMmI) {
    transform(beginProbe(), endProbe(), &intensities[0], boost::mem_fn(&PmMmProbe::getMm));
  }
  else {
    assert(NULL == "Attempt to return other than PM or MM probe type");
  }
  return intensities;
}

/**
 * Returns an array with PM or MM sequences
 *
 * @param pmOrMm pm or mm probe type
 *
 * @return An array with PM or MM sequences of the probeset.
 */
//std::vector<std::string> Probeset::getProbeSequences(const IntensityMode pmOrMm) const
std::vector<std::string> Probeset::getProbeSequences(const PmMmProbeSequenceFunction& getProbeSequence) const
{
  std::vector<std::string> sequences; // (getSize());
  for (PmMmProbePtrVector::const_iterator pi = mProbes.begin(); pi != mProbes.end(); ++pi) {  
    sequences.push_back(getProbeSequence(*(*pi)));
  }
  // I'd like to know I can rewrite this using transform:
  // transform(beginProbe(), endProbe(), &sequences[0], getProbeSequence);
  return sequences;
}



/**
 * Initializes an array of probeset intensities with the same size of
 * the probesets.
 *
 * @param probesets List of probesets
 *
 * @return Pointer to a list of intensity arrays equally sized as the probesets.
 */
ProbesetIntensitiesArrayPtr Probeset::initializeIntensitiesArray(const ProbesetVector& probesets)
{
   ProbesetIntensitiesArrayPtr intensitiesArrayPtr(new ProbesetIntensitiesArray());
   for (ProbesetVector::const_iterator ps = probesets.begin(); ps != probesets.end(); ++ps) {;
     intensitiesArrayPtr->push_back(IntensityArray(ps->getSize()));
   }
   return intensitiesArrayPtr;
}

/**
 * Initializes an array of probeset intensities with actual probe intensities 
 * with the same size of the probesets.
 *
 * @param probesets List of probesets
 * @param pmOrMm pm or mm probe type

 *
 * @return Pointer to a list of intensity arrays equally sized as the probesets.
 */
ProbesetIntensitiesArrayPtr 
Probeset::initializeIntensitiesArray(const ProbesetVector& probesets, const IntensityMode pmOrMm)
{
  ProbesetIntensitiesArrayPtr intensitiesArrayPtr(new ProbesetIntensitiesArray());
  for (ProbesetVector::const_iterator ps = probesets.begin(); ps != probesets.end(); ++ps) {
    intensitiesArrayPtr->push_back(log10(ps->getProbeIntensities(pmOrMm)));
  }
  return intensitiesArrayPtr;
}



/**
 * Constructor which initializes the probeset with a slice
 * of a probe-vector (iterators).Also initializes the corrected probe 
 * intensities with the original PM/MM intensities.
 *
 * @pre The probes have to be in the order given by
 * ProbesetIdAndAlleleAndInterrogationIsSmaller!
 *
 * @param probesetBegin Iterator pointing to the begin of the slice of probes,
 *                      which make up this probeset
 * @param probesetEnd   Iterator pointing to element after the slice of probes
 */
SnpProbeset::SnpProbeset(SnpProbePtrVector::const_iterator probesetBegin, 
                         SnpProbePtrVector::const_iterator probesetEnd,
                         std::string probesetId, const char call) 
  : 
  mProbesetId(probesetId), mCall(call)
{

  // cout << "Probeset " << probesetId << " (" << call << ")." << endl;

  // Define the probe position within the set where second allele begins
  // (probes are sorted and (check this) have always same number of Allele probes!!)
  SnpProbePtrVector::const_iterator currentSecondAlleleBegin = 
    probesetBegin + (probesetEnd - probesetBegin)/2 ;

  // Assert if currentSecondAlleleBegin really hits an allele boundary
  assert((*(currentSecondAlleleBegin-1))->getAllele() != (*currentSecondAlleleBegin)->getAllele());
  // char firstAllele = (*(currentSecondAlleleBegin - 1))->getAllele();
  // char secondAllele = (*currentSecondAlleleBegin)->getAllele();
  // assert(firstAllele != secondAllele);

  // Build allele and cross probes according to call
  if (call == 'H') { // heterozygous case -> only probes of type "allele", none "cross"
    copy(probesetBegin, probesetEnd, back_inserter(mProbes));    
  }
  else if (call == 'A') { // first Part of probeset; use lexocographic sorting of allele here
    copy(probesetBegin, currentSecondAlleleBegin, back_inserter(mProbes));
    copy(currentSecondAlleleBegin, probesetEnd, back_inserter(mProbesCross));
  }
  else if (call == 'B') { // second Part of probeset; use lexocographic sorting of allele here
    copy(probesetBegin, currentSecondAlleleBegin, back_inserter(mProbesCross));
    copy(currentSecondAlleleBegin, probesetEnd, back_inserter(mProbes));
  }
  else {   
    cout << "Invalid call of probeset " << probesetId <<
      " (" << call << "). Exiting." << endl;
    exit(0);
  }
}

/// Plain copy of Probeset::initializeIntensitiesArray
ProbesetIntensitiesArrayPtr SnpProbeset::initializeIntensitiesArray(const SnpProbesetVector& probesets)
{
  ProbesetIntensitiesArrayPtr intensitiesArrayPtr(new ProbesetIntensitiesArray());
  for (SnpProbesetVector::const_iterator ps = probesets.begin(); ps != probesets.end(); ++ps) {;
    intensitiesArrayPtr->push_back(IntensityArray(ps->getSize()));
  }
  return intensitiesArrayPtr;
}

SnpProbePtr SnpProbeset::getProbePtr(const size_t probeIndex) const
{
  return mProbes[probeIndex];
}

size_t SnpProbeset::getSize() const
{
  return mProbes.size();
}

std::vector<std::string> SnpProbeset::getProbeSequences() const
{
  std::vector<std::string> sequences; // (getSize());
  for (SnpProbePtrVector::const_iterator pi = mProbes.begin(); pi != mProbes.end(); ++pi) {  
    sequences.push_back((*pi)->getSequence());
    //sequences.push_back(getProbeSequence(*(*pi)));
  }
  // I'd like to know I can rewrite this using transform:
  // transform(beginProbe(), endProbe(), &sequences[0], getProbeSequence);
  return sequences;
}
