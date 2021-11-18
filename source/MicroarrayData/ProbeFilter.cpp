/**
 * @file ProbeFilter.cpp Defines filters for (PmMm-)Probes.
 * @author $Author$ 
 * @author Mario Fasold
 * @date $Date$
 */
#include "ProbeFilter.hpp"
#include "PmMmProbe.hpp"

// #include <iostream.h>

using namespace larrpack;


/**
 * Returns true iff probes PM intensity is smaller than threshold. If 
 * the probe is a PmMmProbe, it returns true iff PM or MM is smaller 
 * than threshold.
 *
 * @param probe The probe to test.
 */
bool ProbeintensityIsLessThan::operator() (const ProbePtr& probe) const
{ 
  PmMmProbePtr pmMmProbe = boost::shared_static_cast<PmMmProbe>(probe);
  if (pmMmProbe) {
    return pmMmProbe->getPm() < threshold || pmMmProbe->getMm() < threshold;
  }
  else {
    return probe->getPm() < threshold;
  }
}


/**
 * Returns true if the probes probeset_id is in the given array. 
 *
 * @param probe The probe to test.
 */
bool IsInProbesetFilter::operator() (const ProbePtr& probePtr) const
{
  std::vector<std::string>::const_iterator idPtr = find(validProbesetIds.begin(), validProbesetIds.end(), 
                                                        probePtr->getProbesetId());
  if (idPtr == validProbesetIds.end()) {
    return false;
  }
  return true;
}


/**
 * Returns true iff probe sequence contains a substring at a particular
 * position.
 *
 * @param probe The probe to test.
 */
bool ProbesequenceHasSubstring::operator() (const ProbePtr& probe) const
{ 
  return probe->getSequence().substr(mIndex, mSubstring.size()) == mSubstring;
}


/**
 * Returns true iff the probes chromosome equals a given chromsome.
 *
 * @param probe The probe to test.
 */
bool ProbeChromosomeEquals::operator() (const ProbePtr& probe) const
{ 
  return probe->getChromosome() == mChromosome;
}


/**
 * Returns true iff first probe's probeset id is lexicographical
 * smaller than that of the second probe.
 *
 * @param probeLhs Left-hand-side probe object of the comparison.
 * @param probeRhs Richt-hand-side probe object of the comparison.
 *
 */
bool ProbesetIdIsSmaller::operator() (const ProbePtr& probeLhs, const ProbePtr& probeRhs) const
{
  if (probeLhs->getProbesetId() == probeRhs->getProbesetId()) {
    return probeLhs->getLocation() < probeRhs->getLocation();
  }
  return probeLhs->getProbesetId() < probeRhs->getProbesetId();
}

/**
 * Returns true iff (in that order)
 *  1) first probe's probeset id is lexicographical
 *     smaller than that of the second probe.
 *  2) first probe's allele is lexicographical smaller
 *  3) first probe's interrogation position is smaller
 *  4) first probe's strand (f,r) is lexicographical smaller
 *
 *
 * @param probeLhs Left-hand-side probe object of the comparison.
 * @param probeRhs Richt-hand-side probe object of the comparison.
 *
 */
bool ProbesetIdAndAlleleAndInterrogationIsSmaller::operator() (const SnpProbePtr& probeLhs, const SnpProbePtr& probeRhs) const
{
  if (probeLhs->getProbesetId() == probeRhs->getProbesetId()) {
    if (probeLhs->getAllele() == probeRhs->getAllele()) {
      if (probeLhs->getLocation() == probeRhs->getLocation()) {
        return probeLhs->getStrand() < probeRhs->getStrand();      
      }
      return probeLhs->getLocation() < probeRhs->getLocation();      
    }
    return probeLhs->getAllele() < probeRhs->getAllele();
  }
  return probeLhs->getProbesetId() < probeRhs->getProbesetId();
}


