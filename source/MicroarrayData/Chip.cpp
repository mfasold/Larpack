/**
 * @file Chip.cpp Provides a data-structure to hold all chip-specific data
 * @author $Author: jan $ 
 * @date $Date: 2008-06-17 12:53:14 +0200 (Tue, 17 Jun 2008) $
 */
#include <Chip.hpp>

using namespace std;

namespace larrpack
{
/**
 * Constructor that creates the internal probesets using the
 * ProbesetCompositionFunction. 
 *
 * @param probes List of PM/MM-probes
 * @param composeProbesetFunction Function object to create probesets composition.
 */
Chip::Chip(const PmMmProbePtrVector& probes, 
           ProbesetVector probesets)
  : mProbes(probes), mProbesets(probesets)
{}

/**
 * Constructor that extracts all PmMmProbes from a list of probes.
 *
 * @param probes List of probes
 */
Chip::Chip(const ProbePtrVector& probes)
{
	mProbes = extractPmMmProbesFromProbelist(probes);
}

/**
 * Constructor that extracts all PmMmProbes from a list of probes
 * and creates the internal probesets using the
 * ProbesetCompositionFunction.
 *
 * @param probes List of probes
 * @param composeProbesetFunction Function object to create probesets composition.
 */
Chip::Chip(const ProbePtrVector& probes, 
           ProbesetVector probesets)
  : mProbes(extractPmMmProbesFromProbelist(probes)), mProbesets(probesets)
{}

/**
 * Removes some of the probesets using a given filter function.
 *
 * @param filterProbeset Function predicate.
 */
void Chip::filterProbesets(ProbesetPredicate filterProbeset)
{
  mProbesets.erase(remove_if(mProbesets.begin(), mProbesets.end(), not1(filterProbeset)), 
                   mProbesets.end());
}

/**
 * Returns all probes
 *
 * @return A shared pointer to a vector of probes. 
 */
const PmMmProbePtrVector& Chip::getProbes() const
{
	return mProbes;
}

/**
 * Returns all probesets
 *
 * @return Array of probesets.
 */
const vector<Probeset>& Chip::getProbesets() const
{
  return mProbesets;
}

/**
 * Returns all probesets
 *
 * @return The number of probesets of the chip
 */
const size_t Chip::getProbesetCount() const
{
	return mProbesets.size();
}

}
