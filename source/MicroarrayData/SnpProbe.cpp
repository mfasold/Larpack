/**
 * @file SnpProbe.cpp
 * @author $Author: kasper $ 
 * @date $Date: 2008-04-17 16:12:26 +0200 (Thu, 17 Apr 2008) $
 */
#include <string>
#include <cmath>
#include <boost/lexical_cast.hpp>

#include "SnpProbe.hpp"
#include "SequenceUtil.hpp"

using namespace std;
using namespace larrpack;

/**
 * Constructor
 *
 */
SnpProbe::SnpProbe(const std::string& chromosome, const char strand, const std::string& sequence, 
                   const int location, const IntensityType pm, const std::pair<int,int>& pmPosition,
                   const std::string& probesetId, char allele)
  : Probe(chromosome, strand, sequence, location, pm, pmPosition, probesetId), 
    mAllele(allele)
{
}

/**
 * Returns a string representation of the probe.
 *
 * @return string representation of the probe
 */
std::string SnpProbe::toString() 
{
  return mProbesetId + " " +
    boost::lexical_cast<std::string>(mPmPosition.first) + "/" + 
    boost::lexical_cast<std::string>(mPmPosition.second) + " " + 
    boost::lexical_cast<std::string>(mLocation) + " " + 
    + " " + mSequence + " " + mStrand  + " " + 
    mAllele + "\n";
}



/**
 * Extracts all SnpProbes from a list of probes.
 *
 * @param probes The probes.
 *
 * @return A list of all SnpProbes within the probes.
 *
 * @note Copying big vectors may be time-consuming...
 */
SnpProbePtrVector larrpack::extractSnpProbesFromProbelist(const ProbePtrVector& probes)
{
  SnpProbePtrVector snpProbes;
  for (ProbePtrVector::const_iterator pi = probes.begin(); pi != probes.end(); ++pi) {
    // For later: Use dynamic cast in case of virtual classes
    // PmMmProbePtr pmMmProbe = boost::shared_dynamic_cast<PmMmProbe>(*pi);
    SnpProbePtr snpProbe = boost::shared_static_cast<SnpProbe>(*pi); 
    if (snpProbe) {
      snpProbes.push_back(snpProbe);
    }
  }
  return snpProbes;
}
