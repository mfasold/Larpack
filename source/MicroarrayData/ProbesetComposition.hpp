/**
 * @file ProbesetComposition.hpp Defines classes to create probeset compositions from probes. 
 * @author $Author$ 
 * @author Mario Fasold
 * @date $Date$
 */
#ifndef _PROBESETCOMPOSITION_
#define _PROBESETCOMPOSITION_

#include "Probeset.hpp"
#include "PmMmProbe.hpp"
#include "SnpProbe.hpp"
#include "HookcurveAnalyzer.hpp"
#include "DetectKinkPoint.hpp"

namespace larrpack {

/**
 * @class ProbesetCompositor
 *
 * Abstract class that provides some common methods for probeset composition.
 */
class ProbesetCompositor
{
protected:
  /// Creates a list of boundaries based on the distace of probes
  static std::vector<size_t> defineDistanceConstrainedBoundaries(const PmMmProbePtrVector& probes, size_t maxDistance);

  /// Creates a list of boundaries for shuffling based on existing boundaries.
  static std::vector<size_t> defineShufflingBoundaries(const std::vector<size_t>& existingBoundaries, 
                                                       size_t shuffleIndex, size_t maxIndex);
  
  /// Defines a list of boundaries based on the probabilities of the probes to be in the specific domain
  /// of the hookcurve
  static std::vector<size_t> defineProbespecificityBoundaries(const std::valarray<double>& specificityFrequencies, 
                                                              double threshold, size_t movingAverageWindowSize);

  static ProbesetVector createProbesetCompostionFromBoundaries(const PmMmProbePtrVector& probes,
                                                               const std::vector<size_t>& boundaries, 
                                                               size_t optimalProbesPerSet, size_t maxProbesPerSet);
};

/**
 * @class ComposeProbesetsById
 *
 * Groups all probes with the same ProbesetId into one probeset.
 */
class ComposeProbesetsById : public ProbesetCompositor
{
public:
  ProbesetVector operator() (PmMmProbePtrVector&);
};


/**
 * Creates a probeset composition based on following constaints
 *  \li a probeset ideally contains optimalProbesPerSet probes, but not more than maxProbesPerSet
 *  \li the genomic distance between two succint probes within a set does not exceed maxDistance
 *  \li the offset of the probeset is shifted shuffleIndex probes w.r.t. to a composition
 *      satisfied only by the upper two constraints
 *
 * The function creates a list of boundaries, based on the latter two constaints and
 * runs createProbesetCompostionFromBoundaries. The probeset composition is saved internally.
 */
class ComposeProbesetsByShuffling : public ProbesetCompositor
{
public:
  ComposeProbesetsByShuffling(size_t shuffleIndex, int maxDistance, 
                              size_t optimalProbesPerSet, size_t maxProbesPerSet);
  ProbesetVector operator() (const PmMmProbePtrVector&);
private:
  size_t mShuffleIndex;
  int    mMaxDistance;
  size_t mOptimalProbesPerSet;
  size_t mMaxProbesPerSet;
};


/**
 * @class ComposeProbesetsBySpecificityFrequencies
 *
 * Creates a probeset composition based on following constaints
 *  \li a probeset ideally contains optimalProbesPerSet probes, but not more than maxProbesPerSet
 *  \li the genomic distance between two succint probes within a set does not exceed maxDistance
 *  \li no assumed boundary between transcribed and untranscribed region is within a set.
 * 
 * For the last constraint, the routine uses a list containing the "probability of each probe
 * to be transcribed" (We use the "probablity of a probe to bind specific" obtained by 
 * shuffling). The funcion filters those probablities with a moving average and inserts a 
 * boundary every time the moving average crosses the threshold. 
 * 
 * These boundaries are combined with the boundaries from the first two constaints. It then
 * runs createProbesetCompostionFromBoundaries and saves the probeset composition internally.
 *
 */
class ComposeProbesetsBySpecificityFrequencies : public ProbesetCompositor
{
public:
  ComposeProbesetsBySpecificityFrequencies(const std::valarray<double>& specificityFrequencies, 
                                           double threshold, size_t movingAverageWindowSize,
                                           int maxDistance, 
                                           size_t optimalProbesPerSet, size_t maxProbesPerSet);
  ProbesetVector operator() (const PmMmProbePtrVector&);
private:
  const std::valarray<double>& mSpecificityFrequencies;
  double mThreshold;
  size_t mMovingAverageWindowSize;
  size_t mShuffleIndex;
  int    mMaxDistance;
  size_t mOptimalProbesPerSet;
  size_t mMaxProbesPerSet;
};

/**
 * @class ComposeProbesetsByCompleteShufflingMethod
 *
 * Runs the complete method called "probeset shuffling". For each offset 
 * \f$ shuffleIndex = 0 \cdots (optimalProbesPerSet - 1) \f$, a probeset 
 * composition is created by ComposeProbesetsByShuffling. We calculate
 * a hookcurve for each one and thereof the "probablity of a probe to 
 * bind specific". With these probabilities, 
 * ComposeProbesetsBySpecificityFrequencies builds up the final probeset
 * composition.
 *
 */
class ComposeProbesetsByCompleteShufflingMethod : public ProbesetCompositor
{
public:
  ComposeProbesetsByCompleteShufflingMethod(int maxDistance, 
                                            size_t optimalProbesPerSet, 
                                            size_t maxProbesPerSet,
                                            double threshold,
                                            size_t hookcurveMovingAverageWindowSize,
                                            size_t shufflingMovingAverageWindowSize,
                                            DetectKinkPointPtr detectKinkPointFunction);
  ProbesetVector operator() (const PmMmProbePtrVector&);
private:
  int    mMaxDistance;
  size_t mOptimalProbesPerSet;
  size_t mMaxProbesPerSet;

  double mThreshold;
  size_t mHookcurveMovingAverageWindowSize;
  size_t mShufflingMovingAverageWindowSize;
  DetectKinkPointPtr detectKinkPoint;
};


/**
 * Creates probeset compositions for SNP chips based on a genotype call. The genotype
 * call is read in from a tab delimited file as it is provided by affymetrix. For each
 * chip, 3 types of probeset compositions are returned:
 *  \li Containing allelspecific probes
 *  \li Containing cross-allele probes
 *  \li Containing heterocycote probes 
 */
class ComposeProbesetsByGenotypeCall : public ProbesetCompositor
{
public:
  ComposeProbesetsByGenotypeCall(const std::map<std::string,char> genotypeCalls, const Haplotype haplotype);
  ProbesetVector operator() (PmMmProbePtrVector&);
  static SnpProbesetVector composeSnpProbesets(SnpProbePtrVector&, std::map<std::string,char>& genotypeCalls);


  static std::map<std::string,char> readGenotypeCalls(const std::string filename, const std::string chipname);
  static std::map<std::string,char> readGenotypeCallsFromGDAS(const std::string filename, const std::string chipname);
private:
  std::map<std::string,char> mGenotypeCalls;
  Haplotype mHaplotype;
};
 

/**
 * DEBUG CLASS
 * Creates a probeset composition using boundaries from a BED-file
 * (usually known gene annotations)
 */
class ComposeProbesetsFromBedfileSegments : public ProbesetCompositor
{
public:
  ComposeProbesetsFromBedfileSegments(int maxDistance, 
                                      size_t optimalProbesPerSet, 
                                      size_t maxProbesPerSet,
                                      const std::string filename) 
  : mMaxDistance(maxDistance),
    mOptimalProbesPerSet(optimalProbesPerSet), 
    mMaxProbesPerSet(maxProbesPerSet),
    mBedfilename(filename)
  {}
  ProbesetVector operator() (const PmMmProbePtrVector&);
private:
  int    mMaxDistance;
  size_t mOptimalProbesPerSet;
  size_t mMaxProbesPerSet;

  std::string mBedfilename;
};


}
#endif
