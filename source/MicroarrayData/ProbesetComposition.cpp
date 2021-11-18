/**
 * @file ProbesetComposition.cpp Defines functions to create probeset compositions from probes. 
 * @author $Author$
 * @author Mario Fasold
 * @date $Date$
 */ 

#include <string>
#include <algorithm>
#include <iterator>
#include <boost/lexical_cast.hpp>
#include <boost/bind.hpp>
using boost::lexical_cast;

#include "ProbesetComposition.hpp"
#include "ProbeFilter.hpp"
#include "MathUtil.hpp"
#include "StlUtil.hpp"
#include "StringUtil.hpp"
#include "HookcurveAnalyzer.hpp"
#include "SnpProbe.hpp"
#include "ValarrayUtil.hpp"


using namespace std;
using namespace larrpack;


/**
 * Creates a list of boundaries based on the distace of probes.
 *
 * A boundary is an index i within the probe list, meaning that 
 * probe i and probe i+1 are not allowed to be within the same probeset. This function
 * returns a sorted list of boundaries based on the constraint that the genomic distance between 
 * two succint probes i and i+1 exceeds maxDistance. If both probes are on different 
 * chromosomes, a boudary is added as well.
 *
 * @note I'm not yet sure if these functions are too general to be included in this class,
 * therefore they are static.
 *
 * @param probes List of Probes.
 * @param maxDistance Maximum distance between succint probes within a probeset.
 *
 * @return List of boundaries 
 */
std::vector<size_t> ProbesetCompositor::defineDistanceConstrainedBoundaries(const PmMmProbePtrVector& probes, size_t maxDistance)
{
  std::vector<size_t> boundaries;
  // For all probes i in 1..n-1  
  for(size_t i = 1; i < probes.size(); ++i) {
    // If genomic distance between (i-1) and i greater than maxDistance and
    // both are on the same chromosome
    if (((size_t)(probes[i]->getLocation() - probes[i-1]->getLocation()) > maxDistance) ||
        (probes[i]->getChromosome() != probes[i-1]->getChromosome())) {
      // Add that index i to boundaries 
      boundaries.push_back(i-1);
    }
  }
  return boundaries;
}

/**
 * Returns a sorted list of boundaries to do shuffling based on previously created list of
 * boundaries. This is posible since shuffling, in principle, just introduces new
 * boundaries some indexes after existing boundaries. This implies that defineShufflingBoundaries
 * must be run after all other boundary-calculating functions.
 *
 * @param existingBoundaries List of existing boundaries.
 * @param shuffleIndex Offest of the probes w.r.t. non shuffled probesets.
 * @param maxIndex Maximal index allowed.All shuffling boundaries are less than this one. 
 * 
 * @return List of boundaries 
 */
std::vector<size_t> ProbesetCompositor::defineShufflingBoundaries(const std::vector<size_t>& existingBoundaries, 
                                                                  size_t shuffleIndex, size_t maxIndex)
{
  vector<size_t> shufflingeBoundaries;
  
  //  Mind the first shuffling index (there is no boundary at position -1)
  if (existingBoundaries[0] > shuffleIndex) {
    shufflingeBoundaries.push_back(shuffleIndex - 1);
  }

  //  Insert one new b. for each existing boundary, in case the next boundary does not come before
  for (size_t i = 0; i < existingBoundaries.size() - 1; ++i) {
    if (existingBoundaries[i+1] - existingBoundaries[i] > shuffleIndex) {
      shufflingeBoundaries.push_back(existingBoundaries[i] + shuffleIndex);
    }
  }

  // If last shuffling Index would be in tha valid range of indexes, insert it
  size_t lastShufflingBoundary = existingBoundaries.back() + shuffleIndex;
  if (lastShufflingBoundary < maxIndex) {
    shufflingeBoundaries.push_back(lastShufflingBoundary);
  }

  return shufflingeBoundaries;
}

/**
 * Returns a sorted list of boundaries to do shuffling based on the probabilities 
 * of the probes to be in the specific domain of the hookcurve. The routine 
 * calculates a moving average over the specificityFrequency and inserts
 * a boundary anytime the averages "intersect" with the threshold value. 
 *
 * @param specificityFrequencies Probability of probe to be in specific domain (1 ~ specific, 0 unspecific)
 * @param threshold The specificity frequency value where to add a boundary.
 * @param movingAverageWindowSize Number of values to average over.
 * 
 * @return List of boundaries 
 */
std::vector<size_t> ProbesetCompositor::defineProbespecificityBoundaries(const std::valarray<double>& specificityFrequencies, 
                                                                         double threshold, size_t movingAverageWindowSize)
{
  // Calculate moving average
  valarray<double> movingAverages = MathUtil::calculateGeneralizedMovingAverage(specificityFrequencies, 
                                                                                movingAverageWindowSize);  
  
  // Loop over averages and insert boundary when threshold crossed
  vector<size_t> boundaries;
  for(size_t i = 0; i < movingAverages.size() - 1; ++i) {
    if ( (movingAverages[i] - threshold <= 0  && movingAverages[i+1] - threshold > 0) ||
         (movingAverages[i] - threshold >= 0  && movingAverages[i+1] - threshold < 0) ) {
      boundaries.push_back(i);      
    }
  }
  
  return boundaries;
}

/**
 * Creates a probeset composition based on following constraints
 *  - a probeset ideally contains optimalProbesPerSet probes, but not more than maxProbesPerSet
 *  - if i in boundaries, no probeset will contain both probes i and i+1
 *
 * The routine needs the last probe index to be containd in the boundaries.
 * The probeset composition is saved as a member variable.
 *
 * @pre The boundaries list must not contain duplicates.
 *
 * @param probes List of PM/MM probes to be composed into probesets.
 * @param boundaries Ascending sorted list of boundaries (ref)
 * @param optimalProbesPerSet Number of probes optimally contained in a probeset.
 * @param maxProbesPerSet Number of probes maximal contained in a probeset.
 *
 * @ref defineDistanceConstrainedBoundaries for a definition of "boundaries".
 */
ProbesetVector ProbesetCompositor::createProbesetCompostionFromBoundaries(const PmMmProbePtrVector& probes,
                                                                          const vector<size_t>& boundaries, 
                                                                          size_t optimalProbesPerSet, size_t maxProbesPerSet)
{
  ProbesetVector probesets;

  // First probeset begins at probe 0
  size_t currentSetsFirstProbeIndex = 0;

  // Begin with first boundary
  vector<size_t>::const_iterator nextBoundary = boundaries.begin();

  // Loop as long as the current probeset start is valid
  while (currentSetsFirstProbeIndex < probes.size()) {
    //  Set number of probes to add to optimalProbesPerSet
    size_t probesInCurrentSet = optimalProbesPerSet;

    // If next boundary is less than maxProbesPerSet far away
    if ( *nextBoundary < currentSetsFirstProbeIndex + maxProbesPerSet) {

      // Add only probes until that boundary (set number)
      probesInCurrentSet = *nextBoundary - currentSetsFirstProbeIndex + 1;

      // Move to next boundary
      ++nextBoundary;
    }

    //    cout << currentSetsFirstProbeIndex + maxProbesPerSet << "\t" << probesInCurrentSet << "\t" << *nextBoundary <<endl;

    // Add the number of probes just defined to the last probeset
    probesets.push_back( Probeset(larrpack::getVectorSlice(probes, 
                                                           currentSetsFirstProbeIndex,
                                                           size_t(currentSetsFirstProbeIndex + probesInCurrentSet - 1))));

    // Increase currentSetsFirstProbeIndex
    currentSetsFirstProbeIndex += probesInCurrentSet;
  }  

  return probesets;
}



/**
 * Creates probeset composition based on the probeset id property of the probes.
 *
 * @param probes List of probes
 * 
 * @return List of probesets.
 */ 
ProbesetVector ComposeProbesetsById::operator() (PmMmProbePtrVector& probes) 
{
  ProbesetVector probesets;

  // Sort probes to make sure one set's probes are adjacent to each other
  // @note Order in Probesequencefile is not always in desired order!
  sort(probes.begin(), probes.end(), ProbesetIdIsSmaller());


  // for_each(probes.begin(), probes.end(), test_out);
  
  // Stores the id of the current probeset
  string currentProbesetId = probes.front()->getProbesetId();

  // and store iterator where that probeset begun
  PmMmProbePtrVector::const_iterator currentProbesetBegin = probes.begin();

  
  // Iterate over the list of probes
  for(PmMmProbePtrVector::const_iterator currentProbe = probes.begin(); 
      currentProbe != probes.end(); ++currentProbe) {
	  
	// @bug Can we assume that probes of one probeset are always listed successively ?
	  
    // If probeset id of the new probe equals the current probeset, goto next probe
    if ((*currentProbe)->getProbesetId() == currentProbesetId) {
      continue;
    }
    
    // Otherwise create a new probeset with all the probes
    probesets.push_back( Probeset(currentProbesetBegin, currentProbe, currentProbesetId));

    // and update iterator and id
    currentProbesetBegin = currentProbe;
    currentProbesetId = (*currentProbe)->getProbesetId();
  }

  // Create loop leftover probeset
  if (probes.end() - currentProbesetBegin > 1) {
    probesets.push_back( Probeset(currentProbesetBegin, probes.end(), currentProbesetId));    
  }  
  return probesets;
}

/**
 * Constructor to set parameters for probeset shuffling.
 *
 * @see ComposeProbesetsByShuffling
 *
 * @param shuffleIndex Offest of the probes w.r.t. non shuffled probesets.
 * @param maxDistance Maximum distance between succint probes within a set.
 * @param optimalProbesPerSet Number of probes optimally contained in a probeset.
 * @param maxProbesPerSet Number of probes maximal contained in a probeset.
 *
 */
ComposeProbesetsByShuffling::ComposeProbesetsByShuffling(size_t shuffleIndex, int maxDistance, 
                                                         size_t optimalProbesPerSet, size_t maxProbesPerSet)
  :   mShuffleIndex(shuffleIndex), mMaxDistance(maxDistance),
      mOptimalProbesPerSet(optimalProbesPerSet), mMaxProbesPerSet(maxProbesPerSet)
{
}

/**
 * Creates the probeset composition.
 *
 * @see ComposeProbesetsByShuffling
 */ 
ProbesetVector ComposeProbesetsByShuffling::operator() (const PmMmProbePtrVector& probes) 
{
  // Test preconditions
  assert(probes.size() > 0);
  assert(mMaxDistance > 1);
  assert(mOptimalProbesPerSet >= 1);
  assert(mMaxProbesPerSet >= 1);

  // Create boundaries for probes that are too distant
  vector<size_t> distanceBoundaries = ProbesetCompositor::defineDistanceConstrainedBoundaries(probes, mMaxDistance);

  // Changing shuffle index can be realised by inserting an additional boundary
  // after every existing one -> create a list of those additional boundaries
  vector<size_t> shufflingeBoundaries;

  // Do shuffling only if the suffle index makes sense
  if (mShuffleIndex < mOptimalProbesPerSet && mShuffleIndex > 0) {
    shufflingeBoundaries = defineShufflingBoundaries(distanceBoundaries, mShuffleIndex, probes.size() - 1);
  }

//   cout << endl;
//   for_each(shufflingeBoundaries.begin(), shufflingeBoundaries.end(), output__ );
//   cout << endl;
  
  //  Merge all boundaries to one list
  vector<size_t> allBoundaries;
  merge(distanceBoundaries.begin(), distanceBoundaries.end(), 
        shufflingeBoundaries.begin(), shufflingeBoundaries.end(), 
        back_inserter(allBoundaries));

  // The next routine needs the last probe index as boundary
  allBoundaries.push_back(probes.size() - 1);

  return createProbesetCompostionFromBoundaries(probes, allBoundaries, mOptimalProbesPerSet, mMaxProbesPerSet);
}



/**
 * Constructor to set parameters for probeset creation.
 *
 * @param specificityFrequencies Probability of probe to be in specific domain (1 ~ specific, 0 unspecific)
 * @param threshold The specificity frequency value where to add a boundary.
 * @param movingAverageWindowSize Number of specificity frequencies to average over.
 * @param maxDistance Maximum distance between succint probes within a set.
 * @param optimalProbesPerSet Number of probes optimally contained in a probeset.
 * @param maxProbesPerSet Number of probes maximal contained in a probeset.
 *
 * @see ComposeProbesetsBySpecificityFrequencies
 */
ComposeProbesetsBySpecificityFrequencies::ComposeProbesetsBySpecificityFrequencies(const std::valarray<double>& specificityFrequencies, 
                                                                                   double threshold, size_t movingAverageWindowSize,
                                                                                   int maxDistance, 
                                                                                   size_t optimalProbesPerSet, size_t maxProbesPerSet)
  :   mSpecificityFrequencies(specificityFrequencies), 
      mThreshold(threshold), 
      mMovingAverageWindowSize(movingAverageWindowSize),
      mMaxDistance(maxDistance),
      mOptimalProbesPerSet(optimalProbesPerSet), mMaxProbesPerSet(maxProbesPerSet)
{}

/**
 * Creates the probeset composition.
 *
 * @todo Try to combine this routine with ComposeProbesetsByShuffling::operator(), since
 *       80% code is similar.
 *
 * @see ComposeProbesetsBySpecificityFrequencies
 */ 
ProbesetVector ComposeProbesetsBySpecificityFrequencies::operator() (const PmMmProbePtrVector& probes) 
{
  // Test preconditions
  assert(probes.size() > 0);
  assert(mMaxDistance > 1);
  assert(mOptimalProbesPerSet >= 1);
  assert(mMaxProbesPerSet >= 1);

  // Create boundaries for probes that are too distant
  vector<size_t> distanceBoundaries = defineDistanceConstrainedBoundaries(probes, mMaxDistance);

  // Changing shuffle index can be realised by inserting an additional boundary
  vector<size_t> probespecificityBoundaries = defineProbespecificityBoundaries(mSpecificityFrequencies, mThreshold, mMovingAverageWindowSize);

  // @debug: export boundaries into a file
  ofstream logFile;
  logFile.open("ShufflingBoundaries.wig");
  logFile << "track type=wiggle_0 name=\"Bed Format\" description=\"ShufflingBoundaries\" visibility=full" << endl;
  size_t currentBound = 0;
  while (currentBound < probespecificityBoundaries.size() - 1) { // while at least two entries left
    const PmMmProbe& p1 = *probes[probespecificityBoundaries[currentBound]];
    const PmMmProbe& p2 = *probes[probespecificityBoundaries[currentBound + 1]];
    if (p1.getChromosome() == p2.getChromosome()) {
      logFile << p1.getChromosome() << "\t" << p1.getLocation() << "\t" << p2.getLocation() 
              << "\t" << mSpecificityFrequencies[currentBound] << endl;
    }
    currentBound += 2;
  }
  logFile.close();

  
//   cout << endl;
//   for_each(shufflingeBoundaries.begin(), shufflingeBoundaries.end(), output__ );
//   cout << endl;
  
  //  Merge all boundaries to one list
  vector<size_t> allBoundaries;
  merge(distanceBoundaries.begin(), distanceBoundaries.end(), 
        probespecificityBoundaries.begin(), probespecificityBoundaries.end(), 
        back_inserter(allBoundaries));

  // Make sure no boundary comes twice -> delete duplicates
  allBoundaries.erase(unique(allBoundaries.begin(), allBoundaries.end()), allBoundaries.end());

  // The next routine needs the last probe index as boundary
  allBoundaries.push_back(probes.size() - 1);

  return createProbesetCompostionFromBoundaries(probes, allBoundaries, mOptimalProbesPerSet, mMaxProbesPerSet);  
}


/**
 * Constructor to set parameters for probeset shuffling method.
 *
 * @param maxDistance Maximum distance between succint probes within a set.
 * @param optimalProbesPerSet Number of probes optimally contained in a probeset.
 * @param maxProbesPerSet Number of probes maximal contained in a probeset.
 * @param threshold The specificity frequency value where to add a boundary.
 * @param hookcurveMovingAverageWindowSize Number of values to average the hookcurve over.
 * @param shufflingMovingAverageWindowSize Number of values to average specificity frequencies over.
 * @param detectKinkPointFunction Function to calculate the kink point.
 *
 */
ComposeProbesetsByCompleteShufflingMethod::ComposeProbesetsByCompleteShufflingMethod(int maxDistance, 
                                                                                     size_t optimalProbesPerSet, 
                                                                                     size_t maxProbesPerSet,
                                                                                     double threshold,
                                                                                     size_t hookcurveMovingAverageWindowSize,
                                                                                     size_t shufflingMovingAverageWindowSize,
                                                                                     DetectKinkPointPtr detectKinkPointFunction)
  :   mMaxDistance(maxDistance),
      mOptimalProbesPerSet(optimalProbesPerSet), mMaxProbesPerSet(maxProbesPerSet),
      mThreshold(threshold), 
      mHookcurveMovingAverageWindowSize(hookcurveMovingAverageWindowSize),
      mShufflingMovingAverageWindowSize(shufflingMovingAverageWindowSize),
      detectKinkPoint(detectKinkPointFunction)
{}

/**
 * Constructor 
 *
 * @param filename Name of the tab-file that contains genotype calls.
 */
ComposeProbesetsByGenotypeCall::ComposeProbesetsByGenotypeCall(const std::map<std::string,char> genotypeCalls,
                                                               const Haplotype haplotype)
  : mGenotypeCalls(genotypeCalls), mHaplotype(haplotype)
{}


// Read genotype calls from file
std::map<std::string,char>  
ComposeProbesetsByGenotypeCall::readGenotypeCalls(const std::string filename, const std::string chipname)
{
  // Define mapping from file info to CALL
  std::map<char,char> crlmmCallMap;
  crlmmCallMap['1'] = 'A';
  crlmmCallMap['2'] = 'H';
  crlmmCallMap['3'] = 'B';

  // Try to open file
  ifstream genotypeFile;
	genotypeFile.open(filename.c_str());
	if (!genotypeFile) {
    cerr << "Fatal error: genotype-call file " <<  filename << " does not exist." << endl;
    exit(0);
  }

  // Get column containing haplotypes of chip "chipname"
	string currentLine;
  getline(genotypeFile, currentLine);
  vector<string> lineTokens;
  lineTokens = stringutil::splitString(currentLine, "\t");
  vector<string>::iterator chipColumnIter = 
    find(lineTokens.begin(), lineTokens.end(), chipname);
  if (chipColumnIter == lineTokens.end()) {
    cerr << "Requested chip '" + chipname + "' does not occur in genotype call file " << endl;
    exit(0);
  }
  size_t chipColumnIndex = chipColumnIter - lineTokens.begin() + 1; // !!!! +1 for CRLMM files 

  // Load alleles into map structure
  map<string,char> genotypeCalls;
	while (!genotypeFile.eof()) {
    getline(genotypeFile, currentLine);
    lineTokens = stringutil::splitString(currentLine, "\t");
    //genotypeCalls[lineTokens[0]] = lineTokens[chipColumnIndex][0];
    genotypeCalls[lineTokens[0]] = crlmmCallMap[ lineTokens[chipColumnIndex][0] ];
  }
  genotypeFile.close();

  return(genotypeCalls);
}

std::map<std::string,char>  
ComposeProbesetsByGenotypeCall::readGenotypeCallsFromGDAS(const std::string filename, const std::string chipname)
{
  // Try to open file
  ifstream genotypeFile;
	genotypeFile.open(filename.c_str());
	if (!genotypeFile) {
    cerr << "Fatal error: genotype-call file " <<  filename << " does not exist." << endl;
    exit(0);
  }

  // Get column containing haplotypes of chip "chipname"
	string currentLine;
  getline(genotypeFile, currentLine);
  vector<string> lineTokens;
  lineTokens = stringutil::splitString(currentLine, "\t");
  vector<string>::iterator chipColumnIter = 
    find(lineTokens.begin(), lineTokens.end(), chipname);
  if (chipColumnIter == lineTokens.end()) {
    cerr << "Requested chip '" + chipname + "' does not occur in genotype call file " << endl;
    exit(0);
  }
  size_t chipColumnIndex = chipColumnIter - lineTokens.begin();

  // Load alleles into map structure
  map<string,char> genotypeCalls;
	while (!genotypeFile.eof()) {
    getline(genotypeFile, currentLine);
    lineTokens = stringutil::splitString(currentLine, "\t");
    genotypeCalls[lineTokens[0]] = lineTokens[chipColumnIndex][0];
    // @bug: last linetoken is processed twice, but that has no actual effect on the map
    // cout << lineTokens[0] << "\t" << lineTokens[chipColumnIndex][0] << endl;
  }
  genotypeFile.close();

  return(genotypeCalls);
}




/// Helper function that compares the snp allele (casting a PmMmProbePtr)
bool pmMmProbeAlleleEquals(const PmMmProbePtr& probe, char c) 
{
  // const SnpProbePtr snpProbe = boost::shared_static_cast<SnpProbe>(probe);
  // if (snpProbe) {
  //   return snpProbe->getAllele() == c; 
  // }
  assert(false);
  return false; // Avoid compiler warnings
}

/// Helper function that gets the snp allele (casting a PmMmProbePtr)
char getPmMmProbeAllele(const PmMmProbePtr& probe) 
{
  //const SnpProbePtr snpProbe = boost::shared_static_cast<SnpProbe>(probe);
  // assert(snpProbe);
  // return snpProbe->getAllele(); 
  assert(false);
}

// Returns the lexicographical sorted pair of alleles occuring in a probeset
std::pair<char, char> getSortedProbesetAlleles(const PmMmProbePtrVector::const_iterator beginProbe, 
                                               const PmMmProbePtrVector::const_iterator endProbe)
{
  char firstAllele = getPmMmProbeAllele(*beginProbe);
  PmMmProbePtrVector::const_iterator currentProbe = beginProbe + 1;
  while (getPmMmProbeAllele(*currentProbe) == firstAllele) { // There must be two alleles
    assert(currentProbe != endProbe);
    ++currentProbe;
  }
  char secondAllele = getPmMmProbeAllele(*currentProbe);
  return make_pair(min(firstAllele, secondAllele), 
                   max(firstAllele, secondAllele));  
}


/**
 * Creates probeset composition based on the probeset id property of the probes.
 *
 * @note Not well tested.
 *
 * @todo Rewrite such that only one haplotype-set is created
 *
 * @param probes List of probes
 * 
 * @return List of probesets.
 */ 
ProbesetVector ComposeProbesetsByGenotypeCall::operator() (PmMmProbePtrVector& probes)
{
  assert(probes.size() > 0);
  
  // Sort probes to make sure one set's probes are adjacent to each other
  sort(probes.begin(), probes.end(), ProbesetIdIsSmaller());

  // Define 3 lists of probesets
  ProbesetVector probesetsAllele, probesetsCross, probesetsHetero;

  PmMmProbePtrVector::const_iterator currentProbesetBegin = probes.begin();
  PmMmProbePtrVector::const_iterator currentProbesetEnd = probes.begin();

  // For all probesets
  while (currentProbesetBegin != probes.end()) {
    // Find end of current probeset
    string currentProbesetId = (*currentProbesetBegin)->getProbesetId();
    while (currentProbesetEnd != probes.end() &&
           (*currentProbesetEnd)->getProbesetId() == currentProbesetId) {
      ++currentProbesetEnd;
    }

    // Get probeset allele on analyzed chip
    char currentAllele = mGenotypeCalls[currentProbesetId];

    // Put probes into ...
    if (currentAllele == 'H') { // ... a single probeset if set is heterocygote
      probesetsHetero.push_back( Probeset(currentProbesetBegin, currentProbesetEnd, currentProbesetId));
    }
    else {  // .. or into different probesets depending on probe allele
      // Get base b that corresponds to present allele of current probeset
      char currentAlleleBase; 
      // @note We use the implicit knowledge that allele A is the lexicographic smallest base!
      if (currentAllele == 'A') {
        currentAlleleBase = getSortedProbesetAlleles(currentProbesetBegin, currentProbesetEnd).first;
      }
      else {
        currentAlleleBase = getSortedProbesetAlleles(currentProbesetBegin, currentProbesetEnd).second;
      }

      // Sort probes into bins of 'allele' (has that base b) and 'cross' (does not have base b)
      PmMmProbePtrVector probesCross, probesAllele;      
      remove_copy_if(currentProbesetBegin, currentProbesetEnd,  // remove_copy_if equals !copy_if
                     back_inserter(probesAllele), boost::bind(getPmMmProbeAllele, _1) != currentAlleleBase);
      probesetsAllele.push_back( Probeset(probesAllele));
      remove_copy_if(currentProbesetBegin, currentProbesetEnd,  // remove_copy_if equals !copy_if
                     back_inserter(probesCross), boost::bind(getPmMmProbeAllele, _1) == currentAlleleBase); 
      probesetsCross.push_back( Probeset(probesCross));
    }
    
    // Update iterator and id
    currentProbesetBegin = currentProbesetEnd;    
  }

  // Return desired probeset
  if (mHaplotype == kAllele) {
    return probesetsAllele; 
  }
  else if (mHaplotype == kCross) {
    return probesetsCross; 
  }
  else
    return probesetsHetero; 
}


/**
 * Creates probeset composition based on the probeset id and allele type of the probes.
 *
 * @note Not well tested.
 *
 * @param probes List of SNP probes
 * 
 * @return List of probesets.
 */ 
SnpProbesetVector ComposeProbesetsByGenotypeCall::composeSnpProbesets(SnpProbePtrVector& probes, std::map<std::string,char>& mGenotypeCalls)
{
  assert(probes.size() > 0);
  
  // Sort probes to make sure one set's probes are adjacent to each other
  sort(probes.begin(), probes.end(), ProbesetIdAndAlleleAndInterrogationIsSmaller());

  // Define 3 lists of probesets
  SnpProbesetVector probesets;

  SnpProbePtrVector::const_iterator currentProbesetBegin = probes.begin();
  SnpProbePtrVector::const_iterator currentProbesetEnd = probes.begin();

  // FOR ME: PRINT ORDERED SNP PROBES !!!!!!!!!!!!
  // transform(probes.begin(), probes.end(), ostream_iterator<string>(cout), boost::mem_fn(&SnpProbe::toString));


  // For all probesets
  while (currentProbesetBegin != probes.end()) {
    // Find end of current probeset
    string currentProbesetId = (*currentProbesetBegin)->getProbesetId();
    while (currentProbesetEnd != probes.end() &&
           (*currentProbesetEnd)->getProbesetId() == currentProbesetId) {
      ++currentProbesetEnd;
    }

   

    // Get probeset allele on analyzed chip (if available)
    // char currentCall = mGenotypeCalls[currentProbesetId];
    std::map<std::string,char>::const_iterator it = mGenotypeCalls.find(currentProbesetId);
    if ( it != mGenotypeCalls.end()) { // If call available     
      probesets.push_back( SnpProbeset(currentProbesetBegin, currentProbesetEnd, 
                                       currentProbesetId, it->second)); // currentCall));
    }
    else { // What if call not available?? And why is it not available??
      cout << currentProbesetId << " is missing in genotype call file. Removing." << endl;
    }


  //   // Put probes into ...
  //   if (currentAllele == 'H') { // ... a single probeset if set is heterocygote
  //     probesetsHetero.push_back( Probeset(currentProbesetBegin, currentProbesetEnd, currentProbesetId));
  //   }
  //   else {  // .. or into different probesets depending on probe allele
  //     // Get base b that corresponds to present allele of current probeset
  //     char currentAlleleBase; 
  //     // @note We use the implicit knowledge that allele A is the lexicographic smallest base!
  //     if (currentAllele == 'A') {
  //       currentAlleleBase = getSortedProbesetAlleles(currentProbesetBegin, currentProbesetEnd).first;
  //     }
  //     else {
  //       currentAlleleBase = getSortedProbesetAlleles(currentProbesetBegin, currentProbesetEnd).second;
  //     }

  //     // Sort probes into bins of 'allele' (has that base b) and 'cross' (does not have base b)
  //     SnpProbePtrVector probesCross, probesAllele;      
  //     remove_copy_if(currentProbesetBegin, currentProbesetEnd,  // remove_copy_if equals !copy_if
  //                    back_inserter(probesAllele), boost::bind(getSnpProbeAllele, _1) != currentAlleleBase);
  //     probesetsAllele.push_back( Probeset(probesAllele));
  //     remove_copy_if(currentProbesetBegin, currentProbesetEnd,  // remove_copy_if equals !copy_if
  //                    back_inserter(probesCross), boost::bind(getPmMmProbeAllele, _1) == currentAlleleBase); 
  //     probesetsCross.push_back( Probeset(probesCross));
  //   }
    
    // Update iterator and id
    currentProbesetBegin = currentProbesetEnd;    
  }
  return probesets;
}

/**
 * Debug function to print pm values and specificity frequencies
 * into a file
 *
 */
void exportIntensitiesAndFrequencies(const PmMmProbePtrVector& probes, const std::valarray<double>& specificityFrequencies)
{
  ofstream logFile;
  logFile.open("ClusteringMethods.txt");
  logFile << "Chromosome\tLocation\tlog PM\tDeltaLogI\tShuffle Frequency" << endl;
  for(size_t i = 0; i < probes.size(); ++i) {
    PmMmProbe p = *probes[i];
    logFile << p.getChromosome() << "\t" << p.getLocation() << "\t" << log10(p.getPm()) << "\t" << p.getDeltaLogI() << "\t" << specificityFrequencies[i] << endl;
  }
  logFile.close();
}

/**
 * Creates the probeset composition.
 *
 *
 * @note Internally creates another HookcurveAnalyzer to use the other shuffling functions.
 *  
 */ 
ProbesetVector ComposeProbesetsByCompleteShufflingMethod::operator() (const PmMmProbePtrVector& probes) 
{
  valarray<double> specificityFrequencies(0.0, probes.size());
  int validShufflingResults = 0;

  for (size_t shuffleIndex = 0; shuffleIndex < 5; ++shuffleIndex) {
    cout << "Beginning a shuffle" << endl;
    // -----------------------------------------------------------------------------------------
    ProbesetCompositionFunction composeProbeset = ComposeProbesetsByShuffling(shuffleIndex, mMaxDistance, 
                                                                              mOptimalProbesPerSet, mMaxProbesPerSet);
    Chip aChip = Chip(probes, composeProbeset(probes));
    HookModel simpleHookModel = HookModel();
    HookcurveAnalyzer aArray = HookcurveAnalyzer(aChip, simpleHookModel/*, (*mIntensityArraysPtr)*/);

    // Calculate kink point and document results to a file
    detectKinkPoint->setExportFilename(string("testplot_") + lexical_cast<std::string>(shuffleIndex) + string(".dat"));
    cout << "Calculating hookcurve and kink point." << endl;
    //IntensityArrayTable& iTable = (*mIntensityArraysPtr);
//     ProbesetIntensitiesArray sumLogIs   =  ( *(*mIntensityArraysPtr)[kPmI] + *(*mIntensityArraysPtr)[kMmI]  ) * 0.5;
//     ProbesetIntensitiesArray deltaLogIs =    *(*mIntensityArraysPtr)[kPmI] - *(*mIntensityArraysPtr)[kMmI];

    ProbesetIntensitiesArrayPtr logPmIs = Probeset::initializeIntensitiesArray(aChip.getProbesets(), kPmI);
    ProbesetIntensitiesArrayPtr logMmIs = Probeset::initializeIntensitiesArray(aChip.getProbesets(), kMmI);
    ProbesetIntensitiesArray sumLogIs   =  (*logPmIs + *logMmIs) * 0.5;
    ProbesetIntensitiesArray deltaLogIs =  (*logPmIs - *logMmIs);

    IntensityMappingPtr hookcurve = aArray.getAveragedHookcurvePlot(sumLogIs, deltaLogIs, mHookcurveMovingAverageWindowSize);
    IntensityPair intersectionCoordinates = (*detectKinkPoint)(MathUtil::digitizeCurveNew(*hookcurve));
//    IntensityPair intersectionCoordinates = (*detectKinkPoint)(*hookcurve);

    // Only use these shuffling values if we could calculate the intersection point      
    cout << "Calculated Point: " << intersectionCoordinates.first << endl;
    if (!isnan(intersectionCoordinates.first)) {
      specificityFrequencies += simpleHookModel.markPresentProbes(sumLogIs, intersectionCoordinates.first);
      ++validShufflingResults;
    }
  }
  
  cout << "Finished final run" << endl;
  specificityFrequencies /= (double)validShufflingResults;
  cout << "Beginning final shuffle" << endl;
  
  // DEBUG
  exportIntensitiesAndFrequencies(probes, specificityFrequencies);
//   exit(0);

  ProbesetCompositionFunction composeFinalProbeset = 
    ComposeProbesetsBySpecificityFrequencies(specificityFrequencies, 
                                             mThreshold, 
                                             mShufflingMovingAverageWindowSize, 
                                             mMaxDistance,
                                             mOptimalProbesPerSet, 
                                             mMaxProbesPerSet);
  return composeFinalProbeset(probes);
}

std::vector<size_t> getChromosomeBoundariesFromBedfile(const std::string filename, const std::string chromosome)
{
  // Define line eintries of bedfile
  enum LineentriesOfPmapfile {
    kLeChromsome = 0,
    kLeSegStart,
    kLeSegEnd
  }; 

  // Try to open file
  ifstream bedFile;
	bedFile.open(filename.c_str());
	if (!bedFile) {
    // raise file not found exception
    throw invalid_argument("Could not open bedfile.");
	}

  // Omit file header
	string currentLine;
  getline(bedFile, currentLine); // use std::getLine instead ifstram::getLine
                                  // since it gives back string

	// Read file line-by-line
  std::vector<size_t> boundaries;  
	while (!bedFile.eof()) {
    getline(bedFile, currentLine);
    vector<string> lineTokens = stringutil::splitString(currentLine, "\t");
    
    // cout << lineTokens[kLeChromsome] << endl;
    if (lineTokens.size() > kLeSegEnd && lineTokens[kLeChromsome] == chromosome) {
      boundaries.push_back(lexical_cast<size_t>(lineTokens[kLeSegStart]));
      boundaries.push_back(lexical_cast<size_t>(lineTokens[kLeSegEnd]));
    }
  }
  bedFile.close();
  
  // Sort and remove duplicates
  sort(boundaries.begin(), boundaries.end());
  boundaries.erase(unique(boundaries.begin(), boundaries.end()),
                   boundaries.end());

  return boundaries;
}


/**
 * DEBUG 
 * REads probeset boundaries from a given BED-file
 *
 * @pre Probes must be sorted with respect to chromosomes
 */
ProbesetVector ComposeProbesetsFromBedfileSegments::operator() (const PmMmProbePtrVector& probes) 
{
  // Test preconditions
  assert(probes.size() > 0);

  // Create boundaries for probes that are too distant
  vector<size_t> distanceBoundaries = defineDistanceConstrainedBoundaries(probes, mMaxDistance);

  // For each chromsome, insert boundaries from bedfile
  vector<size_t> bedfileBoundaries;
  string currentChromsome = "";
  vector<size_t> currentChromosomeLocationBoundaries;
  size_t segmentBoundaryIndex = 0;

  for (size_t probeIndex = 0; probeIndex < probes.size() - 1; ++probeIndex) {
    const PmMmProbe& p = *probes[probeIndex];

    // When finding a new chromsome, load new segment boundaries and reset segment index
    if (p.getChromosome() != currentChromsome) {
      currentChromsome = p.getChromosome();
      currentChromosomeLocationBoundaries = getChromosomeBoundariesFromBedfile(mBedfilename, currentChromsome);
      segmentBoundaryIndex = 0;      
    }

    // Skip probe when segment does not exist
    if (currentChromosomeLocationBoundaries.size() == 0) {
      continue;
    }

    // Increment segent index, while its genomic location <= p.location
    while (currentChromosomeLocationBoundaries[segmentBoundaryIndex] <= (size_t)p.getLocation()) {
      ++segmentBoundaryIndex;
    }

    // Insert boundary if next segement boundary < next probes location
    if (currentChromosomeLocationBoundaries[segmentBoundaryIndex] < (size_t)probes[probeIndex + 1]->getLocation()) {
//       cout << currentChromsome << "\t" << p.getLocation() << "\t" << probes[probeIndex + 1]->getLocation()
//            << "\t" << currentChromosomeLocationBoundaries[segmentBoundaryIndex] << endl;
      bedfileBoundaries.push_back(probeIndex);
    }
  }

  ofstream logFile;
  logFile.open("ShufflingBoundaries.wig");
  logFile << "track type=wiggle_0 name=\"Bed Format\" description=\"ShufflingBoundaries\" visibility=full" << endl;
  size_t currentBound = 0;
  while (currentBound < bedfileBoundaries.size() - 1) { // while at least two entries left
    const PmMmProbe& p1 = *probes[bedfileBoundaries[currentBound]];
    const PmMmProbe& p2 = *probes[bedfileBoundaries[currentBound + 1]];
    if (p1.getChromosome() == p2.getChromosome()) {
      logFile << p1.getChromosome() << "\t" << p1.getLocation() << "\t" << p2.getLocation() 
              << "\t1" << endl;
    }
    currentBound += 2;
  }
  logFile.close();


  //  Merge all boundaries to one list
  vector<size_t> allBoundaries;
  merge(distanceBoundaries.begin(), distanceBoundaries.end(), 
        bedfileBoundaries.begin(), bedfileBoundaries.end(), 
        back_inserter(allBoundaries));

  // Make sure no boundary comes twice -> delete duplicates
  allBoundaries.erase(unique(allBoundaries.begin(), allBoundaries.end()), allBoundaries.end());

  // The next routine needs the last probe index as boundary
  allBoundaries.push_back(probes.size() - 1);

  return createProbesetCompostionFromBoundaries(probes, allBoundaries, mOptimalProbesPerSet, mMaxProbesPerSet);  
}
