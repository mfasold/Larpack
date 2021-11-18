/**
 * @file SensitivityProfile.cpp Defines methods to calculate sequence specific intensity increments.
 * @author $Author$ 
 * @author Mario Fasold
 * @author Jan Bruecker
 * @date $Date$
 */
#include <fstream>
#include <valarray>
#include <boost/lexical_cast.hpp>
#include "boost/tuple/tuple.hpp"
#include <boost/progress.hpp>
using boost::lexical_cast;

#include "Matrix.hpp"
#include "SensitivityProfile.hpp"
#include "Probe.hpp"
#include "PmMmProbe.hpp"
#include "Probeset.hpp"
#include "SequenceUtil.hpp"
#include "MathUtil.hpp"
#include "ValarrayUtil.hpp"
#include "StringUtil.hpp"
        
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace larrpack;
using namespace sequenceutil;
// using namespace boost;
namespace su = sequenceutil;




/**
 * Calculates the sensitivity values \f$ \sigma_k \f$.
 *
 * This function is only an alternative interface to computeProfile
 * for backward compatibility,
 *
 * @param probesets List of probesets used to determine the profile.
 * @param getProbeIntensities Probe function used to get the intensities 
 *        of a probeset's probes. (either Probeset::GetPmLogIs or GetMmLogIs).
 * @param getProbeSequence Probe function used to get the sequence (either
 *        PmMmProbe::getSequence or PmMmProbe::getSequenceMm).
 *
 */
void SensitivityProfile::computeProfile(const ProbesetVector& probesets,
                                        const ProbesetIntensitiesArray& probesetIntensitiesArray,
                                        const PmMmProbeSequenceFunction& getProbeSequence,
                                        const size_t probesetLimit)
{
  // Create list of probe sequences
  ProbesetSequencesArray probsetSequencesArray;
  for (ProbesetVector::const_iterator pi = probesets.begin(); pi != probesets.end(); ++pi) {
    probsetSequencesArray.push_back(pi->getProbeSequences(getProbeSequence));
  } 

  computeProfile(probsetSequencesArray, probesetIntensitiesArray, probesetLimit);
}


/**
 * Constructor. Initializes with an empty profile (zero increment).
 *
 * @param sequenceLength Length of the sequence.
 * @param modelRank Size of the base tuples within the current model
 * (mononucleotides = 1, dinucleotides = 2, ...)
 * 
 */
SimpleSensitivityProfile::SimpleSensitivityProfile(size_t sequenceLength, 
                                       unsigned short modelRank) 
  : mSequenceLength(sequenceLength), 
    mModelRank(modelRank),
    mProfile(IntensityArray(0.0, getProfileSize(modelRank, sequenceLength, kNucleotideCount)))
{}

/**
 * Constructor for alternate profiles. Initializes with an empty profile (zero increment).
 *
 * @param sequenceLength Length of the sequence.
 * @param modelRank Size of the base tuples within the current model
 * (mononucleotides = 1, dinucleotides = 2, ...)
 * @param profileSize Size of the profile.
 * 
 * 
 */
SimpleSensitivityProfile::SimpleSensitivityProfile(size_t sequenceLength, 
                                       unsigned short modelRank,
                                       size_t profileSize) 
  : mSequenceLength(sequenceLength), 
    mModelRank(modelRank),
    mProfile(IntensityArray(0.0, profileSize))
{}


/**
 * Returns the size of an sequence affinity profile.
 * The size of the internal profile is computed as
 * \f$ 4^modelRank * (sequenceLength - modelRank + 1) \f$,
 * assumig we have 4 bases. The latter term arises from the
 * fact that we cannot compute the affnity for the last base 
 * positions in a sequence if we use model ranks larger than 
 * single-base.
 *
 * 
 */
size_t SimpleSensitivityProfile::getProfileSize(unsigned short modelRank, size_t sequenceLength, size_t nucleotideCount)
{
  return (size_t)pow((float)nucleotideCount,modelRank)*(sequenceLength - modelRank + 1);
}

SensitivityProfilePtr SimpleSensitivityProfile::cloneWithNewRank(const size_t modelRank) const
{
  SensitivityProfilePtr newProfile = SensitivityProfilePtr(new SimpleSensitivityProfile(*this, modelRank));
  return newProfile; // @note Unnamed shared_ptr temporaries should be avoided
}


/**
 * Copy-Constructor that takes another profile of arbitrary modelRank 
 * \f$ r_2 \f$ to build a profile of modelRank \f$ r_1 \leq r_2 \f$.
 *
 * For \f$ r_1 = r_2 + 1 \f$ the function computes 
 * \f[ \sigma_k(B) = 1/4 ( \sum_{\tilde{B} \in A,C,G,T} \sigma_{k-1} (B,\tilde{B}) 
 *     + \sum_{\tilde{B} \in A,C,G,T} \sigma_k(\tilde{B},B) )
 *  , \forall 2 \leq k \leq n - 1
 * \f]
 *
 * @bug Apparently, does not work correctly for differences in modelRank >= 2.
 *      While adding up the marginals, only baseIndex is used and therefore only works for
 *      rank 1!! 
 *
 * @note This function uses information about the sequenceToIndex format.
 */
SimpleSensitivityProfile::SimpleSensitivityProfile(const SimpleSensitivityProfile& s, unsigned short modelRank)
  : mSequenceLength(s.mSequenceLength), 
    mModelRank(modelRank),
    mProfile(IntensityArray(0.0, getProfileSize(modelRank, s.mSequenceLength, kNucleotideCount)))
{
  // If modelRank identical, just copy profile and we're done
  if (modelRank == s.mModelRank) {
    mProfile = s.mProfile;
    return;
  }

  // If the modelRank is bigger, run interpolation
  if (modelRank > s.mModelRank) {
    mProfile = interpolateFromProfile(s, modelRank);
    return;
  }  

  // Recursively get the rank of the source profile 
  // to modelRank + 1
  IntensityArray profile(0.0, getProfileSize(modelRank + 1, mSequenceLength, kNucleotideCount));
  size_t validSequenceIndices;
  if (modelRank + 1 == s.mModelRank) {
    profile = s.getProfile();
    validSequenceIndices = s.getValidSequenceCount();
  }
  else {
    SimpleSensitivityProfile sensitivityProfileIncrement = SimpleSensitivityProfile(s, modelRank + 1);
    profile = sensitivityProfileIncrement.getProfile();
    validSequenceIndices = sensitivityProfileIncrement.getValidSequenceCount();
  }

  // Add up marginal totals
  for(size_t baseIndex = 0; baseIndex < kNucleotideCount; ++baseIndex) {
    for(size_t k = 1; k < getValidSequenceCount() - 1; ++k) { // k = sequenceIndex
      mProfile[baseIndex*getValidSequenceCount() + k] = 
        IntensityArray(profile[slice(baseIndex*validSequenceIndices + k - 1, 
                                               kNucleotideCount, 
                                               kNucleotideCount*validSequenceIndices)]).sum() +
        IntensityArray(profile[slice(baseIndex*validSequenceIndices*kNucleotideCount + k, 
                                               kNucleotideCount, validSequenceIndices)]).sum();
      mProfile[baseIndex*getValidSequenceCount() + k] /= 4;
    }
      
    // Special cases for position 1 and n
    mProfile[baseIndex*getValidSequenceCount() + 0] = 0.25*
      IntensityArray(profile[slice(baseIndex*validSequenceIndices*kNucleotideCount + 0, 
                                             kNucleotideCount, validSequenceIndices)]).sum();
    mProfile[baseIndex*getValidSequenceCount() + getValidSequenceCount() - 1] = 0.25* 
      IntensityArray(profile[slice(baseIndex*validSequenceIndices + validSequenceIndices - 1, 
                                             kNucleotideCount,
                                             kNucleotideCount*validSequenceIndices)]).sum();
  }
}


/**
 * Computes the increment for sequence triple "sequence" at position pos
 *
 * @todo test function!!
 */
IntensityType SimpleSensitivityProfile::getSubTripleIncrement(const std::string& subTuple, size_t pos) const
{
  IntensityType sensitivity = 0.0;

  size_t subLength = 3;
  assert(subTuple.length() == subLength); 

  assert(mModelRank == 2); // currently only works for models of rank 2

  size_t completeLowerRankTuplesCount = subLength - mModelRank + 1; // how many tuples of lower rank suit in triple
  // Add increment for tuples that are fully in the subsequence    
  //for (size_t posIndex = pos; posIndex <= pos + subLength - mModelRank; ++posIndex) {
  for (size_t posIndex = 0; posIndex < completeLowerRankTuplesCount; ++posIndex) {
    // sensitivity += mProfile[sequenceutil::baseToIndex(sequence[posIndex])*mSequenceLength + posIndex];
    sensitivity += mProfile[getProfileIndex(subTuple.substr(posIndex,mModelRank), pos + posIndex)] / completeLowerRankTuplesCount;
  }


  // Add increment for borders with overlap 1
  // @optional: implement fast, but difficult to understand slice variant (that only works for r=2)
  if (pos > 0) { // if left margin available
    IntensityType marginalSensitivity = 0.0;
    for(size_t baseIndex = 0; baseIndex < kNucleotideCount; ++baseIndex) { // for all bases X in A,C,G,T
      marginalSensitivity += mProfile[getProfileIndex(su::indexToBase(baseIndex) + subTuple.substr(0, 1), 
                                                      pos - 1)];             // add XB
    }
    sensitivity += marginalSensitivity / 4 ; // normalize marginals by number of contributions
  }

  if (pos < getValidSequenceCount() - subLength + 1) { // if right margin available
    IntensityType marginalSensitivity = 0.0;
    for(size_t baseIndex = 0; baseIndex < kNucleotideCount; ++baseIndex) { // for all bases X in A,C,G,T
      marginalSensitivity += mProfile[getProfileIndex(subTuple.substr(subLength - 1, 1) + su::indexToBase(baseIndex), 
                                                      pos + subLength - 1)]; // add BX
    }
    sensitivity += marginalSensitivity / 4 ; // normalize marginals by number of contributions
  }

  // for modelRank3 add here more marginals!!

  return sensitivity;
}

/**
 * Interpolates a more sophisticated (bigger modelRank) profile
 * 
 * @param s Source profile
 * @param targetRank modelRank of target profile
 */
IntensityArray SimpleSensitivityProfile::interpolateFromProfile(const SimpleSensitivityProfile& s, unsigned short targetRank)
{
  unsigned short sourceRank = s.mModelRank;
  assert(targetRank > sourceRank);

  const size_t seqLength = s.getSequenceLength();
  IntensityArray interpolatedProfile(0.0, getProfileSize(targetRank, seqLength, kNucleotideCount));

  // Recursively run function until targetRank - 1 == s.mModelRank
  IntensityArray profile(0.0, getProfileSize(targetRank - 1, seqLength, kNucleotideCount));
  size_t validSequenceIndices;
  if (targetRank - 1 == s.mModelRank) {
    profile = s.getProfile();
    validSequenceIndices = s.getValidSequenceCount();
  }
  else {
    SimpleSensitivityProfile sensitivityProfileIncrement = SimpleSensitivityProfile(s, targetRank - 1);
    profile = sensitivityProfileIncrement.getProfile();
    validSequenceIndices = sensitivityProfileIncrement.getValidSequenceCount();
  }

  // Compute interpolation for rank difference 1
  for (size_t tupleIndex = 0; tupleIndex < su::getBasetupleCount(targetRank); ++tupleIndex) {
    for (size_t seqIndex = 0; seqIndex < su::getValidSequenceCount(seqLength, targetRank); 
         ++seqIndex) {
      string tuple = su::indexToSequence(tupleIndex, targetRank);
      IntensityType contribution = 0;
      contribution += profile[getProfileIndex(tuple.substr(0,sourceRank), seqIndex, seqLength, sourceRank)];
      contribution += profile[getProfileIndex(tuple.substr(1,sourceRank), seqIndex + 1, seqLength, sourceRank)];
      for (size_t overhangCount = 1; overhangCount < sourceRank; ++overhangCount) {
        for (size_t tupleFraction = 0; tupleFraction < su::getBasetupleCount(overhangCount); ++tupleFraction) {
          if (seqIndex >= sourceRank - overhangCount) { // @todo Check this // If i can determine tuples ending at
            contribution += profile[getProfileIndex(su::indexToSequence(tupleFraction, overhangCount) 
                                                    + tuple.substr(0,sourceRank - overhangCount), 
                                                    seqIndex - overhangCount, seqLength, sourceRank)] / 4;            
          }    
//         if (seqIndex > 0) {
//           contribution += profile[getProfileIndex(indexToSequence(tupleFraction, 1) + tuple.substr(0,1), 
//                                                   seqIndex - 1, seqLength, sourceRank)] / 4;
//         }
          if (seqIndex < su::getValidSequenceCount(seqLength, targetRank) - (sourceRank - overhangCount)) {
            contribution += profile[getProfileIndex(tuple.substr(targetRank - overhangCount, overhangCount) 
                                                    + su::indexToSequence(tupleFraction, overhangCount), 
                                                    seqIndex + targetRank - overhangCount, seqLength, sourceRank)] / 4;          
          }    
        }
          
//         if (seqIndex < su::getValidSequenceCount(seqLength, targetRank) - 1) {
//           contribution += profile[getProfileIndex(tuple.substr(2,1) + indexToSequence(tupleFraction, 1), 
//                                                   seqIndex + 2, seqLength, sourceRank)] / 4;            
//         }
      }

//       cout << "tupleIndex: " << tupleIndex << endl;
//       cout << "seqIndex: " << seqIndex << endl;
//       cout << "getProfileIndex(tupleIndex, seqIndex, seqIndex, targetRank): " << getProfileIndex(tupleIndex, seqIndex, seqLength, targetRank) << endl;
      interpolatedProfile[getProfileIndex(tupleIndex, seqIndex, seqLength, targetRank)] = contribution / 2;
    }        
  }
  return interpolatedProfile;
}
    
/**
 * Returns the sequence specific increment 
 * \f$ \delta S = \sum_{k=1}^{sequenceLength} \sigma_k( sequence_k ) \f$
 *
 * @param sequence The nucleotide sequence.
 */
IntensityType SimpleSensitivityProfile::getSequenceIncrement(const std::string& sequence) const
{
  assert(sequence.size() == mSequenceLength);

  IntensityType sensitivity = 0.0;
  for (size_t posIndex = 0; posIndex < getValidSequenceCount(); ++posIndex) {
    // sensitivity += mProfile[sequenceutil::baseToIndex(sequence[posIndex])*mSequenceLength + posIndex];
    sensitivity += mProfile[sequenceutil::sequenceToIndex(sequence.substr(posIndex,mModelRank))*getValidSequenceCount()
                            + posIndex];
  }
  return sensitivity;
}

/**
 * Returns the sigma for a base/oligonucleotide index at a given position.
 * 
 * @param baseIndex Numerical value representing the base/oligonucleotide
 * @param position Position in the probe
 *
 * @return Sigma
 */
IntensityType SimpleSensitivityProfile::getSigma(size_t baseIndex, size_t position) const
{
  return mProfile[baseIndex*getValidSequenceCount() + position];
}

/**
 * Sets the effect of the middle base to zero. Currently, it
 * sets the effect off *all* base tuples, which include the middle base,
 * to zero. This is a oversimplification for models with rank > 1.
 *
 */
void SimpleSensitivityProfile::zeroMiddlebase()
{    
  size_t middlebaseIndex = sequenceutil::getMiddlebaseIndex(mSequenceLength);
  // For all k in sigma_k(base tuple) which include the middle base
  for(size_t posIndex = middlebaseIndex - mModelRank + 1; posIndex <= middlebaseIndex; ++posIndex) {
    // For all base tupes in sigma_k(base tuple)
    for(size_t baseIndex = 0; baseIndex < getBasetupleCount(); ++baseIndex) {
      mProfile[baseIndex*getValidSequenceCount() + posIndex] = 0;   
    }
  }
}

/**
 * Returns the rank of the sequence affinity model, e.g. single-base 
 * model (rank 1), dinucleotide model (rank 2) and so forth.
 * 
 * @return Rank of the sequence affinity model.
 */
unsigned short SimpleSensitivityProfile::getModelRank() const
{
  return mModelRank;
}

/**
 * Exports the profile to a file in tabular text format.
 *
 * @param filename Filename.
 */

void SimpleSensitivityProfile::exportToDatafile(const std::string filename) const
{
  ofstream sensitivityProfileFile;
  sensitivityProfileFile.open(filename.c_str());

  // For all positions
  for(size_t pos = 0; pos < mProfile.size(); ++pos) {
    size_t seqPos = (pos % getValidSequenceCount());
    size_t baseIndex = size_t (pos / getValidSequenceCount());
    string baseString = indexToSequence(baseIndex);
    sensitivityProfileFile << seqPos + 1 << "\t" << baseString << "\t" << mProfile[pos] << endl;
  }
//   printValarrayLinebyline(sensitivityProfileFile, mProfile);
  sensitivityProfileFile.close();

  // Create Gnuplot file
  if (true) {
    ofstream gnuplotFile;
    gnuplotFile.open((stringutil::splitString(filename, ".")[0] + ".gnuplot").c_str());

    gnuplotFile << "set term png large size 800,600" << endl;
    gnuplotFile << "set output \"" << stringutil::splitString(filename, ".")[0]+".png\"" << endl;
    gnuplotFile << "set xlabel \"Sequence Index\"" << endl;
    gnuplotFile << "set ylabel \"Sensitivity\"" << endl;
    gnuplotFile << "set grid" << endl;
    gnuplotFile << "show grid" << endl;
    gnuplotFile << "set xrange[0:" << getValidSequenceCount() << "]" << endl;
//     gnuplotFile << "set title \"" << filename << "\"" << endl;
//     gnuplotFile << "set key below" << endl;

    size_t profileIndex = 0;
    gnuplotFile << "plot ";
    size_t baseCount = mProfile.size() / getValidSequenceCount();
    for (size_t i = 0; i < baseCount; ++i) {
      if (i > 0) {
        gnuplotFile  << ", ";
      }
      gnuplotFile  << "\"" << filename << "\" u 3 every ::" << profileIndex << "::" 
                   << (profileIndex + getValidSequenceCount() - 1) << " w linesp lw 2 title '" 
                   << indexToSequence(i) << "'";
      profileIndex += getValidSequenceCount();
    }
    // gnuplotFile << "0 lw 1 lt -1 title ''" << endl;
    gnuplotFile.close();
  }
}


/**
 * Imports the profile from a tabular text format file.
 *
 * @pre The modelRank of the profile to be read must equal mModelRank
 *
 * @todo Support reading of profiles of arbitrary size.
 *
 * @param filename Filename.
 */
void SimpleSensitivityProfile::importFromDatafile(const std::string filename) 
{
  // define line eintries of genomic file
  enum LineentriesOfProfilefile {
    kLeSequenceIndex = 0,
    kLeBaseTuple,
    kLeSequenceContribution
  }; 

  ifstream sensitivityProfileFile;
  sensitivityProfileFile.open(filename.c_str());

	if (!sensitivityProfileFile) {
    // raise file not found exception
    throw invalid_argument("Could not open profile file " + filename);
	}

  // For all lines in file
  vector<string> lineTokens;
	string currentLine;
	while (!sensitivityProfileFile.eof()) {
    // Get and type-cast tokens
    getline(sensitivityProfileFile, currentLine);

    if (currentLine.empty() || currentLine[0] == '#' ) {
      continue;
    }
    lineTokens = stringutil::splitString(currentLine, "\t");
    size_t seqIndex = lexical_cast<size_t>(lineTokens[kLeSequenceIndex]) - 1;
    IntensityType contribution = lexical_cast<IntensityType>(lineTokens[kLeSequenceContribution]);

    // Make sure preconditions are fulfilled
    assert(lineTokens[kLeBaseTuple].size() == mModelRank);

    // Enter data to existing profile
    mProfile[sequenceutil::sequenceToIndex(lineTokens[kLeBaseTuple]) * getValidSequenceCount() + seqIndex] 
      = contribution; 
  }

  sensitivityProfileFile.close();
}



/**
 *  Helper function: internal loop for the sequence profile computation. Sourced out to allow
 *  parallelization.
 *
 */
void SimpleSensitivityProfile::probeIterationHelper(math::Matrix<double>& delta, valarray<double>& theta, 
                                                    const ProbesetSequencesArray& probsetSequences, 
                                                    const size_t& probesetIndex, 
                                                    const ProbesetIntensitiesArray& probesetIntensitiesArray,
                                                    const size_t& modelNucleotides, 
                                                    const size_t& validSequencePositions, 
                                                    size_t profileSize)
{

  math::Matrix<double> probesetPSWM = 
    sequenceutil::getPositionSpecificWeightMatrix(probsetSequences[probesetIndex], mModelRank);
    
  // This is not nicely done... we need to filter out all probes smaller 2 (they might produce problems in the log scale)
  size_t count = 0; // Count saves how many probes of a probeset are > 2
  for (size_t probeIndex = 0; probeIndex < probsetSequences[probesetIndex].size(); ++probeIndex)
    if (MathUtil::gExp10(probesetIntensitiesArray[probesetIndex][probeIndex]) > 2)
      count++;
  IntensityArray intensitiesTmp(0.0, count); // IntensitiesTmp holds all intensities with unlogged values > 2
  size_t idx = 0;
  for (size_t probeIndex = 0; probeIndex < probsetSequences[probesetIndex].size(); ++probeIndex){
    if (MathUtil::gExp10(probesetIntensitiesArray[probesetIndex][probeIndex]) < 2)
      continue;
    else {
      intensitiesTmp[idx] = probesetIntensitiesArray[probesetIndex][probeIndex];
      idx++;
    }
  }
  if (intensitiesTmp.size() == 0){ // If all probes of the probeset have intensities < 2 quit for this probeset,
    return; // Leave the probeset loop
  }
    
  IntensityType average = computeAverageIntensity(intensitiesTmp);
  // Get intensity values
  IntensityArray intensities = probesetIntensitiesArray[probesetIndex]; // getProbeIntensities(currentProbeset);
    
  // Store all kroeneckers for one position to an array to speed things up    
  valarray<unsigned short> kroeneckers((unsigned short)0, modelNucleotides); // Init all with 0
    
  // Calculate the experimental sensitivities for the whole probeset array
  IntensityArray yExp = intensities - average;
    
    
    
  //     if (considerFrequencies) { // Here we use the frequencies
  //cout << "True" << endl;
  // Holds the relative frequencies of each base at each sequence position
  // P(B_i, Rho_i)_set (as in Preibisch, P.18)
  // math::Matrix<double> probesetPSWM = currentProbeset->getPositionSpecificWeightMatrix(mModelRank);
  // Average log(PM) intensity of the probeset
    
  // Creates matrix delta and theta sorted by bases first
  for (size_t probeIndex = 0; probeIndex < probsetSequences[probesetIndex].size(); ++probeIndex) {
    if (MathUtil::gExp10(probesetIntensitiesArray[probesetIndex][probeIndex]) < 2){
      continue;
    }  
    string sequence = probsetSequences[probesetIndex][probeIndex];
      
    // First we fill up this valarray, this will save runtime in the end.
    IntensityArray kroeMinusP(0.0, validSequencePositions * modelNucleotides);
    for (size_t pos = 0; pos < validSequencePositions; ++pos) {
      // Set kroenecker of the only one occuring base tuple to one
      size_t tupleIndexBase = sequenceToIndex(sequence.substr(pos, mModelRank));
      kroeneckers[tupleIndexBase] = 1;

      for (size_t base = 0; base < modelNucleotides; ++base) {
        // This way we have precalculated all the (kroenecker - P)s
        kroeMinusP[base * validSequencePositions + pos] = kroeneckers[base] - probesetPSWM[base][pos];
      }
      kroeneckers[tupleIndexBase] = 0; // Set it back to zero.
    }
    // After this precalculation we now can fill up the matrix:
      
    for (size_t pos1 = 0; pos1 < validSequencePositions; ++pos1) {
      for (size_t base1 = 0; base1 < modelNucleotides; ++base1) {
        // Fill up theta
        theta[base1 * validSequencePositions + pos1] += 
          yExp[probeIndex] * kroeMinusP[base1 * validSequencePositions + pos1];
        for (size_t pos2 = 0; pos2 < validSequencePositions; ++pos2) { 
          // Now only compute values for the lower triangle, that is
          // either base2 < base1 or (base2 == base1 and pos2 <= pos1)
          for (size_t base2 = 0; base2 < base1; ++base2) {
            delta[base1 * validSequencePositions + pos1][base2 * validSequencePositions + pos2] += 
              kroeMinusP[base1 * validSequencePositions + pos1] 
              * kroeMinusP[base2 * validSequencePositions + pos2];
          } // base2

            // We pull this case base1 == base2 out of the inner loop to save runtime
          size_t base2 = base1;
          if (pos2 <= pos1) {
            delta[base1 * validSequencePositions + pos1][base2 * validSequencePositions + pos2] += 
              kroeMinusP[base1 * validSequencePositions + pos1] 
              * kroeMinusP[base2 * validSequencePositions + pos2];
          }
        } // pos2 
      } // base1
    } // pos1
  } // for probes of currentProbeset
  
  //   ofstream thetaTmpBefore;
  //   thetaTmpBefore.open("thetaInnerCall.txt");
  //   thetaTmpBefore << theta << endl;
  //   thetaTmpBefore.close();
  //   return make_tuple(delta,theta);
}


/**
 * Calculates the sensitivity values \f$ \sigma_k \f$ as in Preibisch, Page 19. 
 * The function first computes matrix (12) and then solves it using SVD (with
 * zeroing).
 *
 * @param probesetSequencesArray List of probe sequences for all probesets, 
 *        used to determine the profile.
 * @param probesetIntensitiesArray List of probe intensities for all probesets
 * @param probesetLimit Number of probesets to be used for profile computation.
 *        They are sampled ot of all probesets.
 * 
 * @note The resulting profile is internally indexed in the form
 * [baseIndex*valdidSequencePositions + sequencePosition]
 *
 * @note The function is very time consuming and therefore 
 * optimized.
 *
 * @return Sensitivity profile (\f$ \sigma_k \f$). 
 */
void SimpleSensitivityProfile::computeProfile(const ProbesetSequencesArray& probesetSequencesArray,
                                              const ProbesetIntensitiesArray& probesetIntensitiesArray,
                                              const size_t probesetLimit)
{
  // delta * sigma = theta
  // Initialize matrix delta and valarray theta with zeros.
  // the matrix and the valarray will be used to calculate the valarray sigma by singular value decomposition.  
  const size_t profileSize = getProfileSize(mModelRank, kProbeSequenceLength, kNucleotideCount);
  math::Matrix<double> delta(profileSize, profileSize, 0.0);
  valarray<double> theta(0.0, profileSize);

  // Define some counters for fast access
  const size_t modelNucleotides = getBasetupleCount();
  const size_t validSequencePositions = getValidSequenceCount();

  cout << "modelNucleotides: " << modelNucleotides <<  flush;
  cout << ", profileSize: " << profileSize << flush;
  cout << ", number of probesets: " << probesetSequencesArray.size() << endl;


  // Iterate over probesets:
  size_t counter2 = 0; // This counter is to limit the number probesets for the profile calculation.
  size_t ithProbeset;
  if (probesetLimit > 0) {
    ithProbeset = probesetSequencesArray.size() / probesetLimit;
    if (ithProbeset == 0) {
      ithProbeset = 1;
    }
  }
  else { ithProbeset = 1; }

  
  ProbesetSequencesArray filteredProbesetSequencesArray; // We use only every ith probeset.
  ProbesetIntensitiesArray filteredProbesetIntensitiesArray;
  for (size_t psIdx = 0; psIdx < probesetSequencesArray.size(); ++psIdx) {
    if (! (counter2 % ithProbeset == 0)) {
      counter2++;
    }
    else {
      filteredProbesetSequencesArray.push_back(probesetSequencesArray[psIdx]);
      filteredProbesetIntensitiesArray.push_back(probesetIntensitiesArray[psIdx]);
      counter2++;
    }
  }
  //   ithProbeset = 1;
  if (ithProbeset > 1) {
    cout << "We use every " << ithProbeset << "th for sequenceprofile calculation."  
         << " ~ " <<  filteredProbesetSequencesArray.size() << " probesets." << flush;
  }
  
  //cout << endl;
  size_t probesetCount = filteredProbesetSequencesArray.size();
  size_t probesetIndex;
  
  // Calculation of sensitivity profiles is parallized if openmp is available
#ifdef _OPENMP
  cout << "Parallel Process. Running " << omp_get_max_threads() << " threads" << endl;
  // Since reduction doesn't work on matrices and valarrays we sort of build our own reduction by
  // assigning vectors in the size of the [max thread number] and adding the single entries when we leave the
  // parallel region.
  //  math::Matrix<double> m(profileSize, profileSize, 0.0);
  //valarray<double> a(0.0, profileSize);
  vector< math::Matrix<double> > tmpDeltas;
  vector< valarray<double> > tmpThetas;
  //  vector< math::Matrix<double> > tmpDeltas(omp_get_max_threads(), m);
  //  vector< valarray<double> > tmpThetas(omp_get_max_threads(), a);
  for (size_t i = 0; i < omp_get_max_threads(); i++) {
    //	math::Matrix<double> m(profileSize, profileSize, 0.0);
    //	tmpDeltas.push_back(m);  // terminate called after throwing an instance of 'std::bad_alloc'
    tmpDeltas.push_back(math::Matrix<double> (profileSize, profileSize, 0.0));
    //    valarray<double> a(0.0, profileSize);
    //    tmpThetas.push_back(a);
    tmpThetas.push_back(valarray<double> (0.0, profileSize));
  }
#endif
  boost::progress_display show_progress(probesetCount);
  
#ifdef _OPENMP
#pragma omp parallel for shared(tmpDeltas, tmpThetas) /*shared(theta, delta) firstprivate(deltaTmp, thetaTmp)*/ private(probesetIndex)
#endif
  for (probesetIndex = 0; probesetIndex < probesetCount; ++probesetIndex) {
#ifdef _OPENMP
    // If parallel each thread writes to its own matrix/vector to avoid race conditions.
    probeIterationHelper(tmpDeltas[omp_get_thread_num()], tmpThetas[omp_get_thread_num()], filteredProbesetSequencesArray, probesetIndex, filteredProbesetIntensitiesArray, modelNucleotides, validSequencePositions, profileSize);
#else
    // If serial, calculations can be done directly on delta and theta
    probeIterationHelper(delta, theta, filteredProbesetSequencesArray, probesetIndex, filteredProbesetIntensitiesArray, modelNucleotides, validSequencePositions, profileSize);
#endif
#ifdef _OPENMP
#pragma omp critical
#endif
    {
      ++show_progress;
    }
    
  }//for probesets
  //   cout << endl;
  
  
#ifdef _OPENMP
  // Add the thread specific deltas to the global data
  for (size_t i = 0; i < tmpDeltas.size(); i++) {
    assert(tmpDeltas.size() == tmpThetas.size());
    delta = delta + tmpDeltas[i];
    theta = theta + tmpThetas[i];
  }
  tmpDeltas.clear();
  tmpThetas.clear();
#endif
  
  // Calculate upper triangle of the matrix
  for (size_t row = 0; row < profileSize - 1; ++row) {
    for (size_t column = row; column < profileSize; ++column) {
      delta[row][column] = delta[column][row];
    }
  }

  // @debug: Dump matrix to file
  //   ofstream logFile;
  //   logFile.open("affinityDelta.txt");
  //   logFile << delta << endl;
  //   logFile.close();
  //   logFile.open("affinityTheta.txt");
  //   logFile <<  theta << endl;
  //   logFile.close();

  // Solve equation
  mProfile = MathUtil::solveEquationsSvd(delta, theta/*, 1.0*/);
}


/**
 * Checks if any probe within a probeset contains
 * one of a list of subsequences
 *
 * @todo Recode this in stl philosophy
 */
bool containsProbeWithSubsequences(const Probeset& ps, const std::vector<std::string> subsequences,
                                   const PmMmProbeSequenceFunction& getProbeSequence)
{
  // Loop through probes and subsequences
  for (PmMmProbePtrVector::const_iterator currentProbe = ps.beginProbe(); 
       currentProbe < ps.endProbe(); ++currentProbe) {
    for (std::vector<std::string>::const_iterator s = subsequences.begin(); s < subsequences.end(); ++s) {
      // Return true if one of the subsequences found in on of the probes
      if (getProbeSequence(*(*currentProbe)).find(*s) != string::npos) {
        return true;
      }
    }
  }
  return false;
}

// /**
//  * Testing function.
//  *
//  * @todo Remove.
//  * @todo Create Probeset*Ptr*Vector to avoid copying of vector elements
//  */
// void SimpleSensitivityProfile::
// computePartlyInterpolatedProfile(const ProbesetVector& probesets, 
//                                  const ProbesetIntensitiesArray& probesetIntensitiesArray,
//                                  const PmMmProbeSequenceFunction& getProbeSequence,
//                                  const std::vector<std::string> subsequences,
//                                  const size_t sourceRank,
//                                  const size_t probesetLimit,
//                                  const size_t probesetSelectLimit)
// {
//   // Devide probesets into two disjunct sets
//   // @note i don't use algorithms, as probesets and intensities must be filtered in parallel
//   ProbesetVector probesetsSelect, probesetsB;
//   ProbesetIntensitiesArray intensitiesSelect, intensitiesB;
//   for (size_t probesetIndex = 0; probesetIndex < probesets.size(); ++probesetIndex) {
//     const Probeset& currentProbeset = probesets[probesetIndex];
//     if (containsProbeWithSubsequences(currentProbeset, subsequences, getProbeSequence) 
//         && probesetsSelect.size() < probesetSelectLimit) {
//       probesetsSelect.push_back(currentProbeset);
//       intensitiesSelect.push_back(probesetIntensitiesArray[probesetIndex]);
//     }
//     else {
//       probesetsB.push_back(currentProbeset);
//       intensitiesB.push_back(probesetIntensitiesArray[probesetIndex]);
//     }
//   }

//   // Compute small profile for unselected probesets
//   SimpleSensitivityProfile profileSmall(mSequenceLength, sourceRank);
//   profileSmall.computeProfile(probesetsB, intensitiesB, getProbeSequence, probesetLimit);
//   // @debug write to file
//   profileSmall.exportToDatafile("profileSmall.dat");
  
//   // Interpolate higher-rank profile for unselected probesets
//   SimpleSensitivityProfile profileInterpolated(profileSmall, mModelRank);
//   profileInterpolated.exportToDatafile("profileInterpolated.dat");

//   // Compute higher-rank profile for selected probesets
//   SimpleSensitivityProfile profileSelect(mSequenceLength, mModelRank);
//   profileSelect.computeProfile(probesetsSelect, intensitiesSelect, getProbeSequence, probesetLimit);
//   profileSelect.exportToDatafile("profileSelect.dat");

//   // Create a mixed profile weighted with number of probesets in each set
//   IntensityType proportionSelect = (IntensityType) probesetsSelect.size() / probesets.size();
//   mProfile = proportionSelect*profileSelect.getProfile() 
//     + (1-proportionSelect)*profileInterpolated.getProfile();
  
//   // @debug print numbers
//   cout << "probesets.size(): " << probesets.size() << endl;
//   cout << "probesetsSelect.size(): " << probesetsSelect.size() << endl;
//   cout << "proportionSelect: " << proportionSelect << endl;
//   exportToDatafile("profileMixes.dat");
// }

/**
 * Computes the sum of squares error of the optimization problem (12) using 
 * the given profile \f$ \sigma_k \f$. 
 *
 * @param probesets List of probesets used to determine the profile.
 * @param getProbeIntensities Probe function used to get the intensities 
 *        of a probeset's probes. (either Probeset::GetPmLogIs or GetMmLogIs).
 * @param getProbeSequence Probe function used to get the sequence (either
 *        PmMmProbe::getSequence or PmMmProbe::getSequenceMm).
 *
 * @return Sum of squares error.
 *
 */
IntensityType 
SimpleSensitivityProfile::computeSumOfSquares(const ProbesetVector& probesets,
                                        const ProbesetIntensitiesArray& probesetIntensitiesArray,
                                        const PmMmProbeSequenceFunction& getProbeSequence)
{
  IntensityType sumOfSquares = 0;

  for (size_t probesetIndex = 0; probesetIndex < probesets.size(); ++probesetIndex) {
    const Probeset& currentProbeset = probesets[probesetIndex];


    // Get intensity values
    IntensityArray intensities = probesetIntensitiesArray[probesetIndex]; // getProbeIntensities(currentProbeset);
    IntensityType average = computeAverageIntensity(intensities);
        
    // Calculate the experimental sensitivities for the whole probeset array
    IntensityArray yExp = intensities - average;

    // 
    for (size_t probeIndex = 0; probeIndex < currentProbeset.getSize(); ++probeIndex) {
      string sequence = getProbeSequence(currentProbeset.getProbe(probeIndex));
      
      //IntensityType error = getSequenceIncrement(sequence) - yExp[probeIndex];
      IntensityType error = yExp[probeIndex] - getSequenceIncrement(sequence);
      sumOfSquares += error * error;
    }
  }
  return sumOfSquares;
}

/**
 * Computes the sum of squares error of the optimization problem (12) for sequences
 * containing a particular subsequece using  the given profile \f$ \sigma_k \f$. 
 *
 * @param probesets List of probesets used to determine the profile.
 * @param getProbeIntensities Probe function used to get the intensities 
 *        of a probeset's probes. (either Probeset::GetPmLogIs or GetMmLogIs).
 * @param getProbeSequence Probe function used to get the sequence (either
 *        PmMmProbe::getSequence or PmMmProbe::getSequenceMm).
 *
 */
void
SimpleSensitivityProfile::exportConditionalSumOfSquares(const ProbesetVector& probesets,
                                                  const ProbesetIntensitiesArray& probesetIntensitiesArray,
                                                  const PmMmProbeSequenceFunction& getProbeSequence,
                                                  const size_t baseCount,
                                                  const std::string filename,
                                                  const SensitivityProfilePtr profile)
{
//  IntensityArray sseProfile(0.0, getProfileSize(baseCount, mSequenceLength, kNucleotideCount));
//  std::valarray<size_t> frequencyProfile((size_t) 0, getProfileSize(baseCount, mSequenceLength, kNucleotideCount));
	
	 // Statistic 1: SumOfScqares for specific seqence positions
	size_t profileSize = getProfileSize(baseCount, kProbeSequenceLength, kNucleotideCount);
	IntensityArray sseProfile(0.0, profileSize);
	std::valarray<size_t> frequencyProfile((size_t) 0, profileSize);
	
	// Statistic 2: SumOfSquares counting if tuple occurs anywhere within sequence
	IntensityArray sseProfileSum(0.0, sequenceutil::getBasetupleCount(baseCount));
	std::valarray<size_t> frequencyProfileSum((size_t) 0, sequenceutil::getBasetupleCount(baseCount));
	
	// Statistic 3: Total error
	IntensityType totalError = 0;
	size_t probeCount = 0;
	
	// Statistic 4: Raw residuals
	IntensityArray sseRawProfile(0.0, getProfileSize(baseCount, kProbeSequenceLength, kNucleotideCount));

  // Statistic 6!: 
  vector< vector<IntensityType> > intensityVarianceProfile(getProfileSize(baseCount, kProbeSequenceLength, kNucleotideCount));


	                                              
	ofstream allResidualsFile; // @debug
	allResidualsFile.open("allResiduals.dat");  // @debug
	
	
	for (size_t probesetIndex = 0; probesetIndex < probesets.size(); ++probesetIndex) {
		const Probeset& currentProbeset = probesets[probesetIndex];

    // Get intensity values
    IntensityArray intensities = probesetIntensitiesArray[probesetIndex]; // getProbeIntensities(currentProbeset);
    IntensityType average = computeAverageIntensity(intensities);
        
    // Calculate the experimental sensitivities for the whole probeset array
    IntensityArray yExp = intensities - average;

    // 
    for (size_t probeIndex = 0; probeIndex < currentProbeset.getSize(); ++probeIndex) {
      // Compute probe's sum of squares error 
      string sequence = getProbeSequence(currentProbeset.getProbe(probeIndex));  
      //IntensityType error = getSequenceIncrement(sequence) - yExp[probeIndex];
      IntensityType error = yExp[probeIndex] - profile->getSequenceIncrement(sequence);
//       IntensityType error = MathUtil::gExp10(average +getSequenceIncrement(sequence)) 
//         - MathUtil::gExp10(intensities[probeIndex]);

      //error = error*error; // squared error term
      allResidualsFile << sequence << "\t" << error*error << "\t" << error  << "\t" 
                       << intensities[probeIndex]  << "\t" 
                       << yExp[probeIndex]  << "\t" 
//                        << average - profile->getSequenceIncrement(sequence)  << "\t" 
//                        << average + profile->getSequenceIncrement(sequence)  << "\t" 
//                        << exp10(average - profile->getSequenceIncrement(sequence))  << "\t" 
                       << endl;

      // Add error if respective subsequence occurs
 //     size_t maxSeqIndex = sequenceutil::getValidSequenceCount(mSequenceLength, baseCount);
      size_t maxSeqIndex = sequenceutil::getValidSequenceCount(kProbeSequenceLength, baseCount);
      for (size_t seqIndex = 0; seqIndex < maxSeqIndex; ++seqIndex) {
        size_t profileIndex = sequenceutil::sequenceToIndex(sequence.substr(seqIndex, baseCount))*maxSeqIndex 
          + seqIndex;
        //sseProfile[profileIndex] += error;
        sseProfile[profileIndex] += error*error; // squared error term
        sseRawProfile[profileIndex] += error;
        frequencyProfile[profileIndex] += 1;
        intensityVarianceProfile[profileIndex].push_back(yExp[probeIndex]);
      }
      
       // Add error if respective subsequence occurs at at all (Statistic 2)
      for (size_t baseIndex = 0; baseIndex < sequenceutil::getBasetupleCount(baseCount); ++baseIndex) {
        if (sequence.find(sequenceutil::indexToSequence(baseIndex, baseCount)) != string::npos) {
          sseProfileSum[baseIndex] += error*error;
          frequencyProfileSum[baseIndex] += 1;            
        }
      }
      
      // Increment total error
      totalError += error*error;
      probeCount++;
      
    }
  }
  
	allResidualsFile.close();  // @debug
	
  // @debug
  // Save the sum of squares "profile" to file
  ofstream sensitivityProfileFile;
  sensitivityProfileFile.open(filename.c_str());
  
  // For all positions
  for(size_t pos = 0; pos < sseProfile.size(); ++pos) {
//    size_t seqPos = (pos % sequenceutil::getValidSequenceCount(mSequenceLength, baseCount));
//    size_t baseIndex = size_t (pos / sequenceutil::getValidSequenceCount(mSequenceLength, baseCount));
	size_t seqPos = (pos % sequenceutil::getValidSequenceCount(kProbeSequenceLength, baseCount));
	size_t baseIndex = size_t (pos / sequenceutil::getValidSequenceCount(kProbeSequenceLength, baseCount));
	  
    string baseString = sequenceutil::indexToSequence(baseIndex, baseCount);
    sensitivityProfileFile << seqPos + 1 << "\t" << baseString << "\t" << sseProfile[pos] 
                           << "\t" << frequencyProfile[pos] 
                           << "\t" << sseProfile[pos] / frequencyProfile[pos]
                           << "\t" << sseRawProfile[pos]
                           << "\t" << sseRawProfile[pos] / frequencyProfile[pos] 
                           << "\t" << MathUtil::calculateMean(intensityVarianceProfile[pos])
                           << "\t" << MathUtil::calculateStandardDeviation(intensityVarianceProfile[pos])
                           << endl;
  }
  sensitivityProfileFile.close();
  
  // Save the sum of squares "profile" to file (per tuple)
  sensitivityProfileFile.open((filename + "_sum").c_str());
  for(size_t baseIndex = 0; baseIndex < sseProfileSum.size(); ++baseIndex) {
    size_t seqPos = 0;
    string baseString = sequenceutil::indexToSequence(baseIndex, baseCount);
    sensitivityProfileFile << seqPos << "\t" << baseString << "\t" << sseProfileSum[baseIndex] 
                           << "\t" << frequencyProfileSum[baseIndex] 
                           << "\t" << sseProfileSum[baseIndex] / frequencyProfileSum[baseIndex] 
                           << "\t" << probeCount << endl;
  }
  // Write total error into this file
  sensitivityProfileFile << "0\tTOTAL\t" << totalError 
                         << "\t" << probeCount 
                         << "\t" <<  totalError / probeCount
                         << "\t" << probeCount << endl;
  
  sensitivityProfileFile.close();
  
  
}
