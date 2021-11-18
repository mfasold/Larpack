/**
 * @file SubsetSensitivityProfile.hpp Defines methods for sequence profiles of a subset of nucleotides
 * @author $Author $ 
 * @author Mario Fasold
 * @date $Date: 2008-04-17 10:12:26 -0400 (Thu, 17 Apr 2008) $
 */
#include <fstream>
#include <numeric>
#include <boost/mem_fn.hpp>
#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>

#include "MathUtil.hpp"
#include "SubsetSensitivityProfile.hpp"
#include "ValarrayUtil.hpp"

#include <iterator>


using namespace std;
using namespace larrpack;
using namespace sequenceutil;

/**
 * Constructor. Initializes with an empty profile (zero increment).
 *
 * @param relevantBasetuples The base tuples (e.g. ["GGGG","CCCC"]) 
 *        the contributions of which the profile shall model 
 * @param sequenceLength Length of the sequence.
 * @param modelRank Size of the base tuples within the current model
 * (mononucleotides = 1, dinucleotides = 2, ...)
 * 
 */
SubsetSensitivityProfile::SubsetSensitivityProfile(std::vector<std::string> relevantBasetuples,
                                                   size_t sequenceLength, 
                                                   unsigned short modelRank) 
  : SimpleSensitivityProfile(sequenceLength, modelRank, 
                       sequenceutil::getValidSequenceCount(sequenceLength, modelRank) * relevantBasetuples.size()),
    mRelevantBasetuples(relevantBasetuples)
{}

/**
 * Returns the sequence specific increment 
 * \f$ \delta S = \sum_{k=1}^{sequenceLength} \sigma_k( sequence_k ) \f$
 *
 * @param sequence The nucleotide sequence.
 */
IntensityType SubsetSensitivityProfile::getSequenceIncrement(const std::string& sequence) const
{
  assert(sequence.size() == mSequenceLength);

  IntensityType sensitivity = 0.0;
  // Look for all tuples
  for (size_t tupleIndex = 0; tupleIndex < mRelevantBasetuples.size(); ++tupleIndex) {
    size_t beginIndex = 0;
    // ...on all sequence positions
    while (1) {
      size_t findIndex = sequence.find(mRelevantBasetuples[tupleIndex], beginIndex);
      if (findIndex == string::npos) {
        break;
      }

      // Tuple found -> add contribution
      sensitivity += mProfile[tupleIndex*getValidSequenceCount() + findIndex];  //@note Knowledge about profile design used
      beginIndex = findIndex + 1;
    }
  }
  return sensitivity;
}

/**
 * Calculates the sensitivity values \f$ \sigma_k \f$ as. 
 * The function first computes matrix (12) and then solves it using SVD (with
 * zeroing).
 *
 * @param probesetSequences List of probe sequences for all probesets, 
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
 * @note This function is 95% copied from SensitivityProfile (evil).
 * @todo Integrate with SensitivityProfile.computeProfile.
 *
 */
void SubsetSensitivityProfile::computeProfile(const ProbesetSequencesArray& probesetSequences,
                                              const ProbesetIntensitiesArray& probesetIntensitiesArray,
                                              const size_t probesetLimit)
{

  // delta * sigma = theta
  // Initialize matrix delta and valarray theta with zeros.
  // the matrix and the valarray will be used to calculate the valarray sigma by singular value decomposition.
  const size_t profileSize = mProfile.size(); // getProfileSize(mModelRank, kProbeSequenceLength, kNucleotideCount);
  math::Matrix<double> delta(profileSize, profileSize, 0.0);
  valarray<double> theta(0.0, profileSize);

  // Define some counters for fast access
  const size_t modelNucleotides = getBasetupleCount();
  const size_t validSequencePositions = getValidSequenceCount();

  cout << "Selective Sequence Profile of " << modelNucleotides;
  cout << " tuple(s) and a total profile size of " << profileSize << "." << endl;

  // Iterate over probesets:
  size_t counter  = 0; // This is to draw progress points
  size_t counter2 = 0; // This counter is to limit the number probesets for the profile calculation.
  size_t ithProbeset;
  if (probesetLimit > 0) {
    ithProbeset = probesetSequences.size() / probesetLimit;
    if (ithProbeset == 0) {
      ithProbeset = 1;
    }
  }
  else { ithProbeset = 1; }
  if (ithProbeset > 1) {
    cout << "number of probesets: " << probesetSequences.size() << " we take every " << ithProbeset << "th for sequenceprofile calculation." << endl;
    cout << "That means ~ " <<  probesetSequences.size()/ithProbeset << " probesets." << endl;
  }
  
  for (size_t probesetIndex = 0; probesetIndex < probesetSequences.size(); ++probesetIndex) {
    const StringVector& currentSequences = probesetSequences[probesetIndex];

    

    //   for (ProbesetVector::const_iterator currentProbeset = probesets.begin(); currentProbeset != probesets.end(); ++currentProbeset){
    if (! (counter2 % ithProbeset == 0)){
      counter2++;
      continue;
    }
    counter2++;
    if ((counter % 1000) == 0) {
      cout << "." << flush;
    }
    counter++;

    //     cout << "profileSize" << endl;
    //     cout << "validSequencePositions" << endl;
    //     cout << "mModelRank: " << mModelRank << endl;
    //     cout << "getBasetupleCount(): " << getBasetupleCount() << endl;
    //     cout << "sequenceToIndex(mRelevantBasetuples[0]): " << sequenceToIndex(mRelevantBasetuples[0]) << endl;
    // @note difference 1 to SensitivityProfile.computeProfile ----------------------------------------
    // adapt PSWM to new base indexes - just copy the required rows
    math::Matrix<double> originalProbesetPSWM = 
      sequenceutil::getPositionSpecificWeightMatrix(currentSequences, mModelRank);

    math::Matrix<double> probesetPSWM(getBasetupleCount(), mSequenceLength);
    for (size_t i = 0; i < getBasetupleCount(); ++i) {
      for (size_t col = 0; col < mSequenceLength; ++col) {
        probesetPSWM[i][col] = originalProbesetPSWM[sequenceToIndex(mRelevantBasetuples[i])][col];
      }    
    }
    // difference 1 end
    
    // This is not nicely done... we need to filter out all probes smaller 2 (they might produce problems in the log scale)
    size_t count = 0; // Count saves how many probes of a probeset are > 2
    for (size_t probeIndex = 0; probeIndex < currentSequences.size(); ++probeIndex)
      if (MathUtil::gExp10(probesetIntensitiesArray[probesetIndex][probeIndex]) > 2)
        count++;
    IntensityArray intensitiesTmp(0.0, count); // IntensitiesTmp holds all intensities with unlogged values > 2
    size_t idx = 0;
    for (size_t probeIndex = 0; probeIndex < currentSequences.size(); ++probeIndex){
      if (MathUtil::gExp10(probesetIntensitiesArray[probesetIndex][probeIndex]) < 2)
        continue;
      else {
        intensitiesTmp[idx] = probesetIntensitiesArray[probesetIndex][probeIndex];
        idx++;
      }
    }
    if (intensitiesTmp.size() == 0){ // If all probes of the probeset have intensities < 2 quit for this probeset,
      continue; // Leave the probeset loop
    }
    
    IntensityType average = computeAverageIntensity(intensitiesTmp);
    //     IntensityType average = computeAverageIntensity(getProbeIntensities(*currentProbeset));
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
    for (size_t probeIndex = 0; probeIndex < currentSequences.size(); ++probeIndex) {
      if (MathUtil::gExp10(probesetIntensitiesArray[probesetIndex][probeIndex]) < 2){
        continue;
      }  
      string sequence = currentSequences[probeIndex];
      
      // First we fill up this valarray, this will save runtime in the end.
      IntensityArray kroeMinusP(0.0, validSequencePositions * modelNucleotides);
      for (size_t pos = 0; pos < validSequencePositions; ++pos) {
        // Set kroenecker of the only one occuring base tuple to one
        // @note difference 2 to SensitivityProfile.computeProfile ---------------------------------------
        size_t tupleIndexBase = sequenceToSubsetIndex(sequence.substr(pos, mModelRank)); 
        if (tupleIndexBase < getBasetupleCount()) {
          kroeneckers[tupleIndexBase] = 1;
        }
        // difference 2 end       

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
          theta[base1 * validSequencePositions + pos1] += yExp[probeIndex] * kroeMinusP[base1 * validSequencePositions + pos1];
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

  }//for probesets
  cout << endl;
  
  // Calculate upper triangle of the matrix
  for (size_t row = 0; row < profileSize - 1; ++row) {
    for (size_t column = row; column < profileSize; ++column) {
      delta[row][column] = delta[column][row];
    }
  }
  // ofstream logFile;
  // logFile.open("affinityDelta2.txt");
  // logFile << delta << endl;
  // logFile.close();


  // logFile.open("affinityTheta2.txt");
  // logFile <<  theta << endl;
  // logFile.close();

  // Solve equation
  mProfile = MathUtil::solveEquationsSvd(delta, theta/*, 1.0*/);
}


/**
 * Constructor. Initializes employed profiles.
 *
 * @param relevantBasetuples The base tuples (e.g. ["GGGG","CCCC"]) 
 *        the contributions of which the profile shall model 
 * @param profiles Profiles that build up the new composite profile.
 */
AdditivePofile::AdditivePofile(std::vector<std::string> relevantBasetuples, 
                               std::vector<SensitivityProfilePtr> profiles) 
  : mRelevantBasetuples(relevantBasetuples), mProfiles(profiles)
{}

/**
 * Constructor. Initializes employed profiles.
 *
 * @param relevantBasetuples The base tuples (e.g. ["GGGG","CCCC"]) 
 *        the contributions of which the profile shall model 
 * @param sequenceLength Length of the sequence.
 * @param modelRank Size of the base tuples within the current model
 * (mononucleotides = 1, dinucleotides = 2, ...)
 */
AdditivePofile::AdditivePofile(std::vector<std::string> relevantBasetuples,
                               size_t sequenceLength, 
                               unsigned short modelRank)
  : mRelevantBasetuples(relevantBasetuples)
{
  // Check if basetuples have same size
  vector<size_t> tupleSizes;
  transform(relevantBasetuples.begin(), relevantBasetuples.end(), back_inserter(tupleSizes),
            boost::mem_fn(&string::size));
  assert(*min_element(tupleSizes.begin(), tupleSizes.end()) == *max_element(tupleSizes.begin(), tupleSizes.end()));

  // Create both profiles and add to vector
  SensitivityProfilePtr profileNSMm(new SimpleSensitivityProfile(sequenceLength, modelRank));  
  SensitivityProfilePtr profileGGG(new SubsetSensitivityProfile(relevantBasetuples, sequenceLength, 
                                                                tupleSizes[0]));
  mProfiles.push_back(profileNSMm);
  mProfiles.push_back(profileGGG);
}

// string getSeq(const PmMmProbePtr& p) 
// {
//   return p->getSequence() + "\n";
// }

/**
 * Computes all profiles incrementally
 *
 * @param probesetSequences List of probe sequences for all probesets, 
 *        used to determine the profile.
 * @param probesetIntensitiesArray List of probe intensities for all probesets
 * @param probesetLimit Number of probesets to be used for profile computation.
 *        They are sampled ot of all probesets.
 *
 */
void AdditivePofile::computeProfile(const ProbesetSequencesArray& probesetSequences,
                                    const ProbesetIntensitiesArray& probesetIntensitiesArray,
                                    const size_t probesetLimit)
{
  ProbesetIntensitiesArray currentIntensities = probesetIntensitiesArray; // @todo: Why copy??
  //ProbesetVector probesetsWithTuple, probesetsNoTuples;
  ProbesetSequencesArray probesetSequencesWithTuple, probesetSequencesNoTuples;
  ProbesetIntensitiesArray intensitiesWithTuple, intensitiesNoTuples;

  cout << " Probesets in total: " << probesetSequences.size() << endl;
  // Divide probesets into two sets, either containing or not containg the tuples
  ContainsBasetuples containsBasetuples(mRelevantBasetuples);
  for (size_t probesetIndex = 0; probesetIndex < probesetSequences.size(); ++probesetIndex) {
    if (containsBasetuples(probesetSequences[probesetIndex])) {
      // transform(probesetSequences[probesetIndex].begin(), probesetSequences[probesetIndex].end(),
      //                 ostream_iterator<string>(cout), getSeq2);
      //       // boost::mem_fn(&PmMmProbe::getSequence));      
      //       cout << endl;
      probesetSequencesWithTuple.push_back(probesetSequences[probesetIndex]);
      intensitiesWithTuple.push_back(probesetIntensitiesArray[probesetIndex]);
    }
    else {
      probesetSequencesNoTuples.push_back(probesetSequences[probesetIndex]);
      intensitiesNoTuples.push_back(probesetIntensitiesArray[probesetIndex]);
    }
  }

  // Compute simple profile for sets not containing the tuples
  cout << "Computing additive profile 0, probesets: " << probesetSequencesNoTuples.size() << endl;
  mProfiles[0]->computeProfile(probesetSequencesNoTuples, intensitiesNoTuples, probesetLimit); 
  cout << "done: " << endl;

  // Compute other profiles using the probes containing the tuples
  for (vector<SensitivityProfilePtr>::iterator currentProfile = mProfiles.begin() + 1; 
       currentProfile != mProfiles.end(); ++currentProfile) {
    // Subtract increments of previous profile from corrected intensities 
    for (size_t probesetIndex = 0; probesetIndex < probesetSequencesWithTuple.size(); ++probesetIndex) {
      for (size_t probeIndex = 0; probeIndex < probesetSequencesWithTuple[probesetIndex].size(); ++probeIndex) {
        intensitiesWithTuple[probesetIndex][probeIndex] -= (*(currentProfile - 1))
          ->getSequenceIncrement(probesetSequencesWithTuple[probesetIndex][probeIndex]);
      }
    }
    // Compute profile with current corrected intensities
    cout << "Computing other additive profile, probesets: " << probesetSequencesWithTuple.size() << endl;    
    (*currentProfile)->computeProfile(probesetSequencesWithTuple, intensitiesWithTuple, probesetLimit);
  }
}


/**
 * Returns the sum over each member profile's increment
 *
 * @param sequence Nucleotide sequence
 */
IntensityType AdditivePofile::getSequenceIncrement(const std::string& sequence) const
{
  vector<IntensityType> increments;
  transform(mProfiles.begin(), mProfiles.end(), std::back_inserter(increments), // Compute all increments
            boost::bind(&SensitivityProfile::getSequenceIncrement, _1, sequence));
//   cout << sequence  << "\t" <<  mProfiles[0]->getSequenceIncrement(sequence) << "\t" << 
//     mProfiles[1]->getSequenceIncrement(sequence) << "\t" <<  accumulate(increments.begin(), increments.end(), 0.0) << endl;
  return accumulate(increments.begin(), increments.end(), 0.0);
  // return mProfiles[0]->getSequenceIncrement(sequence);
//   for (std::vector<std::string>::const_iterator tuple = mRelevantBasetuples.begin();
//        tuple != mRelevantBasetuples.end(); ++tuple) {
//     if (sequence.find(*tuple) != std::string::npos) {
//       return 10;
//     }
//   }
//   return  -1;
}

/**
 * Returns the sigma for a base/oligonucleotide index at a given position.
 *
 * @note Different profiles might use different indices and have different 
 *       lengths. For simplicity, we only return the sigma of the first profile
 * 
 * @param baseIndex Numerical value representing the base/oligonucleotide
 * @param position Position in the probe
 *
 * @return Sigma
 */
IntensityType AdditivePofile::getSigma(size_t baseIndex, size_t position) const
{
//   vector<IntensityType> sigmas;
//   transform(mProfiles.begin(), mProfiles.end(), std::back_inserter(sigmas), // Compute all sigmas
//             boost::bind(&SensitivityProfile::getSigma, _1, baseIndex, position));
//   return accumulate(sigmas.begin(), sigmas.end(), 0.0);
  return mProfiles[0]->getSigma(baseIndex, position);
}

/**
 * Returns the length of the profile 
 *
 * @note Different profiles might use different indices and have different 
 *       lengths. For simplicity, we only refer to the the first profile. 
 * @see AdditivePofile::getSigma
 */
size_t AdditivePofile::getValidSequenceCount() const
{
  return mProfiles[0]->getValidSequenceCount();
}

/**
 * Sets the effect of the middle base of all profiles to zero. 
 *
 */
void AdditivePofile::zeroMiddlebase()
{    
  for_each(mProfiles.begin(), mProfiles.end(), boost::mem_fn(&SensitivityProfile::zeroMiddlebase));
}


/**
 * Exports the profile to a file in tabular text format.
 *
 * @param filename Filename.
 */
void AdditivePofile::exportToDatafile(const std::string filename) const
{
  string filenameEdit = filename;
  filenameEdit.insert(filenameEdit.rfind("."), "0");
  for (size_t i = 0; i < mProfiles.size(); ++i) {    
    filenameEdit.replace(filenameEdit.rfind(".") - 1, 1, boost::lexical_cast<std::string>(i));
    mProfiles[i]->exportToDatafile(filenameEdit);
  }
}
/// @bug Rank changes but basetuple sizes stay the same
SensitivityProfilePtr AdditivePofile::cloneWithNewRank(const size_t modelRank) const
{
  std::vector<SensitivityProfilePtr> newProfiles;
  transform(mProfiles.begin(), mProfiles.end(), back_inserter(newProfiles),
            boost::bind(&SensitivityProfile::cloneWithNewRank, _1, modelRank));
  SensitivityProfilePtr newProfile = SensitivityProfilePtr(new AdditivePofile(mRelevantBasetuples, newProfiles));
  return newProfile;
}

