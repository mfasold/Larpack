/**
 * @file ExpressionMeasure.cpp Defines methods for calculating the final Expression Measures
 * @author $Author: mario $ 
 * @author Jan Bruecker
 * @date $Date: 2008-08-29 16:21:54 +0200 (Fri, 29 Aug 2008) $
 */

#include <ExpressionMeasure.hpp>
#include <boost/lexical_cast.hpp>
#include "HookModel.hpp"

using namespace std;

namespace larrpack
{

ExpressionMeasure::ExpressionMeasure(const Chip& chip, // const
        HookModel& hookModel/*,
        HookcurveAnalyzer& hookcurveModel*/)
: mExpressionMeasures(kExpressionMeasureTypeCount),
  mExpressionMeasuresCalculated(kExpressionMeasureTypeCount, false),
  mChip(chip),
  mHookModel(hookModel)/*,
  mHookcurveModel(hookcurveModel)*/
{
}

//ExpressionMeasure::~ExpressionMeasure()
//{
//}
//
//}


/**
 * Initializes the expression measure
 * 
 */
void ExpressionMeasure::initializeExpressionMeasure(ExpressionMeasureType t) {
  mExpressionMeasures[t] = Probeset::initializeIntensitiesArray(mChip.getProbesets());
  mExpressionMeasuresCalculated[t] = true;
}




/**
 * Calculates S^{PM}_{Set} using the distributions of the NS PM probes.
 *
 * @param nsSubtractedIntensities Holds the log intensities desaturated and substracted by their nonspecific content
 * @param profileType Identifier to the SensitivityProfile that holds the specific profiles
 * @param getProbeSequence Function that gives back the sequence of a PmMmProbe (either the MM or the Pm sequence)
 * @param invLogarithm Function that calculates the inverse Logarithm. (Can either be the invGlog or the invLog)
 * @param expMeasureType If an ExpressionMeasureType is given, the results are also printed to an
 * @return   Vector of Probeset Expression pairs ( A pair holds the probeset identifier and the expression value)
 * ProbesetIntensitiesArray that is accessible by the mCorrectedIntensities vector.
 */
std::vector<ProbesetExpressionPair> ExpressionMeasure::calculateProbesetExpressionMeasure(const ProbesetIntensitiesArray& nsSubtractedIntensities, /*const IntensityMode intensityMode,*/
                                                                                 const ProfileType profileType,
                                                                                 const PmMmProbeSequenceFunction& getProbeSequence,
                                                                                 UnaryIntensityFunction& invLogarithm,
                                                                                 ExpressionMeasureType expMeasureType) /*const*/
{  
  bool writeToArray = false;
  if (expMeasureType != kExpressionMeasureTypeCount) {
    initializeExpressionMeasure(expMeasureType);
    writeToArray = true;
  }
  std::vector<ProbesetExpressionPair> expressionVector; 
  for (size_t probesetIndex = 0; probesetIndex < mChip.getProbesets().size(); ++probesetIndex) {
    const Probeset& currentProbeset = mChip.getProbesets()[probesetIndex];
    IntensityArray probesetExpressionValues(0.0, currentProbeset.getSize());
    for (size_t probeIndex = 0; probeIndex < currentProbeset.getSize(); ++probeIndex) {
      PmMmProbePtr currentProbe = currentProbeset.getProbePtr(probeIndex);
      IntensityType probesetIncrementSpecific;
      probesetIncrementSpecific = mHookModel.getProfile(profileType)->getSequenceIncrement(getProbeSequence(*currentProbe));
      probesetExpressionValues[probeIndex] = nsSubtractedIntensities[probesetIndex][probeIndex] - probesetIncrementSpecific;
    }
    // If you want reasonable values for affycomp use exp10 instead of the inverse glog.
//     if (useExp10) {
//       expressionVector.push_back(ProbesetExpressionPair(currentProbeset.getProbesetId(), exp10(computeAverage(probesetExpressionValues))));
//     }
//     else {
   IntensityType expMeasure = invLogarithm(computeAverageIntensity(probesetExpressionValues));
   expressionVector.push_back(ProbesetExpressionPair(currentProbeset.getProbesetId(), expMeasure));
   if (writeToArray) {
     for (size_t probeIndex = 0; probeIndex < currentProbeset.getSize(); ++probeIndex) {
       (*mExpressionMeasures[expMeasureType])[probesetIndex][probeIndex] = expMeasure;
     }
   }
//     }
  }
  return expressionVector;
}





/**
 * Calculates S^{PM} for each probe using the distributions of the NS PM probes.
 * 
 * @param nsSubtractedIntensities Holds the nonspecific background subtracted log intensities
 * @param profileType Identifier to the SensitivityProfile that holds the specific profiles
 * @param getProbeSequence Function that gives back the sequence of a PmMmProbe (either the MM or the Pm sequence)
 * @param invLogarithm Function that calculates the inverse Logarithm. (Can either be the invGlog or the invLog)
 * @param expMeasureType If an ExpressionMeasureType is given, the results are also printed to an
 * 
 * @return vector holding the expression values (and their identifier)
 *   ProbesetIntensitiesArray that is accessible by the mCorrectedIntensities vector.
 */
std::vector<ProbesetExpressionPair> 
ExpressionMeasure::calculateProbeExpressionMeasures(const ProbesetIntensitiesArray& nsSubtractedIntensities, /*const IntensityMode intensityMode,*/
                                                  const ProfileType profileType,
                                                  const PmMmProbeSequenceFunction& getProbeSequence,
                                                  UnaryIntensityFunction& invLogarithm,
                                                  ExpressionMeasureType expMeasureType) /*const*/
{  
  bool writeToArray = false;
  if (expMeasureType != kExpressionMeasureTypeCount) {
    initializeExpressionMeasure(expMeasureType);//mExpressionMeasures[expMeasureType] = Probeset::initializeIntensitiesArray(mProbesets);
    writeToArray = true;
  }
  std::vector<ProbesetExpressionPair> expressionVector;
 
  for (size_t probesetIndex = 0; probesetIndex < mChip.getProbesets().size(); ++probesetIndex) {
    const Probeset& currentProbeset = mChip.getProbesets()[probesetIndex];
    
//    cout << "size of probeset " << currentProbeset.getSize() << endl;
    for (size_t probeIndex = 0; probeIndex < currentProbeset.getSize(); ++probeIndex) {
      PmMmProbePtr currentProbePtr = currentProbeset.getProbePtr(probeIndex);

      // Compute expression measure by subtracting (specific) sequence dependent incremnt
//      cout << mHookcurveModel.getLogIs(intensityMode, probesetIndex)[probeIndex] - mHookModel.getProfile(profileType)->getSequenceIncrement(getProbeSequence(*currentProbePtr));
      IntensityType expressionMeasure = nsSubtractedIntensities[probesetIndex][probeIndex] /*mHookcurveModel.getLogIs(intensityMode, probesetIndex)[probeIndex]*/
        - mHookModel.getProfile(profileType)->getSequenceIncrement(getProbeSequence(*currentProbePtr));

      // Construct probeId from probesetId (if available, i.e. expression array) 
      // or from chromosome and location (tiling)
      std::string probeId = currentProbeset.getProbesetId();
      if (probeId.empty()) {
        probeId = (*currentProbePtr).getChromosome() + "\t" 
          + boost::lexical_cast<std::string>((*currentProbePtr).getLocation()); 
      }

      IntensityType expMeasure = invLogarithm(expressionMeasure);
      expressionVector.push_back(ProbesetExpressionPair(probeId, expMeasure));
      if (writeToArray) {
        (*mExpressionMeasures[expMeasureType])[probesetIndex][probeIndex] = expMeasure;
      }
    }
  }
  return expressionVector;
}

}


