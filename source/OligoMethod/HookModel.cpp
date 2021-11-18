/**
 * @file HookModel.cpp Defines the model used to describe microarray
 * intensities 
 * @author $Author: jan $ @date $Date: 2008-06-17 12:53:14
 * +0200 (Tue, 17 Jun 2008) $
 */

#include "Probe.hpp"
#include "HookModel.hpp"
#include "HookcurveAnalyzer.hpp"
#include "MathUtil.hpp"
#include "AveragingProcedure.hpp"
#include "StlUtil.hpp"
#include "Probe.hpp"

using namespace std;
using namespace larrpack;

HookModel::HookModel()
:mImax(), mSensitivityProfiles(kProfileTypeCount)
{
	initializeSensitivityProfiles();
}

/**
 * Initializes the sensitivity profiles
 *
 * @todo check if this is necessary

 */
void HookModel::initializeSensitivityProfiles()
{
  mSensitivityProfiles[kNsPm] = SensitivityProfilePtr(new SimpleSensitivityProfile(kProbeSequenceLength));
  mSensitivityProfiles[kNsMm] = SensitivityProfilePtr(new SimpleSensitivityProfile(kProbeSequenceLength));
  mSensitivityProfiles[kSFastDiff] = SensitivityProfilePtr(new SimpleSensitivityProfile(kProbeSequenceLength));
  mSensitivityProfiles[kSGcrmaDiff] = SensitivityProfilePtr(new SimpleSensitivityProfile(kProbeSequenceLength));
}


/**
 * Calculates the maximum intensity \f[ I_max \f]. We use the average of the
 * topX probes (only background corrected) ! with the highest averaged PM intensities.
 *
 * @param topX The topX largest intensities are used to compute Imax
 *
 * @return Maximum intensity.
 */
IntensityType HookModel::estimateImax(const ProbesetIntensitiesArray& probesetIntensitiesArray, const size_t topX)
{
  // Build a sorted list of intensities
  // vector<IntensityType> intensities(mProbes.size());
  // transform(mProbes.begin(), mProbes.end(), intensities.begin(), boost::mem_fn(&PmMmProbe::getPm)); 
  vector<IntensityType> intensities;

  // Get a list of all PMs
  for (ProbesetIntensitiesArray::const_iterator ps = probesetIntensitiesArray.begin();
       ps != probesetIntensitiesArray.end(); ++ps) {
    for (size_t probeIndex = 0; probeIndex < ps->size(); ++probeIndex) {
      intensities.push_back((*ps)[probeIndex]); 
    }
  }

  // Get a list of all MMs
/*  for (ProbesetIntensitiesArray::const_iterator ps = mCorrectedIntensities[kMmI]->begin();
       ps != mCorrectedIntensities[kMmI]->end(); ++ps) {
    for (size_t probeIndex = 0; probeIndex < ps->size(); ++probeIndex) {
      intensities.push_back((*ps)[probeIndex]); 
    }
  }*/
  sort(intensities.begin(), intensities.end());

  // Return average of topX elements
  IntensityArray maxIntensities(&intensities[intensities.size() - topX], topX);
  IntensitiesFunction mean = Mean<IntensityType>();
  
  //cout << "pre Imax: " << mean(maxIntensities) << endl;
  /*mImax = mean(maxIntensities);*/
  return mean(maxIntensities);
}


/**
 * Setter method. Sets a new (log) Imax value
 * 
 * @param iMax (Maximum intensity in log scale)
 */
void HookModel::setImax(IntensityType iMax)
{
  mImax = iMax;
}

/**
 * Getter method. Gets a new (log) Imax value
 * 
 */
const IntensityType& HookModel::getImax() const
{
  return mImax;
}


/**
 * Sets the value of the specific threshold (former intersection point) 
 * (The threshold is updated several times throughout the analysis of the microarray data)
 */
void HookModel::setNsThreshold(IntensityPair sThreshold)
{
	mThreshold = sThreshold;
}


/**
 * Returns the current current intensity threshold of specific hybridisation
 * The threshold is an intensity pair, this is due to the reason that in the hookcurve analysis
 * the thrshold point has a \Sigma and \Delta coordinate  
 * 
 * @return specific threshold (former intersection point)
 */
IntensityPair HookModel::getNsThreshold() const
{
  return mThreshold;
}



/**
 * Gives back the unsaturated intensity. Imax needs to be calculated before.
 *
 * @return The unlogged value of the unsaturated intensity.
 *
 * @param logSaturatedIntensity The saturated intensity as a log value.
 * @param cutOffMargin Desaturation of intensities close to Imax tend to
 *  distort for mathematical reasons. Thus, cutOffMargin specifies the 
 *  maximum desaturation value in percentage of Imax (ie. 0.99).
 * 
 */
IntensityType HookModel::getDesaturatedIntensity(const IntensityType logSaturatedIntensity, 
                                                 const IntensityType cutOffMargin) const
{
  // If the intensities get to close to Imax we obtain very high numbers 
  // since dividing by (almost) zero. We have to avoid
  // this by setting a margin, that is almost 1, but smaller
  IntensityType unloggedImax = exp10(mImax);
  IntensityType unloggedI = exp10(logSaturatedIntensity);
  IntensityType fi = 1.0/unloggedImax;
  fi  = fi * unloggedI;
  if (fi > cutOffMargin) {
    fi = cutOffMargin;
  }
  return (unloggedI  / (1.0 - fi));
}

/**
 * Returns the number of elements within a ProbesetIntensitiesArray (=the number of total probes)
 *
 *
 */
size_t countElementsInProbesetIntensityArray(const ProbesetIntensitiesArray& intensityArray)
{
  vector<size_t> probsetSizes(intensityArray.size());
  transform(intensityArray.begin(), intensityArray.end(), probsetSizes.begin(), boost::mem_fn(&std::valarray<IntensityType>::size));
  return accumulate(probsetSizes.begin(), probsetSizes.end(), 0);           
}

/**
 * Returns a list of values mark_i, i=0..n-1, where for each probe i
 *    mark_i = 1, if SumLogI( Probeset(i) ) > threshold
 *    mark_i = 0, if SumLogI( Probeset(i) ) <= threshold 
 *
 *
 * @param threshold The value that separates the absent vs. present probesets.
 *
 */
std::valarray<double> HookModel::markPresentProbes(const ProbesetIntensitiesArray& intensityArray, IntensityType threshold) const
{
  // Initialize list of marks with 0
  valarray<double> marks(0.0, countElementsInProbesetIntensityArray(intensityArray));
  
  // Loop over all probesets, remembering the index of the first probe 
  // of current probeset in mProbes
  size_t currentFirstProbeIndex = 0;
  // for(ProbesetVector::const_iterator psi = mProbesets.begin(); psi != mProbesets.end(); ++psi) {
  for (size_t probesetIndex = 0; probesetIndex < intensityArray.size(); ++probesetIndex) {
    
    // Set Mark to 1 if SumLogI > theshold
    //if (computeAverageIntensity(psi->getSumLogIs()) > threshold) {
//	  IntensityArray sumLog = sumLogIs[probesetIndex];
    //IntensityArray sumLogIs = getSumLogIs(kPmI, kMmI, probesetIndex); 
    if (computeAverageIntensity(intensityArray[probesetIndex]) > threshold) {  
      slice_array<double> specificityMarksOfCurrentSet = marks[slice(currentFirstProbeIndex, 
    		  intensityArray[probesetIndex].size(), 1)];
      specificityMarksOfCurrentSet = 1;
    }

    currentFirstProbeIndex += intensityArray[probesetIndex].size(); // psi->getSize();
  }

  return marks;
}

/**
 * Calculates the PM and MM increments for a whole probeset.
 * The results are returned in a valarray.
 *
 * @param intensities Uncorrected probe intensities 
 * @param ps The probeset to calculate the increments for.
 * @param calculationType Selects how to handle the profiles.
 *
 * @note If specific profiles are available, the compression factor 
 * is always beeing used (-> don't use this method for final 
 * expression measures)
 *
 * @return increment
 */
IntensityArray HookModel::calculateIncrements(const Probeset& probeset, const IntensityArray& intensities,
                                              ProbeType probeType, IntensityType probesetSumLogI) const
{
  // Get specific accessor types
  ProfileType nsProfile = probeTypeParams[probeType].nsProfile;
  ProfileType sProfile = probeTypeParams[probeType].sProfile;
  PmMmProbeSequenceFunction getSequence  = probeTypeParams[probeType].getSequence;
  
  // If we use the average specific / ns fractions for each probeset, then we need to save
  // the calculated increments to a valarray which can be accessed from the probeset level.
  IntensityArray probesetIncrementsSpecific(0.0, probeset.getSize());
  IntensityArray probesetIncrementsNonspecific(0.0, probeset.getSize());
  IntensityArray increment(0.0, probeset.getSize());

  // For each probe in the set
  for (size_t probeIndex = 0; probeIndex < probeset.getSize(); ++probeIndex) {
    PmMmProbePtr currentProbe = probeset.getProbePtr(probeIndex);

    // Calculate the Nonspecific increments for the current probe
    probesetIncrementsNonspecific[probeIndex] = getProfile(nsProfile)->getSequenceIncrement(getSequence(*currentProbe));
    if (getProfile(sProfile)) { // if specific profile available
      probesetIncrementsSpecific[probeIndex] = getProfile(sProfile)->getSequenceIncrement(getSequence(*currentProbe));
    } // If kAllProfiles     
  } // Probe iterator

  // Calculate the increment for each probe in the probeset by
  // taking into account the fractions of specific to nonspecific portions.
  // The fractions were calculated as a mean value of each probeset.    
  if (getProfile(sProfile)) { // if specific profile available
    IntensityType nonspecificFraction = getNsFraction(probesetSumLogI);

    increment = ((1 - nonspecificFraction) * probesetIncrementsSpecific) 
      + nonspecificFraction * probesetIncrementsNonspecific; 
  } 
  else {
    // In case of kNonspecificProfilesOnly the ns-fraction is 100% 
    increment = probesetIncrementsNonspecific;
  }
  
  return increment;
}


/**
 * Returns the compression to the calculated sensitivity of a probe.
 * 
 * @param logIntensity The logarithm of a probe intensity
 * @return The compression of the probes sensitivity
 * 
 */
IntensityType HookModel::getCompressionFactor(IntensityType logIntensity) const
{
  // mAp and mBp are the compression factor parameters 
  return (getParam("mAp") + (getParam("mBp") * logIntensity)) / logIntensity; 
}
