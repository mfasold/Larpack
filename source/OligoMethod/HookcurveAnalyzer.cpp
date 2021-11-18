/**
 * @file HookcurveAnalyzer.cpp Defines methods for analysis of PmMmProbes and -probesets.
 * @author $Author$ 
 * @author Mario Fasold
 * @date $Date$
 */
#include <string>
#include <vector>
#include <iterator>

#include <iostream>
#include <fstream>
#include <stdexcept>

#include <functional>
#include <algorithm>

#include <gsl/gsl_integration.h>
#include <boost/lexical_cast.hpp>
#include <boost/mem_fn.hpp>
#include "boost/tuple/tuple.hpp"
#include <boost/progress.hpp>
#include <boost/bind.hpp>
#include <boost/any.hpp>

#include "StringUtil.hpp"
#include "StlUtil.hpp"
#include "MathUtil.hpp"
#include "ValarrayUtil.hpp"
#include "SequenceUtil.hpp"        
#include "AveragingProcedure.hpp"

#include "HookcurveAnalyzer.hpp"
#include "Probe.hpp"
#include "PmMmProbe.hpp"
#include "ProbeFilter.hpp"

#include "HookModel.hpp"
#include "ExpressionMeasure.hpp"


#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace larrpack;
using namespace sequenceutil;
// using namespace boost;


/**
 * Constructor that creates the internal probesets using the
 * ProbesetCompositionFunction. 
 * Initializes default (zero increment) sensitivity profile.
 *
 * @param probes List of PM/MM-probes
 * @param composeProbesetFunction Function object to create probesets composition.
 */
HookcurveAnalyzer::HookcurveAnalyzer(const Chip& chip, // const
                                     HookModel& hookModel) 
  : mProbes(chip.getProbes()), 
    mProbesets(chip.getProbesets()),
    mAp(), mBp(MathUtil::getNaN()),
    mChip(chip),
    mHookModel(hookModel)
{
}



// Sorting helper functions
bool isSecondElementSmaller(IntensityPair a, IntensityPair b)
{
  return a.second < b.second;
}

/**
 * Checks the chip for saturation and sets the respective member.
 * @param graph Probesets SumLogI and DeltaLogIs
 * @param intervals Number of intervals the curve shall be digitized to
 * 
 */
bool HookcurveAnalyzer::isSaturated(const IntensityMapping& graph, size_t intervals )
{
  vector<IntensityPair> plot = MathUtil::digitizeCurveRobust(graph , intervals);

  // Get the index of maximum element
  size_t maxIndex = distance(plot.begin(), 
                             max_element(plot.begin(),plot.end(),isSecondElementSmaller));

  // Consider the graph saturating if the maximum is not among the last 3 elements
  if (maxIndex < intervals - 3) {
    return true;
  }
  return false;
}

  

/// @todo Remove or use just one (same function is found in SensitivityProfile)
int _kroneck(size_t base, size_t baseAtPosition){
  if (base == baseAtPosition) {
    return 1;}
  else {
    return 0;}
}

void printM(const math::Matrix<double>& m) {
  for(size_t i = 0; i < m.columnNum(); ++i)  {
    for(size_t j = 0; j < m.rowNum(); ++j)
      cout << '\t' << m(i,j);
    cout << endl;
  }
  cout << endl;
}


/**
 * Calculates the compression factor parameters mAp and mBp 
 * 
 * @param probesetIntensitiesArray The data the compression shall be calculated from.
 * @param getProbeSequence A function that returns the sequence of a probe
 * 
 * @todo : Use the sensitivity corrected intensities.
 *
 */
void HookcurveAnalyzer::calculateCompressionFactorParameters(const ProbesetIntensitiesArray& probesetIntensitiesArray, const PmMmProbeSequenceFunction& getProbeSequence)
{
  double dMaxNs = 1;                 // By definition 1
  IntensityType iMaxNs = mHookModel.getNsThreshold().first; // The intensity value at the intersection point
  // Since PM and MM are approximately the same at that point
  // We can directly use the respective intensity value.
  //estimateImax(); // is done befor in main or OligoController
  const IntensityType& iMax = mHookModel.getImax(); // The maximal intensity
      
  
  // We need a single nucleotide profile to calculate the compression. so check if we already have one,
  // otherwise we create it.(that's exactly what the following constructor does:
  SensitivityProfilePtr sensitivityProfileNSPmSingle = mHookModel.getProfile(kNsPm)->cloneWithNewRank(1);
  size_t modelRank = 1; // We now have a single nucleotide profile
  size_t validSequencePositions = sensitivityProfileNSPmSingle->getValidSequenceCount();
  //sensitivityProfileNSPmSingle->exportToDatafile("singlePM.profile");
 
  // Matrix delta and vector theta
  math::Matrix<IntensityType> delta(4, 4, 0.0);
  IntensityArray theta(0.0, 4);
  
  // Only the 5% probesets wit maximal average PM intensity are used
  vector<IntensityType> intensities;
  for (ProbesetIntensitiesArray::const_iterator ps = probesetIntensitiesArray.begin();
       ps != probesetIntensitiesArray.end(); ++ps) {
    intensities.push_back(computeAverageIntensity(*ps));
  }
  sort(intensities.begin(), intensities.end());
  // The 5% probesets with the highest average intensity are used. the lower bound is given by a cutOff.
  IntensityType cutOff = intensities[((size_t) (0.95 * intensities.size()))];
 
  for (size_t probesetIndex = 0; probesetIndex < mProbesets.size(); ++probesetIndex) {
   const IntensityArray& currentPmIs = probesetIntensitiesArray[probesetIndex];
   IntensityType averagePM = computeAverageIntensity(currentPmIs); // currentProbeset->getPmLogIs()
    //if (probesetIndex == 0) cout << averagePM << endl; 
    if (averagePM > cutOff){     
      //IntensityArray intensitiesPM = currentPmIs; // currentProbeset->getPmLogIs(); /// @todo: replace intensitiesPM with currentPmIs
      IntensityArray yExp = currentPmIs - averagePM;
//      if(0){ // As long as we don't substract frequencies here: use the faster loop (else if)
//             // and if this should ever be implemented take a look at computeSequenceProfiles for how thsi can be done efficently.
//        for (size_t probeIndex = 0; probeIndex < currentPmIs.size(); ++probeIndex) {
//          string sequence = getProbeSequence(mProbesets[probesetIndex].getProbe(probeIndex));
//          
//          // The quadroloop fills the entries of vector theta and matrix delta
//          for (size_t baseIndex = 0; baseIndex < 4; ++baseIndex){
//            for (size_t posIndex = 0; posIndex < validSequencePositions; ++posIndex){
//              theta[baseIndex] += yExp[probeIndex] * sensitivityProfileNSPmSingle->getSigma(baseIndex, posIndex)  
//                *  _kroneck(baseIndex,baseToIndex(sequence[posIndex]));
//              for (size_t baseIndex2 = 0; baseIndex2 < 4; ++baseIndex2){
//                for (size_t posIndex2 = 0; posIndex2 < validSequencePositions; ++posIndex2){
//                  delta[baseIndex][baseIndex2] +=   sensitivityProfileNSPmSingle->getSigma(baseIndex, posIndex) 
//                    * _kroneck(baseIndex,baseToIndex(sequence[posIndex]))
//                    * sensitivityProfileNSPmSingle->getSigma(baseIndex2, posIndex2) 
//                    * _kroneck(baseIndex2,baseToIndex(sequence[posIndex2]));
//                }// For posIndex2
//              }// For baseIndex2
//            }// For posIndex
//          }// ForbaseIndex
//        }// For probes
//      }// If
      
//      else if(1) { // The faster version which does not use frequencies.
        for (size_t probeIndex = 0; probeIndex < currentPmIs.size(); ++probeIndex) {
          string sequence = getProbeSequence(mProbesets[probesetIndex].getProbe(probeIndex));
          for (size_t posIndex = 0; posIndex < validSequencePositions; ++posIndex) {
            theta[sequenceToIndex(sequence.substr(posIndex, modelRank))] 
              += yExp[probeIndex] * sensitivityProfileNSPmSingle->getSigma(sequenceToIndex(sequence.substr(posIndex, modelRank)), posIndex);
            for (size_t posIndex2 = 0; posIndex2 < validSequencePositions; ++posIndex2) {
              delta[sequenceToIndex(sequence.substr(posIndex, modelRank))]
                [sequenceToIndex(sequence.substr(posIndex2, modelRank))] 
                += sensitivityProfileNSPmSingle->getSigma(sequenceToIndex(sequence.substr(posIndex, modelRank)), posIndex) 
                * sensitivityProfileNSPmSingle->getSigma(sequenceToIndex(sequence.substr(posIndex2, modelRank)), posIndex2);
            } //pos2
          } // pos1
        } // For probes
//      }// Else if
      
    } // If > cutOff
  }// For probesets
  
  // The P(A), P(C), P(G) and P(T) defined by the linear system of equations. :
  IntensityArray p = MathUtil::solveEquationsSvd(delta, theta);
  
  // Compute the compression at maximum intensity as in 
  // Preibisch page 48: min(d) = dImax = (P(A) + P(C)) / 2
  IntensityType dImax = (p[0] + p[1]) / 2;   

  // Calculate Ap and Bp and save it to their member variables.
  // (Since we now know to pairs of values for I^P_p and d(I^P_p) we can solve 
  // equation 23 in Preibisch to calulate A^P and B^P.
  mAp = ((dMaxNs - dImax) * iMaxNs * iMax) / (iMax - iMaxNs);
  mBp = dMaxNs - (mAp / iMaxNs);

  mHookModel.setParam(string("mAp"), mAp);
  mHookModel.setParam("mBp", mBp);
}




/**
 * If Pm and Mm probes are available on a Chip we can use both to calculate an Expression Measure by performing an Integration
 * on an bivarite distribution.
 * 
 * @param rawIntensityArrayPm The ProbesetIntensityArray with the "raw" (optical corrected) log probe intensities of the Pm probes
 * @param rawIntensityArrayMm The ProbesetIntensityArray with the "raw" (optical corrected) log probe intensities of the Mm probes
 * @param sensitivityProfilePm     SensitivityProfile for the nonspecific binding of Pm probes.
 * @param sensitivityProfileMm     SensitivityProfile for the nonspecific binding of Mm probes.
 * @param distributionParams       A struct holding information of the joind distribution of Pm and Mm noncpecific log intensities (like correlation etc)
 * @param integrationFunction      Function that describes the integral of the bivariate distribution
 * 
 * @return ProbesetIntensitiesArrayPtr Pointer to an ProbesetIntensityArray holding for each probe the desaturated and ns background subtracted log intensities.
 * 
 */
ProbesetIntensitiesArrayPtr  HookcurveAnalyzer::calculateSpecificPortionsBivariate(
		const ProbesetIntensitiesArray&  rawIntensityArrayPm, const ProbesetIntensitiesArray&  rawIntensityArrayMm, 
		const SensitivityProfile& sensitivityProfilePm, const SensitivityProfile& sensitivityProfileMm,
		const BivariateDistributionParameters& distributionParams,
		double (*integrationFunction) (double, void*))
{
  
	BivariateDistributionParameters distribParams = distributionParams;
	distribParams.leftMarginPm     = distribParams.meanPm   -  5 * distribParams.sigmaPm; // Add the margins of the integration.to the other parameters of the struct. 
	distribParams.rightMarginPm  = distribParams.meanPm  +  5 * distribParams.sigmaPm;
	
//	cout << "\n\nDistributionParams: " << distribParams.leftMarginPm << " / " <<  distribParams.meanPm << " / " << distribParams.rightMarginPm << " / " << distribParams.sigmaPm;
//	cout << " / " << distribParams.sigmaAverage << " / " << distribParams.logB << endl;

	
	MathUtil::SimpleBivariateIntegralDistributionParameters bivariateDistribParams;
	bivariateDistribParams.mean = distribParams.meanPm;//   distributionParametersPm.mean;
	bivariateDistribParams.std =   distribParams.sigmaAverage;//         sigmaAverage;
	bivariateDistribParams.correlationCoefficient =   distribParams.correlation;//      correlation;
	
	
	// Integrate denominator for the bivariate integrals. here we actually need to integrate, but only once, since
	// it is the same for all (simplified) bivariate integrations.
	IntensityType inResult, error;
	gsl_function bivariateDistrib;
	bivariateDistrib.function = &MathUtil::simpleBivariateDistributionIntegral;
	gsl_integration_workspace * wDistrib = gsl_integration_workspace_alloc (500);    
	bivariateDistrib.params = &bivariateDistribParams;
	// @bug The bivariate function has a slightly different extend than the simple normal. 
	// In most cases, however, aPM and bPM suffice
	gsl_integration_qag (&bivariateDistrib, distribParams.leftMarginPm /*distribParamsPm.leftMargin*/ ,  distribParams.rightMarginPm /*   distribParamsPm.rightMargin*/, 0.01,  1e-3, 500, GSL_INTEG_GAUSS61,
			wDistrib, &inResult, &error); 
	gsl_integration_workspace_free (wDistrib);
	
	//cout << "Denominator " << inResult << endl;
  
  ProbesetIntensitiesArrayPtr nsSubtractedIntensities =  Probeset::initializeIntensitiesArray(mProbesets);
  
  size_t probesetIndex;
#ifdef _OPENMP
   cout << "Parallel Process. Running " << omp_get_max_threads() << " threads" << endl;
  // Since reduction doesn't work on matrices and valarrays we sort of build our own reduction by
  // assigning vectors in the size of the max thread number and adding the single entries when we leave the
  // parallel region.
#pragma omp parallel for private(probesetIndex)
#endif
  for (probesetIndex = 0; probesetIndex < mProbesets.size(); ++probesetIndex) {

	  (*nsSubtractedIntensities)[probesetIndex] = integrationHelperBivariate(mProbesets[probesetIndex], /*probesetIndex,*/ 
			  rawIntensityArrayPm[probesetIndex], rawIntensityArrayMm[probesetIndex],
			  sensitivityProfilePm, sensitivityProfileMm,
			  distribParams,
			  inResult, &*integrationFunction);   
  }
  cout << endl;
  return nsSubtractedIntensities;
  
}



/**
 * For stabilizing the multicore calculation the outerloop is put into another function
 * @param currentProbeset the probeset currently processed
 * @param  rawIntensitiesPm IntensityArray holding the Pm log values
 * @param  rawIntensitiesMm IntensityArray holding the Pm log values
 * @param sensitivityProfilePm     SensitivityProfile for the nonspecific binding of Pm probes.
 * @param sensitivityProfileMm     SensitivityProfile for the nonspecific binding of Mm probes.
 * @param distribParams       A struct holding information of the joind distribution of Pm and Mm noncpecific log intensities (like correlation etc)
 * @param integrationFunction      Function that describes the integral of the bivariate distribution
 * 
 * @return IntensityArray the IntensityArray holding for each probe the desaturated and ns background subtracted log intensities.
 */
IntensityArray HookcurveAnalyzer::integrationHelperBivariate(const Probeset& currentProbeset,  
										  const IntensityArray& rawIntensitiesPm, const IntensityArray& rawIntensitiesMm,
										  const SensitivityProfile& sensitivityProfilePm, const SensitivityProfile& sensitivityProfileMm,
										  const BivariateDistributionParameters& distribParams,
                                          const IntensityType& denominatorIntegral, 
                                          double (*integrationFunction) (double, void*)) 
{
	assert(rawIntensitiesMm.size() == rawIntensitiesPm.size());
	/*assert(currentProbeset.getSize() == sensitivityCorrectedIntensitiesPm.size());*/
//	assert(sensitivityCorrectedIntensitiesPm.size() == sensitivityCorrectedIntensitiesMm.size());
	
	 // Setup parameters similar for every probe within set
	MathUtil::SimpleBivariateIntegralParameters paramsSimpleDiff;
	paramsSimpleDiff.mean = distribParams.meanPm;//.mean;
	paramsSimpleDiff.std = distribParams.sigmaAverage;//sigmaAverage;
	paramsSimpleDiff.correlationCoefficient = distribParams.correlation;//correlation;
	paramsSimpleDiff.logB = distribParams.logB;//logB;
	
	IntensityArray nsSubtractedIntensitiesArray(0.0, currentProbeset.getSize());
	  // Iterate over the probes of the current probeset perform the integration on each probe
	  for (size_t probeIndex = 0; probeIndex < currentProbeset.getSize(); ++probeIndex) {

	    // Diff parameters.
	    paramsSimpleDiff.desaturatedIntensityPm = mHookModel.getDesaturatedIntensity(rawIntensitiesPm[probeIndex]);
	    paramsSimpleDiff.desaturatedIntensityMm = mHookModel.getDesaturatedIntensity(rawIntensitiesMm[probeIndex]);
	    paramsSimpleDiff.sensitivityPm = sensitivityProfilePm.getSequenceIncrement(currentProbeset.getProbe(probeIndex).getSequence());
	    paramsSimpleDiff.sensitivityMm = sensitivityProfileMm.getSequenceIncrement(currentProbeset.getProbe(probeIndex).getSequenceMm());

	      gsl_function f;
	      f.params = &paramsSimpleDiff;
	      f.function = integrationFunction;

	      IntensityType intResult, intError;
	      gsl_integration_workspace* workspace = gsl_integration_workspace_alloc (500);    
	      gsl_integration_qag (&f, distribParams.leftMarginPm, distribParams.rightMarginPm, 0.01,  1e-3, 500, GSL_INTEG_GAUSS61,
	                           workspace, &intResult, &intError); 
	      gsl_integration_workspace_free (workspace);
	      nsSubtractedIntensitiesArray[probeIndex] = intResult / denominatorIntegral;  
	    }
	  return nsSubtractedIntensitiesArray;
 }




/**
 * 
 * @param rawIntensityArray holding the raw (only optical background corrected) log intensities
 * @param sensitivityCorrectedIntensityArray A reference to the  ProbesetIntensityArray holding the sensitivity corrected probe intensities
 * @param sensitivityProfile Holding the ns sensitivities.
 * @param nsMean Mean value of the Sensitivity Corrected ns intensities.
 * @param nsDeviation The standard deviation of the ns distribution.
 * @param getProbeSequence Function that returns the sequence of a probe (either the PM or MM sequence) default is PM
 * 
 * @return  pointer to the ProbesetIntensityArray holding the ns background subtracted intensities
 * 
 * 
 */
ProbesetIntensitiesArrayPtr  HookcurveAnalyzer::calculateSpecificPortionsMonovariate(const ProbesetIntensitiesArray&  rawIntensityArray, 
//		const ProbesetIntensitiesArray&  sensitivityCorrectedIntensityArray, 
		const SensitivityProfile& sensitivityProfile,
		const DistributionParameters& distributionParameters, 
		const PmMmProbeSequenceFunction& getProbeSequence)
{
  
	DistributionParameters distribParams = distributionParameters; // For security of the multicore precedure we make a copy
	distribParams.leftMargin   = distributionParameters.mean   -   5 * distributionParameters.sigma; // Get the margins of the integration and add it to the other parameters.
	distribParams.rightMargin = distributionParameters.mean   + 5 * distributionParameters.sigma;
  
// cout << distribParams.leftMargin << " / " << distribParams.mean << " / " << distribParams.rightMargin << " / " << distribParams.sigma << endl;
  
  ProbesetIntensitiesArrayPtr nsSubtractedIntensities =  Probeset::initializeIntensitiesArray(mProbesets);
  
  size_t probesetIndex;
#ifdef _OPENMP
  cout << "Parallel Process. Running " << omp_get_max_threads() << " threads" << endl;
  // Since reduction doesn't work on matrices and valarrays we sort of build our own reduction by
  // assigning vectors in the size of the max thread number and adding the single entries when we leave the
  // parallel region.
#pragma omp parallel for private(probesetIndex)
#endif
  for (probesetIndex = 0; probesetIndex < mProbesets.size(); ++probesetIndex) {

	  (*nsSubtractedIntensities)[probesetIndex] = integrationHelperMonovariate(mProbesets[probesetIndex], /*probesetIndex,*/ rawIntensityArray[probesetIndex], /*sensitivityCorrectedIntensityArray[probesetIndex], */
			  distribParams,
			  sensitivityProfile, getProbeSequence);   
  }
  cout << endl;
  return nsSubtractedIntensities;
  
}







/**
 * 
 * Integration helper function for calculating the nsSubtracted Intensities of one IntensityArray
 * (A private function called only from calculateSpecificPortionsMonovariate.)
 * 
 * @param currentProbeset The current probeset 
 * @param rawIntensities IntensityArray holding the raw (only optical background corrected intensities)
 * @param distribParams Holding parameters of the distribution of the nonspecific log intensities, like std and mean
 * @param sensitivityProfile The nonspecific sensitivity Profile
 * @param getProbeSequence Function that gives back the probe sequence of a probe.
 * 
 * @return IntensityArray holding the ns background subtracted intensities.
 * 
 * 
 */

IntensityArray HookcurveAnalyzer::integrationHelperMonovariate(const Probeset& currentProbeset,  
										  const IntensityArray& rawIntensities, 
										  const DistributionParameters& distribParams,
                                          const SensitivityProfile& sensitivityProfile,
                                          const PmMmProbeSequenceFunction& getProbeSequence) 
{
	assert(currentProbeset.getSize() == rawIntensities.size());
	
  // Setup parameters similar for every probe within set
  MathUtil::SimpleIntegralParameters parameters;
  parameters.mean = distribParams.mean;
  parameters.std     = distribParams.sigma;
  
  IntensityArray nsSubtractedIntensitiesArray(0.0, currentProbeset.getSize());
  //assert(nsSubtractedIntensitiesArray.size() == sensitivityCorrectedIntensities.size());
  // Iterate over the probes of the current probeset perform the integration on each probe
  for (size_t probeIndex = 0; probeIndex < currentProbeset.getSize(); ++probeIndex) {    
    // Pm and MM params
	  parameters.sensitivity 
      = sensitivityProfile.getSequenceIncrement(getProbeSequence(currentProbeset.getProbe(probeIndex))/*.getSequence()*/);
    if (1) {
      parameters.desaturatedIntensity = mHookModel.getDesaturatedIntensity(rawIntensities[probeIndex]);
    }
    else { // This is some debug code. Aim: test effect of probe sequence on saturation.
      // for each probe, Imax = M is set to logM = F1* logM+ F2*dN
      IntensityType globalImax = mHookModel.getImax();
      IntensityType F1 = 0.75;
      IntensityType F2 = 2;
      IntensityType probeImax = F1 * globalImax + F2 * parameters.sensitivity;
      mHookModel.setImax(probeImax);
      parameters.desaturatedIntensity = mHookModel.getDesaturatedIntensity(rawIntensities[probeIndex]);
      mHookModel.setImax(globalImax);     
    }
	  
	  //cout << "Increment: " <<  sensitivityProfile.getSequenceIncrement(getProbeSequence(currentProbeset.getProbe(probeIndex))) << endl;
    // cout << "raw/cor: " << rawIntensities[probeIndex] << "\t" << parameters.desaturatedIntensity << endl;
	  
      gsl_function f;
      f.params =  &parameters;  //any_cast<void*>(params[currentParams.gslParameterSet]);
      f.function = &MathUtil::simpleIntegralFunction;

        IntensityType intResult, intError;
      gsl_integration_workspace* workspace = gsl_integration_workspace_alloc (500);    
      gsl_integration_qag (&f, distribParams.leftMargin, distribParams.rightMargin, 0.01,  1e-3, 500, GSL_INTEG_GAUSS61,
                           workspace, &intResult, &intError); 
      gsl_integration_workspace_free (workspace);
      nsSubtractedIntensitiesArray[probeIndex] = intResult;
    }
  
  return nsSubtractedIntensitiesArray;
  } 



/**
 * Computes a list of probeset ids and the respective R value from the theoretic curve.
 *
 * @param intensityType 
 * @param a Parameter a of theoretic function
 * @param f Parameter f of theoretic function
 * @param offset Offset of sumlogi used to fit theoretic function
 */
//std::vector<ProbesetExpressionPair> 
//HookcurveAnalyzer::calculateProbesetR(const IntensityType nsThresholdSumLogI,
//                                      const IntensityType a, const IntensityType f,
//                                      const IntensityType offset) const
//{
//  std::vector<ProbesetExpressionPair> rValues;
//
//  cout << "DEBUG R: nsThresholdSumLogI: " << nsThresholdSumLogI << endl;
//
//  for (size_t probesetIndex = 0; probesetIndex < mProbesets.size(); ++probesetIndex) {
//    // Get probesets latest corrected sumLogI
//    IntensityType sumLogIProbeset = computeAverageIntensity(getSumLogIs(kSensitivityCorrectedPm,
//                                                               kSensitivityCorrectedMm, 
//                                                               probesetIndex));
//
//    // Compute R if sumLogI is in "present fraction", i.e. it is greater than the kink-point
//    IntensityType rProbeset = 0;
//    if (sumLogIProbeset >= mHookModel.getNsThreshold().first) { // nsThresholdSumLogI) {
//      //      sumLogIProbeset = sumLogIProbeset - nsThresholdSumLogI + offset; // transform hookkurve beginning of theo.function
//      sumLogIProbeset = sumLogIProbeset - mHookModel.getNsThreshold().first; // + offset; // transform hookkurve beginning of theo.function
//      rProbeset = ComputeHookcurveModelSumOfSquares::calculateRpmFromSumlogI(a,f, sumLogIProbeset);  
//    }
//    
//    rValues.push_back(ProbesetExpressionPair(// currentProbeset.getProbesetId(),
//                                             boost::lexical_cast<std::string>(sumLogIProbeset)
//                                             //sumLogIProbeset));
//                                                                   , rProbeset));
//  }
//
//  return rValues;
//}



struct evalResultsHelper {
  IntensityPair value;
  size_t index;
};

bool isHelperStructSmaller(const evalResultsHelper& a, const evalResultsHelper& b)
{
  return a.value < b.value;
}
  

/**
 * Calculates moving average over both the the difference and the sum of the log intensities. 
 * For efficiency, the data is written into variable passed as parameter 
 * via reference. All values which a moving average cannot be computed for, are omitted.
 *
 * @param movingAverageWindowSize Size of the moving average window.
 * @param calculationType Which method should be use to calculate the hookcurve.
 *
 * @return The smoothed hookcurve as Sum/Delta pairs.
 * 
 * @note see getDeltaAndSumIntensities
 * @note Maybe use the more generic Math::calculateMovingAverage method
 * 
 */
IntensityMappingPtr HookcurveAnalyzer::getAveragedHookcurvePlot(const ProbesetIntensitiesArray& sumLogIs, const ProbesetIntensitiesArray& deltaLogIs, const int movingAverageWindowSize) const
{     
  // Get the intensities
  vector<IntensityPair> tmpPlot;

  // for all probesets:
  for (size_t probesetIndex = 0; probesetIndex < mProbesets.size(); ++probesetIndex) {
        
    tmpPlot.push_back( IntensityPair(computeAverageIntensity( sumLogIs[probesetIndex] ),
                                     computeAverageIntensity( deltaLogIs[probesetIndex] ) ));
  } // Probeset Iterator
  
  
  // Sort the point on plot with respect to SumLogI
  // @note Sorting works since operator < for std::pair uses lexicographic
  // comparison, i.e. tests if pair1.first < pair2.first first.
  sort(tmpPlot.begin(), tmpPlot.end());

  IntensityMappingPtr intensityPlot =  MathUtil::movingAv2(tmpPlot, movingAverageWindowSize);  

  return intensityPlot;
}



///**
// * Returns a list of values mark_i, i=0..n-1, where for each probe i
// *    mark_i = 1, if SumLogI( Probeset(i) ) > threshold
// *    mark_i = 0, if SumLogI( Probeset(i) ) <= threshold 
// *
// *
// * @param threshold The SumLogI value that separates the marked vs. unmarked.
// *
// */
//std::valarray<double> HookcurveAnalyzer::getSpecificityMarks(const ProbesetIntensitiesArray& sumLogIs, IntensityType threshold) const
//{
//  // Initialize list of marks with 0
//  valarray<double> specificityMarks(0.0, mProbes.size());
//  
//  // Loop over all probesets, remembering the index of the first probe 
//  // of current probeset in mProbes
//  size_t currentFirstProbeIndex = 0;
//  // for(ProbesetVector::const_iterator psi = mProbesets.begin(); psi != mProbesets.end(); ++psi) {
//  for (size_t probesetIndex = 0; probesetIndex < mProbesets.size(); ++probesetIndex) {
//    
//    // Set Mark to 1 if SumLogI > theshold
//    //if (computeAverageIntensity(psi->getSumLogIs()) > threshold) {
////	  IntensityArray sumLog = sumLogIs[probesetIndex];
//    //IntensityArray sumLogIs = getSumLogIs(kPmI, kMmI, probesetIndex); 
//    if (computeAverageIntensity(sumLogIs[probesetIndex]) > threshold) {  
//      slice_array<double> specificityMarksOfCurrentSet = specificityMarks[slice(currentFirstProbeIndex, 
//                                                                                mProbesets[probesetIndex].getSize(), 1)];
//      specificityMarksOfCurrentSet = 1;
//    }
//
//    currentFirstProbeIndex += mProbesets[probesetIndex].getSize(); // psi->getSize();
//  }
//
//  return specificityMarks;
//}


/**
 * Helper function that returns true if p.first < maxValue.
 *
 * @param p A pair.
 * @param maxValue Maximum value.
 */
bool isPairFirstSmaller(const std::pair<IntensityType, IntensityType> p, IntensityType maxValue) 
{
  return p.first < maxValue;
}

/**
 * Calculates the parameters a and f in the theroetic hookcurve model. 
 * This is done by running gradient descend on the least square errors 
 * using the parameters. For the details, see in Preibisch, P.50. 
 *
 * @param hookcurvePlot The SumLogI-DeltaLogI plot obtained from the data.
 *
 * @return Model parameters a and f, and the offset of the fit
 *
 * @todo Pass optimum
 */
boost::tuple<IntensityType, IntensityType, IntensityType> 
HookcurveAnalyzer::estimateTheoreticFunctionParameters(const IntensityMapping& hookcurvePlot) const
{
  // Digitize plot to 10 equally distributed points
  //IntensityMapping digitizedPlot = MathUtil::digitizeCurve(hookcurvePlot, 20);
	IntensityMapping digitizedPlot = MathUtil::digitizeCurveNew(hookcurvePlot);
  digitizedPlot.erase(remove_if(digitizedPlot.begin(), digitizedPlot.end(), 
                                bind2nd(ptr_fun(isPairFirstSmaller), 0)),
                      digitizedPlot.end());
  printVectorPairsToFile(digitizedPlot, "Hookplot_-_digitized.dat");
  //   for_each(digitizedPlot.begin(), digitizedPlot.end(), bind2nd(ptr_fun(isPairFirstSmaller), 0));
  
  // Initialize current minimal values and parameters
  IntensityArray minimumAAndF(2);
  IntensityType minimumError = numeric_limits<IntensityType>::max();

  // Use following initial f and a (starting point for gradient descent)
  IntensityArray initialAAndF(2);
  initialAAndF[0] = pow(10.0,-1);
  initialAAndF[1] = pow(10.0,-8);
  
  IntensityType bestOffset = 0.0;
  IntensityMapping minimumSSQPlot;
  // Try small shifts of the curve for robustness
  for (IntensityType offset = -0.2; offset <= 0.2; offset += 0.02) {
    //IntensityType offset = 0.0;
    // Change the plot values according to offset
    IntensityMapping currentPlot(digitizedPlot.size());
    transform(digitizedPlot.begin(), digitizedPlot.end(), currentPlot.begin(),
              AddToPair<IntensityPair>(IntensityPair(offset, 0.0)));

    // Compute best a and f with gradient descend
    ComputeHookcurveModelSumOfSquares computeLeastSquares(currentPlot);
    IntensityArray currentAAndF = MathUtil::computeGradientDescent(initialAAndF, 
                                                                   computeLeastSquares); 

    // If least square error for these parameters
    // is smaller than currents best, update minimum
    // value, and minumum a and f
    if (computeLeastSquares(currentAAndF) < minimumError) {
      minimumAAndF = currentAAndF;
      bestOffset = offset;
      minimumSSQPlot = currentPlot;
      minimumError = computeLeastSquares(currentAAndF);
    }
  }

  // This plots the energy landscape:
//   ofstream gnuplotFile;
//   ComputeHookcurveModelSumOfSquares computeLeastSquares(minimumSSQPlot);
//   std::valarray<IntensityType> testPoint(2);
//   gnuplotFile.open("landscape.dat");
//   for (IntensityType a = -1; a <= 0; a+= 0.01){// + a/100) {
//     for (IntensityType ff = -5; ff <= -1; ff+= 0.01){// + ff/100) {
//       testPoint[0] = exp10(a); testPoint[1] = exp10(ff);
//       gnuplotFile << a << "\t" << ff << "\t" << (log10(computeLeastSquares(testPoint))) << endl;
//     }
//   }
//   gnuplotFile.close();
  
  return boost::make_tuple(minimumAAndF[0], minimumAAndF[1], bestOffset);
}



/** 
 * Constructor.
 *
 * @param intensityPlot Graph to fit the theretic model to.
 */
ComputeHookcurveModelSumOfSquares::ComputeHookcurveModelSumOfSquares(IntensityMapping& intensityPlot)
  :mIntensityPlot(intensityPlot) {}

/** 
 * Computes the sum of squares error of the hookcurve model using  
 * parameters a and f to the "correct curve".
 *
 * @param aAndF List {a,f}, where a and f are paramters of the theoretic model.
 */
IntensityType ComputeHookcurveModelSumOfSquares::operator() (IntensityArray& aAndF) const 
{    
  IntensityType a = aAndF[0];
  IntensityType f = aAndF[1];
  IntensityType sumOfSquares = 0;

  // For each point of the graph
  for (IntensityMapping::const_iterator currentPoint = mIntensityPlot.begin();
       currentPoint != mIntensityPlot.end(); ++currentPoint) {

    // Given a and f, compute parameter R_Pm from the actual SumLogI
    IntensityType R = calculateRpmFromSumlogI(a, f, currentPoint->first);

    
    // Compute theoretic DeltaLogI of the hookcurve model using R_Pm (and a,f)
    IntensityType deltaLogITheoretic= log10(/*b * */( (R + 1.0) / (a * R + 1.0) )) 
      - log10( (1.0 + f * (R + 1.0)) / (1.0+ f * (a * R + 1.0)) );

    // Add squared deviation from the actual DeltaLogI
    sumOfSquares += pow(deltaLogITheoretic - currentPoint->second, 2);
  }
  return sumOfSquares;
}

/** 
 * Computes the value \f[ R^PM \f] from the actual SumLogI using the model 
 * parameters a and f. See Preibisch, P.52 for precise documentation
 * of the formulas.
 *
 * @param a Model parameter a (curve maximum).
 * @param f Model parameter f (saturation term).
 * @param sumLogI Actual sumLogI of the curve.
 */
IntensityType ComputeHookcurveModelSumOfSquares::calculateRpmFromSumlogI(IntensityType a, IntensityType f, 
                                                                          IntensityType sumLogI) /*const*/
{
  // Use exact formulas as in Priebisch, P.52 but omitting b^-1 (invb)
  double s = (a + 1.0);
  double u = f * f * a;
  double v = f*a + f + f*f*a + f*f;
  double w = f + f*f + f + 1;

  double p = (pow(10,2*sumLogI) * v - s) / (pow(10,2*sumLogI) * u - a);
  double q = (pow(10,2*sumLogI) * w - 1.0) / (pow(10,2*sumLogI) * u - a);

  double R1 = p/-2.0 + sqrt(pow(p,2)/4 -q);
  // double R2 = p/-2.0 - sqrt(pow(p,2)/4 -q); //  p/-2 is always negative -> R2 always negative 
  // We only look for R > 0
  return R1;
}


IntensityType ComputeHookcurveModelSumOfSquares::calculateRpmFromDeltalogI(IntensityType a, IntensityType f, 
    IntensityType deltaLogI) /*const*/
{
  // Use exact formulas as in Priebisch, P.52
  // but omitting b^-1 (invb)
  double invDeltaLogI = pow(10,deltaLogI); 
  
  double mitternachtA = invDeltaLogI*f*a - f*a;
  double mitternachtB = invDeltaLogI*a + invDeltaLogI*f*a + invDeltaLogI*f - 1 - f - f*a;
  double mitternachtC = invDeltaLogI + invDeltaLogI*f - 1 - f;
  
  double R1 = (-mitternachtB + sqrt(pow(mitternachtB,2) - 4 * mitternachtA * mitternachtC)) / 2 * mitternachtA;
  return R1;
}

/** 
 * Constructor.
 *
 * @param intensityPlot1 First graph to fit the theretic model to.
 * @param intensityPlot1 Second graph to fit the theretic model to.

 */
ComputeSumOfTwoHookcurveModelFits::ComputeSumOfTwoHookcurveModelFits(IntensityMapping& intensityPlot1, 
                                                                     IntensityMapping& intensityPlot2)
  :mComputeFit1(intensityPlot1), mComputeFit2(intensityPlot2) {}

/** 
 * Computes the sum of squares error of both hookcurve models by 
 * returning the sum of both fits.
 * 
 * @param params List {a1,f1, a2}, the paramters of the theoretic model.
 */
IntensityType ComputeSumOfTwoHookcurveModelFits::operator() (IntensityArray& params) const
{ 
  // size_t n = params.size() / 2;
  IntensityArray params1 = params[slice(0,2,1)];
  IntensityArray params2(2); //  = params[slice(0,2,1)];
  params2[0] = params[2];
  params2[1] = params[1]; // f2 = f1
  return mComputeFit1(params1) + mComputeFit2(params2);  
  
}

/**
 * Calculates the parameters a and f in the theroetic hookcurve model. 
 * This is done by running gradient descend on the least square errors 
 * using the parameters. For the details, see in Preibisch, P.50. 
 *
 * @param hookcurvePlot The SumLogI-DeltaLogI plot obtained from the data.
 *
 * @return Model parameters a and f, and the offset of the fit
 *
 * @todo Pass optimum
 */
// vector<IntensityType>
boost::tuple<IntensityType, IntensityType, IntensityType, IntensityType> 
HookcurveAnalyzer::estimateTwoTheoreticFunctionParameters(IntensityMappingPtr hookcurvePlot1,
                                                          IntensityMappingPtr hookcurvePlot2) 
{
  vector<IntensityMappingPtr>  plots; // [] = {hookcurvePlot1, hookcurvePlot2};
  plots.push_back(hookcurvePlot1);
  plots.push_back(hookcurvePlot2);

  // Digitize plots to 10 equally distributed points
  vector<IntensityMapping> digitizedPlots(2);
  for (size_t i = 0; i < 2; ++i) {
    digitizedPlots[i] = MathUtil::digitizeCurve(*plots[i], 20);
 
    digitizedPlots[i].erase(remove_if(digitizedPlots[i].begin(), digitizedPlots[i].end(), 
                                     bind2nd(ptr_fun(isPairFirstSmaller), 0)),
                           digitizedPlots[i].end());
    // printVectorPairsToFile(digtizedPlots[i], "Hookplot_-_digitized.dat");
  }

  
  // Initialize current minimal values and parameters
  IntensityArray minimumParams(3);
  IntensityType minimumError = numeric_limits<IntensityType>::max();

  // Use following initial f and a (starting point for gradient descent)
  IntensityArray initialParams(3);
  initialParams[0] = pow(10.0,-1); initialParams[2] = pow(10.0,-1);
  initialParams[1] = pow(10.0,-8); // initialParams[3] = pow(10.0,-8);
  
  IntensityType bestOffset = 0.0;
  // Try small shifts of the curve for robustness
  for (IntensityType offset = -0.2; offset <= 0.2; offset += 0.02) {
    // Change the plot values according to offset
    vector<IntensityMapping> currentPlots(2);
    for (size_t i = 0; i < 2; ++i) {
      currentPlots[i].resize(digitizedPlots[i].size());
      transform(digitizedPlots[i].begin(), digitizedPlots[i].end(), currentPlots[i].begin(),
              AddToPair<IntensityPair>(IntensityPair(offset, 0.0)));
    }
    // Compute best a and f with gradient descend
    ComputeSumOfTwoHookcurveModelFits computeLeastSquares(currentPlots[0], currentPlots[1]);
    IntensityArray currentMinimum = MathUtil::computeGradientDescent(initialParams, 
                                                                     computeLeastSquares); 

    // If least square error for these parameters
    // is smaller than currents best, update minimum
    // value, and minumum a and f
    if (computeLeastSquares(currentMinimum) < minimumError) {
      minimumParams = currentMinimum;
      bestOffset = offset;
      minimumError = computeLeastSquares(currentMinimum);
    }
  }

  // This plots the energy landscape:
  //   ofstream gnuplotFile;
  //   ComputeHookcurveModelSumOfSquares computeLeastSquares(minimumSSQPlot);
  //   std::valarray<IntensityType> testPoint(2);
  //   gnuplotFile.open("landscape.dat");
  //   for (IntensityType a = -1; a <= 0; a+= 0.01){// + a/100) {
  //     for (IntensityType ff = -5; ff <= -1; ff+= 0.01){// + ff/100) {
  //       testPoint[0] = exp10(a); testPoint[1] = exp10(ff);
  //       gnuplotFile << a << "\t" << ff << "\t" << (log10(computeLeastSquares(testPoint))) << endl;
  //     }
  //   }
  //   gnuplotFile.close();
  
  return boost::make_tuple(minimumParams[0], minimumParams[1], minimumParams[2], bestOffset);
}
