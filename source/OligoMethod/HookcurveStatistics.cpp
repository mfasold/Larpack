
/**
 * @file HookcurveStatistics.cpp Defines functions to compute statistics on microarray and hookcurve data
 * @author $Author: mario $ 
 * @author Mario Fasold
 * @date $Date: 2008-04-04 17:54:25 +0200 (Fri, 04 Apr 2008) $
 */
#include <iostream>
#include <fstream>
#include <boost/tuple/tuple.hpp>
#include <boost/bind.hpp>

#include <algorithm>
#include <boost/assign/list_of.hpp>
#include <functional>
#include <iostream>
#include <map>
#include <string>

#include "HookcurveStatistics.hpp"
#include "MathUtil.hpp"
#include "ValarrayUtil.hpp"
#include "StringUtil.hpp"
#include "StlUtil.hpp"
#include "DataCollector.hpp"

using namespace std;
using namespace larrpack;
using namespace boost;

/**
 * Constructor
 *
 * @param hookcurveModel The hookcurve under investigation.
 */
HookcurveStatistics::HookcurveStatistics(HookModelPtr hookModel, HookcurveAnalyzerPtr hookcurveModel,
                                         ExpressionMeasurePtr exprMeasure,
                                         const ProgramOptions pOptions, const Chip& chip,
                                         const  std::vector< boost::shared_ptr<ProbesetIntensitiesArray> >& intensitiesArray)
  : mHookModel(hookModel), mHookcurveModel(hookcurveModel), mExpressionMeasureCalculater(exprMeasure), options(pOptions), mChip(chip), mIntensitiesArray(intensitiesArray)
{ }

// Helper function for algorithm
// @note boost lambda operators have portability issues
IntensityType logPlus1(const IntensityType x)
{
  return log(x + 1);
}



/**
 * Calculates the affinity and saturation corrected distribution of nonspecific probes.
 * 
 * @todo: shouldn't this better just return the mean and std ? -> No since we also need to calculate correaltions between Pm and Mm
 * @todo: Do we need all the paramters passed to that function?
 * 
 * @param filter Filters out probesets
 * @param getProbeSequence Function that gives back a probes sequence
 * @param probeType 
 * @param filename If given a filename teh distribution is printed to a file.
 * 
 * @return a vector of all sensitivity corrected probes of one type matching the filter (being nonspecific)
// */
//boost::shared_ptr< vector<IntensityType> > HookcurveStatistics::calculateNonspecificDistribution(const IntensityPredicate filter, 
//                                                    const IntensityMode intensityMode,
//                                                    const string filename) // filname is optional
//{
////   ProbesetIntensitiesArray sumLogIs =  (*mCorrectedIntensities[kSensitivityCorrectedPm] + *mCorrectedIntensities[kSensitivityCorrectedMm]) / 2.0;
////   vector<IntensityType> avgSumLogIs(sumLogIs.size());
////   transform(sumLogIs.begin(), sumLogIs.end(), avgSumLogIs.begin(),
////             computeAverageIntensity);
////   ProbesetIntensitiesArray relevantLogIs 
////     = filterParallel(*mCorrectedIntensities[intensityMode], avgSumLogIs, filter);
//  
//// //   for (size_t i = 0; i < relevantLogIs.size(); ++i) {
//// //     cout << avgSumLogIs[i];
//// //     printValarray(cout, relevantLogIs[i]);
//// //     cout << endl;
//// //   }
// 
////   boost::shared_ptr< vector<IntensityType> > intensitiesVectorPtr(new vector<IntensityType>);
//// //   for_each(relevantLogIs.begin(), relevantLogIs.end(),
//// //            boost::bind(backInsertValarrayElements<IntensityType>, _1, *intensitiesVectorPtr)); 
//
////   for (size_t i = 0; i < relevantLogIs.size(); ++i) {
////     backInsertValarrayElements(relevantLogIs[i], *intensitiesVectorPtr);
////     boost::bind(backInsertValarrayElements<IntensityType>, _1, *intensitiesVectorPtr)(relevantLogIs[i]); // why does that not work
////   }
//  
//  boost::shared_ptr< vector<IntensityType> > intensitiesVectorPtr(new vector<IntensityType>);
//
//  for (size_t probesetIndex = 0; probesetIndex < mProbesets.size(); ++probesetIndex) {
//
//    IntensityArray correctedSumLogI = getSumLogIs(kSensitivityCorrectedPm, kSensitivityCorrectedMm, probesetIndex);
//    if (filter(computeAverageIntensity(correctedSumLogI))) {  
//      IntensityArray setIntensities = getLogIs(intensityMode, probesetIndex);
//      backInsertValarrayElements(setIntensities, *intensitiesVectorPtr);
//    }
//  }
//  return intensitiesVectorPtr;
//}


IntensityTypeStatistics HookcurveStatistics::computeGeneralStatistics(IntensityType nsThreshold, IntensityType a, IntensityType f)
{
 // IntensityType nsChip, nsPmChip, nsMmChip, logB;
  //tie(nsChip, nsPmChip, nsMmChip, logB) = mHookcurveModel->getNsChip();
  size_t probesetCount = mIntensitiesArray[kSensitivityCorrectedPm]->size();
  size_t probeCount    = 0;
  size_t nsProbesCount = 0;
  size_t nsProbesetsCount = 0;
  vector<IntensityType> nsPmProbes;
  vector<IntensityType> nsMmProbes;
  vector<IntensityType> nsPmProbesets;
  vector<IntensityType> nsMmProbesets;
  
  vector<IntensityType> nsPmProbesUncorrected;
  vector<IntensityType> nsMmProbesUncorrected;
  vector<IntensityType> nsPmProbesetsUncorrected;
  vector<IntensityType> nsMmProbesetsUncorrected;
  //   vector<IntensityType> nsSumLogIsCorrected;

  vector<IntensityType> rValuesGreater05;

  for (size_t probesetIndex = 0; probesetIndex < probesetCount; ++probesetIndex) {
    //Probeset& currentProbeset = mHookcurveModel->getProbesets()[probesetIndex];
    //probeCount += currentProbeset.getSize();
    size_t probesetSize = (*mIntensitiesArray[kSensitivityCorrectedPm])[probesetIndex].size();
    probeCount += probesetSize;
    

    
    //IntensityArray sumLogIs = mHookcurveModel->getSumLogIs(kSensitivityCorrectedPm, kSensitivityCorrectedMm, probesetIndex);
    IntensityArray sumLogIs = 0.5 * ( (*mIntensitiesArray[kSensitivityCorrectedPm])[probesetIndex] 
                                                      + (*mIntensitiesArray[kSensitivityCorrectedMm])[probesetIndex] );
    IntensityType sumLogIsAverage = computeAverageIntensity(sumLogIs);
    if (sumLogIsAverage < nsThreshold) {
      valarray<IntensityType> correctedPms = (*mIntensitiesArray[kSensitivityCorrectedPm])[probesetIndex];
      valarray<IntensityType> correctedMms = (*mIntensitiesArray[kSensitivityCorrectedMm])[probesetIndex];
      valarray<IntensityType> uncorrectedPms = (*mIntensitiesArray[kPmI])[probesetIndex];
      valarray<IntensityType> uncorrectedMms = (*mIntensitiesArray[kMmI])[probesetIndex];
      nsPmProbesets.push_back(computeAverageIntensity(correctedPms));
      nsMmProbesets.push_back(computeAverageIntensity(correctedMms));
      nsPmProbesetsUncorrected.push_back(computeAverageIntensity(uncorrectedPms));
      nsMmProbesetsUncorrected.push_back(computeAverageIntensity(uncorrectedMms));
      // nsSumLogIsCorrected.push_back(options.averagingFunction(sumLogIs);
      
      backInsertValarrayElements(correctedPms, nsPmProbes);
      backInsertValarrayElements(correctedMms, nsMmProbes);
      backInsertValarrayElements(uncorrectedPms, nsPmProbesUncorrected);
      backInsertValarrayElements(uncorrectedMms, nsMmProbesUncorrected);

      ++nsProbesetsCount;
      nsProbesCount += probesetSize;
    }    
    else { // sumLogIsAverage >= nsThreshold
      IntensityType rPm = 
        ComputeHookcurveModelSumOfSquares::calculateRpmFromSumlogI(a, f, sumLogIsAverage - nsThreshold);
      if (rPm > 0.5) {
        rValuesGreater05.push_back(rPm);
      }
    }
  }

  IntensityTypeStatistics stats; 
  stats["probesetCount"] = probesetCount;
  stats["probeCount"] = probeCount;

  // Save mean and std statistics
  stats["nsMeanProbesPmCorrected"] = MathUtil::calculateMean(nsPmProbes);
  stats["nsMeanProbesMmCorrected"] = MathUtil::calculateMean(nsMmProbes);
  stats["nsMeanProbessetsPmCorrected"] = MathUtil::calculateMean(nsPmProbesets);
  stats["nsMeanProbessetsMmCorrected"] = MathUtil::calculateMean(nsMmProbesets);
  stats["nsMeanProbessetsDiffCorrected(logB)"] = stats["nsMeanProbessetsPmCorrected"] - stats["nsMeanProbessetsMmCorrected"];
  stats["nsStdProbesPmCorrected"] = MathUtil::calculateStandardDeviation(nsPmProbes);
  stats["nsStdProbesMmCorrected"] = MathUtil::calculateStandardDeviation(nsMmProbes);
  stats["nsStdProbessetsPmCorrected"] = MathUtil::calculateStandardDeviation(nsPmProbesets);
  stats["nsStdProbessetsMmCorrected"] = MathUtil::calculateStandardDeviation(nsMmProbesets);
  stats["nsCorrelationProbesCorrected"] = MathUtil::calculateCorrelation(nsPmProbes, nsMmProbes);
  stats["nsCorrelationProbesetsCorrected"] = MathUtil::calculateCorrelation(nsPmProbesets, nsMmProbesets);
  
  stats["nsMeanProbesPmUncorrected"] = MathUtil::calculateMean(nsPmProbesUncorrected);
  stats["nsMeanProbesMmUncorrected"] = MathUtil::calculateMean(nsMmProbesUncorrected);
  stats["nsMeanProbessetsPmUncorrected"] = MathUtil::calculateMean(nsPmProbesetsUncorrected);
  stats["nsMeanProbessetsMmUncorrected"] = MathUtil::calculateMean(nsMmProbesetsUncorrected);
  stats["nsStdProbesPmUncorrected"] = MathUtil::calculateStandardDeviation(nsPmProbesUncorrected);
  stats["nsStdProbesMmUncorrected"] = MathUtil::calculateStandardDeviation(nsMmProbesUncorrected);
  stats["nsStdProbessetsPmUncorrected"] = MathUtil::calculateStandardDeviation(nsPmProbesetsUncorrected);
  stats["nsStdProbessetsMmUncorrected"] = MathUtil::calculateStandardDeviation(nsMmProbesetsUncorrected);
  stats["nsCorrelationProbesUncorrected"] = 
    MathUtil::calculateCorrelation(nsPmProbesUncorrected, nsMmProbesUncorrected);
  stats["nsCorrelationProbesetsUncorrected"] = 
    MathUtil::calculateCorrelation(nsPmProbesetsUncorrected, nsMmProbesetsUncorrected);


  // Save other statistics
  stats["rGreaterThan05Mean"] = MathUtil::calculateMean(rValuesGreater05);
  vector<IntensityType> rValuesLoggedPlus1(rValuesGreater05.size());
  transform(rValuesGreater05.begin(), rValuesGreater05.end(), rValuesLoggedPlus1.begin(), logPlus1); 
  stats["rGreaterThan05LogPlus1Mean"] = MathUtil::calculateMean(rValuesLoggedPlus1);

  //stats["saturationImax"] = 1 / f;
  //stats["saturationF"] = f;
  //stats["saturationA"] = a;
//   stats["firstPolynom"] = .toString();
//   stats["secondPolynom"] = .toString();

  return stats;
}







IntensityTypeStatistics HookcurveStatistics::computeGeneralStatisticsPmOnly(IntensityType nsThreshold)
{
  size_t probesetCount = mIntensitiesArray[kSensitivityCorrectedPm]->size();
  size_t probeCount    = 0;
  size_t nsProbesCount = 0;
  size_t nsProbesetsCount = 0;
  vector<IntensityType> nsPmProbes;
  vector<IntensityType> nsPmProbesets;
  
  vector<IntensityType> nsPmProbesUncorrected;
  vector<IntensityType> nsPmProbesetsUncorrected;;


  for (size_t probesetIndex = 0; probesetIndex < probesetCount; ++probesetIndex) {
    size_t probesetSize = (*mIntensitiesArray[kSensitivityCorrectedPm])[probesetIndex].size();
    probeCount += probesetSize;
    

    IntensityArray sumLogIs = 0.5 * ( (*mIntensitiesArray[kSensitivityCorrectedPm])[probesetIndex] 
                                                      + (*mIntensitiesArray[kSensitivityCorrectedPseudoMm])[probesetIndex] );
    IntensityType sumLogIsAverage = computeAverageIntensity(sumLogIs);
    if (sumLogIsAverage < nsThreshold) {
      valarray<IntensityType> correctedPms = (*mIntensitiesArray[kSensitivityCorrectedPm])[probesetIndex];
      valarray<IntensityType> uncorrectedPms = (*mIntensitiesArray[kPmI])[probesetIndex];
      nsPmProbesets.push_back(computeAverageIntensity(correctedPms));
      nsPmProbesetsUncorrected.push_back(computeAverageIntensity(uncorrectedPms));
      
      backInsertValarrayElements(correctedPms, nsPmProbes);
      backInsertValarrayElements(uncorrectedPms, nsPmProbesUncorrected);

      ++nsProbesetsCount;
      nsProbesCount += probesetSize;
    }    
  }

  IntensityTypeStatistics stats; 
  stats["probesetCount"] = probesetCount;
  stats["probeCount"] = probeCount;

  // Save mean and std statistics
  stats["nsMeanProbesPmCorrected"] = MathUtil::calculateMean(nsPmProbes);
  stats["nsMeanProbessetsPmCorrected"] = MathUtil::calculateMean(nsPmProbesets);
  stats["nsStdProbesPmCorrected"] = MathUtil::calculateStandardDeviation(nsPmProbes);
  stats["nsStdProbessetsPmCorrected"] = MathUtil::calculateStandardDeviation(nsPmProbesets);
  
  stats["nsMeanProbesPmUncorrected"] = MathUtil::calculateMean(nsPmProbesUncorrected);
  stats["nsMeanProbessetsPmUncorrected"] = MathUtil::calculateMean(nsPmProbesetsUncorrected);
  stats["nsStdProbesPmUncorrected"] = MathUtil::calculateStandardDeviation(nsPmProbesUncorrected);
  stats["nsStdProbessetsPmUncorrected"] = MathUtil::calculateStandardDeviation(nsPmProbesetsUncorrected);

  return stats;
}











/**
 * Computes statistics related to the actual hookcurve and to the current kink-point.
 *
 * @param nsThreshold X-Coordinate of kink point
 * @param suffix String to be appended to the variables in the statig, e.g."Uncorrected"
 */
IntensityTypeStatistics HookcurveStatistics::computeCurveStatistics(IntensityType nsThreshold, std::string suffix) const
{
  IntensityTypeStatistics stats; 

  ProbesetIntensitiesArray sumLogIVector   =  (*(mIntensitiesArray[kSensitivityCorrectedPm]) + *(mIntensitiesArray[kSensitivityCorrectedMm])) * 0.5;
  ProbesetIntensitiesArray deltaLogIVector =   *(mIntensitiesArray[kSensitivityCorrectedPm]) - *(mIntensitiesArray[kSensitivityCorrectedMm]);
  
  // Compute and save hookcurve statistics
  IntensityMappingPtr hookcurvePlot = 
    mHookcurveModel->getAveragedHookcurvePlot(sumLogIVector, deltaLogIVector, options.hookcurveMovingAverageWindowSize);
  valarray<IntensityType> deltaLogIs(hookcurvePlot->size());
  transform(hookcurvePlot->begin(), hookcurvePlot->end(), &deltaLogIs[0], // Populate valarray with deltaLogIs
            boost::bind(&IntensityMapping::value_type::second, _1));        // = SGI select2nd analoga 
  stats["hookcurveWidth" + suffix] = hookcurvePlot->back().first - hookcurvePlot->front().first;
  stats["hookcurveNsRangeWidth" + suffix] = nsThreshold - hookcurvePlot->front().first;
  stats["hookcurveHeight" + suffix] = deltaLogIs.max() - deltaLogIs.min();
  stats["hookcurveKinkPoint" + suffix] = nsThreshold;

  // Compute some counters depending on the kink-point
  // Count uncorrent ns-probes 
  size_t probesetCount = mIntensitiesArray[kSensitivityCorrectedPm]->size();
  cout << "probesetCount: " << probesetCount << endl;
  size_t probeCount    = 0;
  size_t nsProbesCount = 0;
  size_t nsProbesetsCount = 0;
  
  for (size_t probesetIndex = 0; probesetIndex < probesetCount; ++probesetIndex) {
    size_t probesetSize = (*mIntensitiesArray[kSensitivityCorrectedPm])[probesetIndex].size();
    probeCount += probesetSize;
    IntensityArray sumLogIs = 0.5 * ( (*mIntensitiesArray[kSensitivityCorrectedPm])[probesetIndex] 
                                      + (*mIntensitiesArray[kSensitivityCorrectedMm])[probesetIndex] );
    //IntensityArray sumLogIs = mHookcurveModel->getSumLogIs(kPmI, kMmI, probesetIndex);
    IntensityType sumLogIsAverage = computeAverageIntensity(sumLogIs);
    if (sumLogIsAverage < nsThreshold) {
      ++nsProbesetsCount;
      nsProbesCount += probesetSize;
    }
  }
  stats["probesetNsCount" + suffix] = nsProbesetsCount;  
  stats["probesetNsFraction" + suffix] = (IntensityType) nsProbesetsCount / probesetCount;  
  stats["probeNsCount" + suffix] = nsProbesCount;
  stats["probeNsFraction" + suffix] = (IntensityType) nsProbesCount / probeCount;  

  return stats;
}




/**
 * @Todo: This is now more or less the same function as for PM/MM Chips -> Fit this to pmOnly relevant values.
 * 
 * 
 */
IntensityTypeStatistics HookcurveStatistics::computeCurveStatisticsPmOnly(IntensityType nsThreshold, std::string suffix) const
{
  IntensityTypeStatistics stats; 

  ProbesetIntensitiesArray sumLogIVector   =          (*(mIntensitiesArray[kSensitivityCorrectedPm]) + *(mIntensitiesArray[kSensitivityCorrectedPseudoMm])) * 0.5;
  ProbesetIntensitiesArray deltaLogIVector =           *(mIntensitiesArray[kSensitivityCorrectedPm])   - *(mIntensitiesArray[kSensitivityCorrectedPseudoMm]);
  
  // Compute and save hookcurve statistics
  IntensityMappingPtr hookcurvePlot = 
    mHookcurveModel->getAveragedHookcurvePlot(sumLogIVector, deltaLogIVector, options.hookcurveMovingAverageWindowSize);
  valarray<IntensityType> deltaLogIs(hookcurvePlot->size());
  transform(hookcurvePlot->begin(), hookcurvePlot->end(), &deltaLogIs[0], // Populate valarray with deltaLogIs
            boost::bind(&IntensityMapping::value_type::second, _1));        // = SGI select2nd analoga 
  stats["hookcurveWidth" + suffix] = hookcurvePlot->back().first - hookcurvePlot->front().first;
  stats["hookcurveNsRangeWidth" + suffix] = nsThreshold - hookcurvePlot->front().first;
  stats["hookcurveHeight" + suffix] = deltaLogIs.max() - deltaLogIs.min();
  stats["hookcurveKinkPoint" + suffix] = nsThreshold;

  // Compute some counters depending on the kink-point
  // Count uncorrent ns-probes 
  size_t probesetCount = mIntensitiesArray[kSensitivityCorrectedPm]->size();
  size_t probeCount    = 0;
  size_t nsProbesCount = 0;
  size_t nsProbesetsCount = 0;
  for (size_t probesetIndex = 0; probesetIndex < probesetCount; ++probesetIndex) {
    size_t probesetSize = (*mIntensitiesArray[kSensitivityCorrectedPm])[probesetIndex].size();
    probeCount += probesetSize;
    IntensityArray sumLogIs = 0.5 * ( (*mIntensitiesArray[kSensitivityCorrectedPm])[probesetIndex] 
                                      + (*mIntensitiesArray[kSensitivityCorrectedPseudoMm])[probesetIndex] );
    IntensityType sumLogIsAverage = computeAverageIntensity(sumLogIs);
    if (sumLogIsAverage < nsThreshold) {
      ++nsProbesetsCount;
      nsProbesCount += probesetSize;
    }
  }
  stats["probesetNsCount" + suffix] = nsProbesetsCount;  
  stats["probesetNsFraction" + suffix] = (IntensityType) nsProbesetsCount / probesetCount;  
  stats["probeNsCount" + suffix] = nsProbesCount;
  stats["probeNsFraction" + suffix] = (IntensityType) nsProbesCount / probeCount;  
  return stats;
}


/**
 * Writes to a given file a table of probe related statistics
 *
 * @note For genechips probesets are ordered alphabetically by their id when the probesets are defined.
 * 
 */
void HookcurveStatistics::writeProbeDiagnosis(const IntensityType a, const IntensityType f,  const string& fileName, const string& delimiter)
{
  string id = "Probeset_id";
   string probeNr = "Probe_nr_global";
  string probesetCount = "Probeset_count";
  if (options.chipType == kChiptypeTiling) {
    id = "Chromosome";
    probeNr = "location";
    probesetCount = "pseudo_probeset_count";
  }
  
  vector<Probeset> probesets = mChip.getProbesets();
  
  // For gene chips: To assure each probeset always obtains the same number
  // We assign the numbers for the probesets by a lexicographic order of the probeset ids.
  // The numbers for the probes are assigned by the order of occurrence on the chip (x * max(x) + y).
  // For tiling arrays this is not valid.
  map< pair<size_t, size_t>, size_t> coordinateToNumberMapping;
  vector< pair<size_t, size_t> > coordinates;
  if (options.chipType == kChiptypeGenechip) {
    for (size_t probesetIndex = 0; probesetIndex < probesets.size(); ++probesetIndex) {
      for (size_t probeIndex = 0; probeIndex < probesets[probesetIndex].getSize(); ++probeIndex) {
        Probe currentProbe = probesets[probesetIndex].getProbe(probeIndex);
        coordinates.push_back(currentProbe.getPositionPm());
      }
    }
    //   sort(ids.begin(), ids.end());  // For genechips: Probesets are already in this order. we do not need to sort by the ids.
    sort(coordinates.begin(), coordinates.end());
    for (size_t coordinateIdx = 0; coordinateIdx < coordinates.size(); ++coordinateIdx) {
      coordinateToNumberMapping[coordinates[coordinateIdx]] = coordinateIdx;
    }
  } // If
  

  
  string d = delimiter;
  string columns[] = {id, probesetCount, probeNr, "Probe_nr_local", /*"InterrogationNb",*/ "logPM", "logMM", "setAverage_logPM", 
                      "setAverage_logMM", "logPM_corr", "logMM_corr", "setAverage_logPM_corr", "setAverage_logMM_corr", 
                      "yPM", "yMM", "R_PM", "ExpressionMeasure_PM_Exp10_Probe", "ExpressionMeasure_Diff_Exp10_Probe", "ExpressionMeasure_GCRMA-Diff_Exp10_Probe", 
                      "sequence", "center_base", "marker"};
  //valarray<string> columns(columnsA, 17);
  size_t sizeArray = 21;  //  Attention this must be set to the length of the columns array !!!!
  ofstream table;
//   table.precision(4);
  table.open(fileName.c_str());
  table << "#"; // This makes the first line outcommented for gnuplot.
  for (size_t i = 0; i < sizeArray; ++i) {
    table << columns[i] << d;
  }
  table << "\n";
  //size_t probesetCounter = 0;
//  size_t numberProgressPoints = 30;
//  size_t drawPoint = mChip.getProbesets().size() / numberProgressPoints;
  for (size_t probesetIndex = 0; probesetIndex < probesets.size(); ++probesetIndex) {
//    if ((probesetIndex % drawPoint) == 0) { // Draw progress points.
//    }
    string currentProbesetId  = probesets[probesetIndex].getProbesetId();
    
    //  Probeset& currentProbeset = mHookcurveModel->getProbesets()[probesetIndex];
    IntensityType avPm       = computeAverageIntensity((*mIntensitiesArray[kPmI])[probesetIndex]);
    IntensityType avMm       = computeAverageIntensity((*mIntensitiesArray[kMmI])[probesetIndex]);
    IntensityType avPmcorr   = computeAverageIntensity((*mIntensitiesArray[kSensitivityCorrectedPm])[probesetIndex]);
    IntensityType avMmcorr   = computeAverageIntensity((*mIntensitiesArray[kSensitivityCorrectedMm])[probesetIndex]);
    // The sensitivities are calculated here "backwards" from the original and the corrected intensities.
    // this has two advantages 1. It's easier than to calculate them again
    //                         2. We do not need to think about how the sensitivities were calculated (ns only or s)
    IntensityArray yPms       = (*mIntensitiesArray[kPmI])[probesetIndex] - (*mIntensitiesArray[kSensitivityCorrectedPm])[probesetIndex];
    IntensityArray yMms       = (*mIntensitiesArray[kMmI])[probesetIndex] - (*mIntensitiesArray[kSensitivityCorrectedMm])[probesetIndex];
    
    IntensityArray sumLogIs = 0.5 * ( (*mIntensitiesArray[kSensitivityCorrectedPm])[probesetIndex] 
                                                              + (*mIntensitiesArray[kSensitivityCorrectedMm])[probesetIndex] );
    
    valarray<IntensityType> emPmArray, emFastDiffArray, emGcrmaDiffArray;
    if (mExpressionMeasureCalculater->isEmPresent(kEmPmSingleIntegrationExp10))
      emPmArray = mExpressionMeasureCalculater->getEms(kEmPmSingleIntegrationExp10, probesetIndex);
    if (mExpressionMeasureCalculater->isEmPresent(kEmDeltaSingleIntegrationExp10))
      emFastDiffArray = mExpressionMeasureCalculater->getEms(kEmDeltaSingleIntegrationExp10, probesetIndex);
    if (mExpressionMeasureCalculater->isEmPresent(kEmDeltaGcrmaLikeIntegrationExp10))
      emGcrmaDiffArray = mExpressionMeasureCalculater->getEms(kEmDeltaGcrmaLikeIntegrationExp10, probesetIndex);
    
    
    // Maps to each globally assigned probe number a local inner probeset number
    vector< pair<size_t, size_t> > globalToLocalMapperVec;
    if (options.chipType == kChiptypeGenechip) {
      for (size_t probeIndex = 0; probeIndex < probesets[probesetIndex].getSize(); ++probeIndex) {
        PmMmProbe currentProbe = probesets[probesetIndex].getProbe(probeIndex);
        globalToLocalMapperVec.push_back(pair<size_t, size_t>(coordinateToNumberMapping[currentProbe.getPositionPm()], probeIndex));
      }
      sort(globalToLocalMapperVec.begin(), globalToLocalMapperVec.end());
    }
    
    for (size_t i = 0; i < probesets[probesetIndex].getSize(); ++i) {
      size_t probeIndex;
      if (options.chipType == kChiptypeGenechip) {
        probeIndex = globalToLocalMapperVec[i].second;
      }
      else {
        probeIndex = i;
      }
      
      PmMmProbe currentProbe = probesets[probesetIndex].getProbe(probeIndex);
      if (options.chipType == kChiptypeTiling) {
        currentProbesetId = currentProbe.getChromosome();
      }
      
      IntensityType logPm     = (*mIntensitiesArray[kPmI])[probesetIndex][probeIndex];
      IntensityType logMm     = (*mIntensitiesArray[kMmI])[probesetIndex][probeIndex];
      IntensityType logPmcorr = (*mIntensitiesArray[kSensitivityCorrectedPm])[probesetIndex][probeIndex];
      IntensityType logMmcorr = (*mIntensitiesArray[kSensitivityCorrectedMm])[probesetIndex][probeIndex];
      
//       string emPm = "NA";
//       string emDiff = "NA";
//       string emGcrma = "NA";
//       if (mHookcurveModel->isEmPresent(kEmPmSingleIntegrationExp10))
//         emPm = boost::lexical_cast<std::string>(mHookcurveModel->getEms(kEmPmSingleIntegrationExp10, probesetIndex)[probeIndex]);
//       if (mHookcurveModel->isEmPresent(kEmDeltaSingleIntegrationExp10))
//         emDiff = boost::lexical_cast<std::string>(mHookcurveModel->getEms(kEmDeltaSingleIntegrationExp10, probesetIndex)[probeIndex]); 
//       if (mHookcurveModel->isEmPresent(kEmDeltaGcrmaLikeIntegrationExp10))
//         emGcrma = boost::lexical_cast<std::string>(mHookcurveModel->getEms(kEmDeltaGcrmaLikeIntegrationExp10, probesetIndex)[probeIndex]);
      
      IntensityType rPm = ComputeHookcurveModelSumOfSquares::calculateRpmFromSumlogI(a, f, 
                                                                                     sumLogIs[probeIndex] - mHookModel->getNsThreshold().first);
      // Now write everything to the file
      char centerBase = currentProbe.getSequence()[currentProbe.getSequence().size()/2]; // Integer calculation of the center position
      
      // Write a line to the file:
      size_t probeCount = coordinateToNumberMapping[currentProbe.getPositionPm()];
      if (options.chipType == kChiptypeTiling) {
        probeCount = currentProbe.getLocation();
      }
//       table.precision(5);
      table << currentProbesetId << d // The probeset id (For tiling arrays the chromosome)
            << probesetIndex << d      // A number for each probeset (rising in correspondence to the alphabetical order of the probesetids
            << probeCount << d  // A global number for each probe (rising by position on the chip)
            << i     << d // Local number (same order as global number but within a probeset)
            << logPm << d  
            << logMm << d 
            << avPm  << d 
            << avMm  << d 
            << logPmcorr << d // sensitivity corrected values
            << logMmcorr << d 
            << avPmcorr  << d 
            << avMmcorr  << d 
            << yPms[probeIndex] << d // The sequence dependency (sensitivity)
            << yMms[probeIndex] << d 
            << rPm << d;
      if (mExpressionMeasureCalculater->isEmPresent(kEmPmSingleIntegrationExp10_Probe))
        table << mExpressionMeasureCalculater->getEms(kEmPmSingleIntegrationExp10_Probe, probesetIndex)[probeIndex] << d;
      else table << "NA" << d;
      if (mExpressionMeasureCalculater->isEmPresent(kEmDeltaSingleIntegrationExp10_Probe))
        table << mExpressionMeasureCalculater->getEms(kEmDeltaSingleIntegrationExp10_Probe, probesetIndex)[probeIndex] << d; 
      else table << "NA" << d;
      if (mExpressionMeasureCalculater->isEmPresent(kEmDeltaGcrmaLikeIntegrationExp10_Probe))
        table << mExpressionMeasureCalculater->getEms(kEmDeltaGcrmaLikeIntegrationExp10_Probe, probesetIndex)[probeIndex] << d;
      else table << "NA" << d;
      
      table << currentProbe.getSequence() << d
            << centerBase << d 
            << "" << endl;
    }      
  }
}





void HookcurveStatistics::writeProbeDiagnosisPmOnly(const string& fileName, const string& delimiter)
{
  string id = "Probeset_id";
   string probeNr = "Probe_nr_global";
  string probesetCount = "Probeset_count";

  
  vector<Probeset> probesets = mChip.getProbesets();
  
  map< pair<size_t, size_t>, size_t> coordinateToNumberMapping;
  vector< pair<size_t, size_t> > coordinates;
  if (options.chipType == kChiptypeGenechip) {
    for (size_t probesetIndex = 0; probesetIndex < probesets.size(); ++probesetIndex) {
      for (size_t probeIndex = 0; probeIndex < probesets[probesetIndex].getSize(); ++probeIndex) {
        Probe currentProbe = probesets[probesetIndex].getProbe(probeIndex);
        coordinates.push_back(currentProbe.getPositionPm());
      }
    }
    sort(coordinates.begin(), coordinates.end());
    for (size_t coordinateIdx = 0; coordinateIdx < coordinates.size(); ++coordinateIdx) {
      coordinateToNumberMapping[coordinates[coordinateIdx]] = coordinateIdx;
    }
  } // If
  
  
  string d = delimiter;
  string columns[] = {id, probesetCount, probeNr, "Probe_nr_local", /*"InterrogationNb",*/ "logPM", "setAverage_logPM", 
                      "logPM_corr", "setAverage_logPM_corr", "yPM", "ExpressionMeasure_PM_Exp10_Probe" 
                      "sequence", "center_base", "marker"};
  //valarray<string> columns(columnsA, 17);
  size_t sizeArray = 12;  //  Attention this must be set to the length of the columns array !!!!
  ofstream table;
//   table.precision(4);
  table.open(fileName.c_str());
  table << "#"; // This makes the first line outcommented for gnuplot.
  for (size_t i = 0; i < sizeArray; ++i) {
    table << columns[i] << d;
  }
  table << "\n";
  
  for (size_t probesetIndex = 0; probesetIndex < probesets.size(); ++probesetIndex) {
    string currentProbesetId  = probesets[probesetIndex].getProbesetId();
    
    IntensityType avPm       = computeAverageIntensity((*mIntensitiesArray[kPmI])[probesetIndex]);
    IntensityType avPmcorr   = computeAverageIntensity((*mIntensitiesArray[kSensitivityCorrectedPm])[probesetIndex]);
    // The sensitivities are calculated here "backwards" from the original and the corrected intensities.
    // this has two advantages 1. It's easier than to calculate them again
    //                         2. We do not need to think about how the sensitivities were calculated (ns only or s)
    IntensityArray yPms       = (*mIntensitiesArray[kPmI])[probesetIndex] - (*mIntensitiesArray[kSensitivityCorrectedPm])[probesetIndex];
    
    IntensityArray sumLogIs = 0.5 * ( (*mIntensitiesArray[kSensitivityCorrectedPm])[probesetIndex] 
                                                              + (*mIntensitiesArray[kSensitivityCorrectedPseudoMm])[probesetIndex] );

    valarray<IntensityType> emPmArray;
    if (mExpressionMeasureCalculater->isEmPresent(kEmPmSingleIntegrationExp10))
      emPmArray = mExpressionMeasureCalculater->getEms(kEmPmSingleIntegrationExp10, probesetIndex);
  
    // Maps to each globally assigned probe number a local inner probeset number
    vector< pair<size_t, size_t> > globalToLocalMapperVec;
    if (options.chipType == kChiptypeGenechip) {
      for (size_t probeIndex = 0; probeIndex < probesets[probesetIndex].getSize(); ++probeIndex) {
        PmMmProbe currentProbe = probesets[probesetIndex].getProbe(probeIndex);
        globalToLocalMapperVec.push_back(pair<size_t, size_t>(coordinateToNumberMapping[currentProbe.getPositionPm()], probeIndex));
      }
      sort(globalToLocalMapperVec.begin(), globalToLocalMapperVec.end());
    }
    
    for (size_t i = 0; i < probesets[probesetIndex].getSize(); ++i) {
      size_t probeIndex;
      if (options.chipType == kChiptypeGenechip) {
        probeIndex = globalToLocalMapperVec[i].second;
      }
      else {
        probeIndex = i;
      }
      
      PmMmProbe currentProbe = probesets[probesetIndex].getProbe(probeIndex);
      if (options.chipType == kChiptypeTiling) {
        currentProbesetId = currentProbe.getChromosome();
      }
      
      IntensityType logPm     = (*mIntensitiesArray[kPmI])[probesetIndex][probeIndex];
      IntensityType logPmcorr = (*mIntensitiesArray[kSensitivityCorrectedPm])[probesetIndex][probeIndex];
      

      // Now write everything to the file
      char centerBase = currentProbe.getSequence()[currentProbe.getSequence().size()/2]; // Integer calculation of the center position
      
      // Write a line to the file:
      size_t probeCount = coordinateToNumberMapping[currentProbe.getPositionPm()];
      if (options.chipType == kChiptypeTiling) {
        probeCount = currentProbe.getLocation();
      }
//       table.precision(5);
      //table << currentProbesetId << d // The probeset id (For tiling arrays the chromosome)
      table << currentProbe.getProbesetId() << d
            << probesetIndex << d      // A number for each probeset (rising in correspondence to the alphabetical order of the probesetids
            << probeCount << d  // A global number for each probe (rising by position on the chip)
            << i     << d // Local number (same order as global number but within a probeset)
            << logPm << d   
            << avPm  << d  
            << logPmcorr << d // sensitivity corrected values 
            << avPmcorr  << d 
            << yPms[probeIndex] << d; // The sequence dependency (sensitivity)
      if (mExpressionMeasureCalculater->isEmPresent(kEmPmSingleIntegrationExp10_Probe))
        table << mExpressionMeasureCalculater->getEms(kEmPmSingleIntegrationExp10_Probe, probesetIndex)[probeIndex] << d;
      else table << "NA" << d;
      
      table << currentProbe.getSequence() << d 
            << centerBase << d
            << "" << endl;
    }
  }
}





// struct evalResultsHelper {
//   IntensityPair value;
//   size_t index;
// };
bool isPairFirstSmaller(const std::pair<IntensityPair, size_t>& a, const std::pair<IntensityPair, size_t>& b)
{
  return a.first < b.first;
}

IntensityMappingPtr HookcurveStatistics::getAveragedHookcurveInProbesetOrder(IntensityMode pmType, 
                                                                             IntensityMode mmType, 
                                                                             IntensityType defaultIntensity) const
{
  const vector<Probeset>& probesets = mChip.getProbesets();

  // Get the intensities
  vector<IntensityPair> tmpPlot;

  // Get (corrected) intensities for all probesets:
  for (size_t probesetIndex = 0; probesetIndex < probesets.size(); ++probesetIndex) { 
	  
	  IntensityArray correctedSumLogI = 0.5 * ( (*mIntensitiesArray[pmType])[probesetIndex] 
	                                                                + (*mIntensitiesArray[mmType])[probesetIndex] );
//	  
//    IntensityArray correctedSumLogI = 
//      mHookcurveModel->getSumLogIs(pmType, mmType, probesetIndex);
//    IntensityArray correctedDeltaLogI = 
//      mHookcurveModel->getDeltaLogIs(pmType, mmType, probesetIndex);

	  IntensityArray correctedDeltaLogI = (*mIntensitiesArray[pmType])[probesetIndex] 
	       	                                                  - (*mIntensitiesArray[mmType])[probesetIndex];
	  
    tmpPlot.push_back( IntensityPair(computeAverageIntensity(correctedSumLogI),
                                     computeAverageIntensity(correctedDeltaLogI) ));
  } 
  
  // Fill a helper vector with plot pairs and indices...
  vector<pair<IntensityPair, size_t> > mappingHelper;
  for(size_t i = 0; i < tmpPlot.size(); ++i) {
    mappingHelper.push_back(make_pair(tmpPlot[i], i));
  }  

  // and sort it the same way as the plot, to have sorted indices
  sort(mappingHelper.begin(), mappingHelper.end(), isPairFirstSmaller);

  // Create a vector in original probeset order containing the new index after sorting
  vector<pair<IntensityPair, size_t> > evaluateResultsHelper(mappingHelper.size());
  for(size_t i = 0; i < mappingHelper.size(); ++i) {
    evaluateResultsHelper[mappingHelper[i].second] =  make_pair(mappingHelper[i].first, i);
  }  

  // Sort the point on plot with respect to SumLogI
  // @note Sorting works since operator < for std::pair uses lexicographic
  // comparison, i.e. tests if pair1.first < pair2.first first.
  sort(tmpPlot.begin(), tmpPlot.end());

  // Compute the moving averaged plot
  IntensityMappingPtr intensityPlot =
	  MathUtil::movingAv2(tmpPlot, options.hookcurveMovingAverageWindowSize);
    //MathUtil::calculateMovingAverage(tmpPlot, options.hookcurveMovingAverageWindowSize);  

  // Calculate first Point which the moving average can be calculated for
  size_t firstAveragedDatapoint = options.hookcurveMovingAverageWindowSize/2 + 1 - 1; // integer division
  
  // Change order of averaged probesets
  IntensityMappingPtr probesetHookvalues = IntensityMappingPtr(new IntensityMapping(probesets.size()));
  for (size_t probesetIndex = 0; probesetIndex < probesets.size(); ++probesetIndex) {
    // If probeset index within the range of indexes a average could be computed for
    if (firstAveragedDatapoint <= evaluateResultsHelper[probesetIndex].second &&
        evaluateResultsHelper[probesetIndex].second < intensityPlot->size() + firstAveragedDatapoint) {
      (*probesetHookvalues)[probesetIndex] = 
        (*intensityPlot)[evaluateResultsHelper[probesetIndex].second - firstAveragedDatapoint];
    }
    // otherwise set it to an default value
    else {
      (*probesetHookvalues)[probesetIndex] = make_pair(defaultIntensity, defaultIntensity);
    }
  }

  return probesetHookvalues;
}

/**
 * Prints the difference to sumLogI for each probe
 *
 *
 */
void HookcurveStatistics::printProbeDegration(std::string filename,
                                              IntensityMode pmType, 
                                              IntensityMode mmType) const
{
  ofstream hookplotFile;
  hookplotFile.open(filename.c_str());

  const vector<Probeset>& probesets = mChip.getProbesets();
  for (size_t probesetIndex = 0; probesetIndex < probesets.size(); ++probesetIndex) {
    IntensityArray pmLogI = (*mIntensitiesArray[pmType])[probesetIndex];
    IntensityType avgPmLogI = computeAverageIntensity(pmLogI);


    for (size_t probeIndex = 0; probeIndex < probesets[probesetIndex].getSize(); ++probeIndex) {
      const Probe& p = probesets[probesetIndex].getProbe(probeIndex);

      hookplotFile << p.getProbesetId() << "\t" << p.getLocation() << "\t" 
                   << p.getSequence() << "\t" 
                   << pmLogI[probeIndex] << "\t" 
                   << avgPmLogI << "\t" 
                   << endl;
    }
  }
    
  hookplotFile.close();
}


/**
 * Prints a special SumLogI plot to a file (3' and 5' bias)
 *
 * @todo Clean this debug function
 */
void HookcurveStatistics::print3PrimeBiasAveragedHookcurve(std::string filename,
                                                           IntensityMode pmType, 
                                                           IntensityMode mmType,
                                                           IntensityType nsThreshold) const
{
  const size_t biasProbeCount = 3;
  const size_t movingAverageSize = 600;
  const size_t minimumMovingAverageSize = 100; // used 20 for extreme values


  const vector<Probeset>& probesets = mChip.getProbesets();

  ProbesetIntensitiesArray correctedSumLogIs = (*mIntensitiesArray[pmType] + *mIntensitiesArray[mmType]) / 2.0;

  // Get (corrected) sumLog intensities of 3' and 5' probes for all probesets:
  typedef map<string, vector<IntensityPair> > MapOfPlots;
  MapOfPlots plots;
  for (size_t probesetIndex = 0; probesetIndex < probesets.size(); ++probesetIndex) { 
    // Sort out relevant probe intensities
    IntensityArray sumLogI = correctedSumLogIs[probesetIndex];
    //IntensityType sumLogI = correctedSumLogIs[probesetIndex];
    IntensityArray pmI = (*mIntensitiesArray[pmType])[probesetIndex];
  
	  // IntensityArray correctedSumLogI = 0.5 * ( (*mIntensitiesArray[pmType])[probesetIndex] 
	  //                                                               + (*mIntensitiesArray[mmType])[probesetIndex] );
    // Compute the average of the first and last 4 (=biasProbeCount) elements
    size_t n = probesets[probesetIndex].getSize();
    // if (n==11) { // @debug: If size is 11
      plots["lowerZange"].push_back
        (IntensityPair(computeAverageIntensity(sumLogI), 
                       computeAverageIntensity(pmI[slice(0,biasProbeCount,1)])
                       - computeAverageIntensity(pmI)));
      plots["upperZange"].push_back
        (IntensityPair(computeAverageIntensity(sumLogI), 
                       computeAverageIntensity(pmI[slice(n-biasProbeCount,biasProbeCount,1)]) 
                       - computeAverageIntensity(pmI)));
      plots["middleZange"].push_back
        (IntensityPair(computeAverageIntensity(sumLogI), 
                       computeAverageIntensity(pmI[slice(n/2 -2,biasProbeCount,1)]) 
                       - computeAverageIntensity(pmI)));
      // }
  }
  // Sort the points on the plot with respect to SumLogI
  // @note Sorting works since operator < for std::pair uses lexicographic
  // comparison, i.e. tests if pair1.first < pair2.first first.
  for (MapOfPlots::iterator plot = plots.begin(); plot != plots.end(); ++plot) {
    sort(plot->second.begin(), plot->second.end());
  }
  
  // Compute the moving averaged plot
  typedef map<string, IntensityMappingPtr > MapOfAveragedPlots;
  MapOfAveragedPlots averagedPlots;
  for (MapOfPlots::const_iterator plot = plots.begin(); plot != plots.end(); ++plot) {
    // averagedPlots[plot->first] = MathUtil::calculateMovingAverage(plot->second, 
    //                                                               // options.hookcurveMovingAverageWindowSize);
    //                                                               //3 * options.hookcurveMovingAverageWindowSize);
    //                                                               500);
    // @todo Do not use constants here.
    averagedPlots[plot->first] = MathUtil::movingAv2(plot->second, movingAverageSize, minimumMovingAverageSize);
  }

  // Print "Zange" to file
  ofstream hookplotFile;
  hookplotFile.open(filename.c_str());

  hookplotFile << "# SumLogI";
  for (MapOfPlots::const_iterator plot = plots.begin(); plot != plots.end(); ++plot) {
    hookplotFile << "\t" << plot->first;
  }
  // DeltaLogI\tSumLogI(Small Location)\tSumLogI(Big Location)" << endl;
  hookplotFile << endl;

  const IntensityMappingPtr firstPlot = averagedPlots.begin()->second;
  for (size_t plotIndex = 0; plotIndex < firstPlot->size(); ++plotIndex) {
      hookplotFile << (*firstPlot)[plotIndex].first;
      for (MapOfAveragedPlots::const_iterator plot = averagedPlots.begin(); plot != averagedPlots.end(); ++plot) {
        hookplotFile << "\t" << (*plot->second)[plotIndex].second;
    }
    hookplotFile << endl;
  }
  hookplotFile.close();

  // Write a degradation measure: robust tongs opening
  averagedPlots["degradationHook"] = MathUtil::movingAv2(plots["upperZange"], movingAverageSize, minimumMovingAverageSize); // intialize
  for (size_t i = 0; i < averagedPlots["degradationHook"]->size(); ++i) {
    (*averagedPlots["degradationHook"])[i].second = 
      (*averagedPlots["upperZange"])[i].second - (*averagedPlots["lowerZange"])[i].second;
  }
  
  IntensityMapping digitizedZange = MathUtil::digitizeCurve(*averagedPlots["degradationHook"], 8);
  IntensityType bias = max_element(digitizedZange.begin(), digitizedZange.end(), 
                                   boost::bind(&IntensityMapping::value_type::second,_1) < boost::bind(&IntensityMapping::value_type::second,_2))->second;
  DataCollector::instance().insert("DegradationTongsOpening", bias);

  // Write a secind distance measure as SumLogI(3'Probes)/SumLogI(all Probes)
  IntensityType p3Sum = 0;
  IntensityType p5Sum = 0;
  size_t nsProbeCount = 0;
  for (size_t probesetIndex = 0; probesetIndex < probesets.size(); ++probesetIndex) { 
    IntensityType sumLogI = computeAverageIntensity(correctedSumLogIs[probesetIndex]);
    if (sumLogI > nsThreshold + 0.3) {
      // Compute the average of the first and last 3 (=biasProbeCount) elements
      size_t n = probesets[probesetIndex].getSize();
      
      p5Sum += (*mIntensitiesArray[pmType])[probesetIndex][n-1];
      p3Sum += (*mIntensitiesArray[pmType])[probesetIndex][0];    
      ++nsProbeCount;
    }
  }

  DataCollector::instance().insert("DegradationDk", exp10(p3Sum/nsProbeCount) / exp10(p5Sum/nsProbeCount));

//   IntensityMapping digitizedZange = MathUtil::digitizeCurve(*averagedPlots["upperZange"], 8);
//   IntensityType bias = max_element(digitizedZange.begin(), digitizedZange.end(), 
//                                    boost::bind(&IntensityMapping::value_type::second,_1) < boost::bind(&IntensityMapping::value_type::second,_2))->second;
//   DataCollector::instance().insert("ReverseTranscriptionBias", bias);
  
}

// Debug function to compute the gc-content of a probe
double getProbesetGcContent(const Probeset &ps) {
  size_t n = ps.getSize();
  size_t gccount = 0;
  for(size_t i = 0; i < n; ++i) {
    string seq = ps.getProbe(i).getSequence();
    for(size_t j = 0; j < seq.size(); ++j) {
      if (seq[j] == 'G' || seq[j] == 'C') 
        ++gccount;
    }    
    //     cout << seq  << "\t" << gccount << endl;
  }
  return (double)gccount/(n*25.0);
}

/**
 * Writes to a given file a table of probe related statistics
 * 
 * 
 */
void HookcurveStatistics::writeProbesetDiagnosis(const IntensityType a, const IntensityType f, const string& fileName, const string& delimiter)
{
  string d = delimiter;
  string na = "NA";
  // Those fields that are always valid:
  string columnsA[] = {"Probeset_id", "Probeset_nr", "setAverage_logPM", "setAverage_logMM", 
                       "setAverage_logPM_corr", "setAverage_logMM_corr", "xHook_raw", "yHook_raw", "xHook_corr", "yHook_corr", "R_PM", 
                       "ExpressionMeasure_PM_Exp10", "ExpressionMeasure_Diff_Exp10", "ExpressionMeasure_GCRMA-Diff_Exp10", "marker", "presence"};
  size_t arraySize = 16;
  ofstream table;
  table.open(fileName.c_str());
//   table.precision(5);
  table << "#"; // This makes the first line outcommented for gnuplot.
  for (size_t i = 0; i < arraySize; ++i) { // Write the descriptor line
    table << columnsA[i] << d;
  }
  table << "\n";
  const vector<Probeset>& probesets = mChip.getProbesets();


  IntensityMappingPtr uncorrectedHookcurve = 
    getAveragedHookcurveInProbesetOrder(kPmI, kMmI);

  IntensityMappingPtr correctedHookcurve = 
    getAveragedHookcurveInProbesetOrder(kSensitivityCorrectedPm, kSensitivityCorrectedMm);

  size_t probesetCounter = 0;
  for (size_t probesetIndex = 0; probesetIndex < probesets.size(); ++probesetIndex) {
    IntensityArray sumLogIs = 0.5 * ( (*mIntensitiesArray[kSensitivityCorrectedPm])[probesetIndex] 
                                      + (*mIntensitiesArray[kSensitivityCorrectedMm])[probesetIndex] );
//    IntensityArray sumLogIs   = mHookcurveModel->getSumLogIs(kSensitivityCorrectedPm, kSensitivityCorrectedMm, probesetIndex);
    IntensityType avSumLogI = computeAverageIntensity(sumLogIs);
    IntensityType rPm = ComputeHookcurveModelSumOfSquares::calculateRpmFromSumlogI(a, f, 
                                                                                   avSumLogI - mHookModel->getNsThreshold().first);
    if (avSumLogI < mHookModel->getNsThreshold().first) {
      rPm = 0;
    }
    
    string currentProbesetId    = probesets[probesetIndex].getProbesetId();
    IntensityType avPm         = computeAverageIntensity((*mIntensitiesArray[kPmI])[probesetIndex]);
    IntensityType avMm         = computeAverageIntensity((*mIntensitiesArray[kMmI])[probesetIndex]);
    IntensityType avPmcorr     = computeAverageIntensity((*mIntensitiesArray[kSensitivityCorrectedPm])[probesetIndex]);
    IntensityType avMmcorr     = computeAverageIntensity((*mIntensitiesArray[kSensitivityCorrectedMm])[probesetIndex]);
       
    /// @note If no hookcurve mapping is defined (which happens to the smoothingWindowsize/2 smallest and biggest probesets) they are
    /// mapped onto 0/0  (if not defined differently) 
    //     table.precision(5);
    table << currentProbesetId << d << probesetCounter << d
          << avPm << d << avMm << d << avPmcorr << d << avMmcorr << d  
          << (*uncorrectedHookcurve)[probesetIndex] << d << (*correctedHookcurve)[probesetIndex] << d // Position of the probeset in the hookcurve
          << rPm << d;
    if (mExpressionMeasureCalculater->isEmPresent(kEmPmSingleIntegrationExp10))
      table << computeAverageIntensity(mExpressionMeasureCalculater->getEms(kEmPmSingleIntegrationExp10, probesetIndex)) << d;
    else table << "NA" << d;
    if (mExpressionMeasureCalculater->isEmPresent(kEmDeltaSingleIntegrationExp10))
      table << computeAverageIntensity(mExpressionMeasureCalculater->getEms(kEmDeltaSingleIntegrationExp10, probesetIndex)) << d;
    else table << "NA" << d;
    if (mExpressionMeasureCalculater->isEmPresent(kEmDeltaGcrmaLikeIntegrationExp10))
      table << computeAverageIntensity(mExpressionMeasureCalculater->getEms(kEmDeltaGcrmaLikeIntegrationExp10, probesetIndex)) << d;
    else table << "NA" << d;
    table << getProbesetGcContent(probesets[probesetIndex]) << d;
    if (avSumLogI < mHookModel->getNsThreshold().first) {
      table << "absent";
    } else {
      table << "present";
    }
    table << endl;

    ++probesetCounter;      
  }
  table.close();
}

/**
 * Writes statistics to a file that contains the intra-probeset variance measures
 *
 *
 */
void HookcurveStatistics::writeProbesetVariance(const string fileName, const string d) const
{
  ofstream table;
  table.open(fileName.c_str());
  table << "#Probeset_id" << d << "SD_PM" << d << "MAD_PM" << d << "SD_MM" << d << "MAD_MM" << d << "SD_EM_PmExp10" << d << "MAD_EM_PmExp10" << endl;

  const vector<Probeset>& probesets = mChip.getProbesets();
  for (size_t probesetIndex = 0; probesetIndex < probesets.size(); ++probesetIndex) {
    table << probesets[probesetIndex].getProbesetId() << d 
          << MathUtil::calculateStandardDeviation((*mIntensitiesArray[kSensitivityCorrectedPm])[probesetIndex]) << d
          << MathUtil::calculateStandardDeviation((*mIntensitiesArray[kSensitivityCorrectedMm])[probesetIndex]) << d
          << MathUtil::calculateMAD((*mIntensitiesArray[kSensitivityCorrectedPm])[probesetIndex]) << d
          << MathUtil::calculateMAD((*mIntensitiesArray[kSensitivityCorrectedMm])[probesetIndex]) << d;
    if (mExpressionMeasureCalculater->isEmPresent(kEmPmSingleIntegrationExp10_Probe)) {      
      table << MathUtil::calculateStandardDeviation(log10(mExpressionMeasureCalculater->getEms(kEmPmSingleIntegrationExp10_Probe, probesetIndex))) << d
            << MathUtil::calculateMAD(log10(mExpressionMeasureCalculater->getEms(kEmPmSingleIntegrationExp10_Probe, probesetIndex)));    
    }
    table << endl;
  }
  table.close();
}



void HookcurveStatistics::writeProbesetDiagnosisPmOnly(const string& fileName, const string& delimiter)
{
  string d = delimiter;
  string na = "NA";
  // Those fields that are always valid:
  // string columnsA[] = {"Probeset_id", "Probeset_nr", "setAverage_logPM",  
  //                      "setAverage_logPM_corr", "xHook_raw", "yHook_raw", "xHook_corr", "yHook_corr", 
  //                      "ExpressionMeasure_PM_Exp10",  "marker"};
  // size_t arraySize = 10;

  // Note: I also print many empty/dummy values - this file should be compatible with PmMm Variant!
  string columnsA[] = {"Probeset_id", "Probeset_nr", "setAverage_logPM", "setAverage_logMM", 
                       "setAverage_logPM_corr", "setAverage_logMM_corr", 
                       "xHook_raw", "yHook_raw", "xHook_corr", "yHook_corr", "R_PM", 
                       "ExpressionMeasure_PM_Exp10", "ExpressionMeasure_Diff_Exp10", "ExpressionMeasure_GCRMA-Diff_Exp10", 
                       "marker", "presence"};
  size_t arraySize = 16;


  ofstream table;
  table.open(fileName.c_str());
//   table.precision(5);
  table << "#"; // This makes the first line outcommented for gnuplot.
  for (size_t i = 0; i < arraySize; ++i) { // Write the descriptor line
    table << columnsA[i] << d;
  }
  table << "\n";
  const vector<Probeset>& probesets = mChip.getProbesets();


  IntensityMappingPtr uncorrectedHookcurve = 
    getAveragedHookcurveInProbesetOrder(kPmI, kPseudoMmI);

  IntensityMappingPtr correctedHookcurve = 
    getAveragedHookcurveInProbesetOrder(kSensitivityCorrectedPm, kSensitivityCorrectedPseudoMm);

  size_t probesetCounter = 0;
  for (size_t probesetIndex = 0; probesetIndex < probesets.size(); ++probesetIndex) {    
    string currentProbesetId    = probesets[probesetIndex].getProbesetId();
    IntensityType avPm         = computeAverageIntensity((*mIntensitiesArray[kPmI])[probesetIndex]);
    IntensityType avPmcorr     = computeAverageIntensity((*mIntensitiesArray[kSensitivityCorrectedPm])[probesetIndex]);
    IntensityType avPseudoMm         = computeAverageIntensity((*mIntensitiesArray[kPseudoMmI])[probesetIndex]);
    IntensityType avPseudoMmcor     = computeAverageIntensity((*mIntensitiesArray[kSensitivityCorrectedPseudoMm])[probesetIndex]);
    
    // Take care that these have been calculated and printed to the respective array.

    
    /// @note If no hookcurve mapping is defined (which happens to the smoothingWindowsize/2 smallest and biggest probesets) they are
    /// mapped onto 0/0  (if not defined differently) 
//     table.precision(5);
    table << currentProbesetId << d << probesetCounter << d
          << avPm  << d << avPseudoMm<< d << avPmcorr << d << avPseudoMmcor << d
          << (*uncorrectedHookcurve)[probesetIndex] << d << (*correctedHookcurve)[probesetIndex] << d; // Position of the probeset in the hookcurve
    table << na << d; // R_PM not yet available for PM-only
    if (mExpressionMeasureCalculater->isEmPresent(kEmPmSingleIntegrationExp10))
      table << computeAverageIntensity(mExpressionMeasureCalculater->getEms(kEmPmSingleIntegrationExp10, probesetIndex)) << d;
    else table << "NA" << d;
    // @note print dummy values
    table << na << d; // ExpressionMeasure_Diff_Exp10
    table << na << d; // ExpressionMeasure_GCRMA-Diff_Exp10   
    table << getProbesetGcContent(probesets[probesetIndex]) << d;

    // @todo Check if corrent, only inserted temporarily
    if (avPmcorr < mHookModel->getNsThreshold().first) { 
        table << "absent";
      } else {
        table << "present";
      }
    table << endl;
    
    ++probesetCounter;      
  }
}



// void HookcurveStatistics::computeGGGStatistics(const IntensityType nsThreshold) {
//   string motifs[] = {"GGG", "CCC", "AAA", "TTT", 
//                      "GGG1", "CCC1", "AAA1", "TTT1"
//                      "GGG12", "CCC12", "AAA12", "TTT12"}
//   // mising: GG
//   const size_t motifCount = 3*4;
  
//   vector<size_t> motifCounter(motifCount, 0);
//   vector<IntensityType> motifLogISum(motifCount, 0.0);

//   vector<Probeset> probesets = mChip.getProbesets();
//   ProbesetIntensitiesArray correctedSumLogIs = (*mIntensitiesArray[pmType] + *mIntensitiesArray[mmType]) / 2.0;

//   for (size_t probesetIndex = 0; probesetIndex < probesets.size(); ++probesetIndex) {    
//     IntensityType sumLogI = computeAverageIntensity(correctedSumLogIs[probesetIndex]);
//     if (sumLogI > nsThreshold + 0.3) {

//       for (size_t probeIndex = 0; probeIndex < probesets[probesetIndex].getSize(); ++probeIndex) {
//         string seqeuence = probesets[probesetIndex].getProbe(probeIndex).getSequence();
//         IntensityType logI = (*mIntensitiesArray[pmType])[probesetIndex][probeIndex];

//         size_t startIndex = 0;
//         for (size_t i = 0*4; i < 1*4; ++i) {
//           if (sequence.find(motifs[i]) != std::string::npos) {
//             motifCounter[i] = motifCounter[i] + 1;
//             motifLogISum[i] = motifLogISum[i] + logI;
//           }
//         }

//         startIndex = 1*4;
//         for (size_t i = 0*4; i < 1*4; ++i) {
//           if (sequence.find(motifs[i]) == 1) {
//             motifCounter[i + startIndex] = motifCounter[i + startIndex] + 1;
//             motifLogISum[i + startIndex] = motifLogISum[i + startIndex] + logI;
//           }
//         }

//         startIndex = 2*4;
//         for (size_t i = 0*4; i < 1*4; ++i) {
//           if (sequence.find(motifs[i]) == 12) {
//             motifCounter[i + startIndex] = motifCounter[i + startIndex] + 1;
//             motifLogISum[i + startIndex] = motifLogISum[i + startIndex] + logI;
//           }
//         }
//       }
        
//     }
//   }

//   ofstream oFile;
//   string targetFilename = "motifStats2.dat";
//   oFile.open(targetFilename.c_str());
//   for (size_t i = 0; i < 3*4 +3; ++i) {
//     oFile << motifs[i] << "\t" << motifCounter[i] << "\t" << motifLogISum[i] << "\t" << (motifLogISum[i]/(IntensityType)motifCounter[i]) << endl;
//   }
//   oFile.close();

//   // x1 = GGG1 / (AAA1 + TTT1) * 2
//   // x1 = CCC12 / (AAA12 + TTT12) * 2
//   // GGG = GGG_1 / GGG


// }
