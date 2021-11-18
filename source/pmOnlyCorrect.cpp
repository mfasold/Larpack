/**
 * @file pmOnlyCorrect.cpp -
 * @author $Author: jan $
 * @authot Jan Bruecker 
 * @date $Date: 2009-01-28 12:29:54 +0100 (Wed, 28 Jan 2009) $
 * 
 */
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator> 

#include "boost/tuple/tuple.hpp"
#include <boost/lexical_cast.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/function.hpp>

#include "OligoController.hpp"
#include "DataCollector.hpp"
#include "MicroarrayExperiment.hpp"
#include "Probeset.hpp"
#include "HookcurveAnalyzer.hpp"
#include "StringUtil.hpp"
#include "MathUtil.hpp"
#include "BackgroundSubtraction.hpp"
#include "StlUtil.hpp"
#include "ValarrayUtil.hpp"
#include "ProbeFilter.hpp"
#include "ProbesetComposition.hpp"
#include "SensitivityProfile.hpp"
#include "DetectKinkPoint.hpp"

#include "HookcurveStatistics.hpp"
#include "Chip.hpp"
#include "HookModel.hpp"
#include "ExpressionMeasure.hpp"

#include <functional>
#include <numeric>
#include <valarray>


using namespace std;
using namespace larrpack;
using namespace stringutil;
using namespace boost;
using boost::lexical_cast;


size_t gcContant(string seq)
{
	return stringutil::countSubstr("G", seq) + stringutil::countSubstr("C", seq);
}

vector< vector<size_t> > getGcToIndexMap(const vector<Probeset>& probesets)
{
	vector< vector<size_t> > gcMap(26);
	for (size_t i = 0; i < probesets.size(); ++i) {
		if (probesets[i].getProbesetId().substr(0, 15) == "BackgroundProbe") {
			size_t gc = gcContant(probesets[i].getProbe(0).getSequence()); // Is the same for all probes of that probeset
			//cout << "gc " << gc << endl;
			gcMap[gc].push_back(i); // Maps the index of the probeset to the GC contant
		}
	}
	return gcMap;
}


int main (int argc, char *argv[])
{  
	
//	string test = "BackgroundProbe1234455";
//	cout << test.substr(0, 15) << endl;
//	//cout << "kachen an Stelle " << test.find("kachen") << endl;
//	exit(0);
	
  // Read options from File/Commandline
  ProgramOptions options = parseCommandlineOptions(argc, argv); 
  
  // Even if set different, for pmOnly Arrays some options have to be set to fixed levels
  options.chipType = kChiptypeExon; // If we call pmOnlyCorrect we surely want to calculate an pm Only chip, even if we don't say so in the command options.
  options.requiredProfileTypes.clear(); // We only have pm Profiles
  options.requiredProfileTypes.push_back(kEmProfilePm);
  options.expressionMeasureIdentifier = "1"; // We can't calculate other profile types in case of Pm Only arrays.
  
  OligoController controller(options);
  controller.logInitialOptions();
    
  cout << "Importing Experiment and applying background correction" << endl;
  PmMmProbePtrVector probes = controller.importProbeData();
  ChipPtr theChip = controller.createChip(probes);
  
  
  vector< vector<size_t> > gcToIndexMap = getGcToIndexMap(theChip->getProbesets());
//  for (size_t i = 0; i < gcToIndexMap.size(); ++i) {
//	  cout << i << ":   ";
//	  for (size_t  j = 0; j < gcToIndexMap[i].size(); ++j) {
//		  cout << gcToIndexMap[i][j] << "    ";
//	  }
//	  cout << endl;
//  }
  
  HookModelPtr theModel = controller.createHookModel();
  controller.initializeIntensityArraysPmOnly(gcToIndexMap);
  HookcurveAnalyzerPtr hookcurveModel = controller.createHookcurveModel();
  ExpressionMeasurePtr exprCalculater = ExpressionMeasurePtr(new ExpressionMeasure(*theChip, *theModel));
//
  cout << "There are "  << theChip->getProbes().size() 
       << " probes and " << theChip->getProbesets().size() << " probesets."  << endl;
  
  //controller.inverseWashingOfIntensities();
//  
//  // Debug: Use only probesets with certain size
//  if (options.probesetSizeFilter > 0) {
//    cout << "Waning: Debug filter applied.Using only probesets of size " 
//         << options.probesetSizeFilter << endl;
//    theChip->filterProbesets(ProbeCountEquals(options.probesetSizeFilter));
//  }
//
  
  // Create hook curve
  IntensityMappingPtr hookcurvePlot = controller.computeHookcurvePmOnly(kRawIntensitiesOnly, hookcurveModel, gcToIndexMap);
  
  // @Debug:
//  IntensityMappingPtr gcPlot = controller.calculateGCContant(kPmI, kPseudoMmI, hookcurveModel);
//  printVectorPairsToFile(*gcPlot, "gcContant-Primary");
//  printVectorPairsToFile(MathUtil::digitizeCurveNew(*gcPlot), "gcContantDigitized-Primary");
  
  
  // Compute kink point
  IntensityType nsThreshold = options.uncorrectedKinkPoint;
  if (options.computeUncorrectedKinkPoint) {
    nsThreshold = 
      controller.computeKinkCoordinates("Hookplot-Primary.dat", 
                                        MathUtil::digitizeCurveNew(*hookcurvePlot), theModel).first;
  }
  
  
//  IntensityType nsThreshold = 7.0;
  cout << "Uncorrected kink point found at: " << nsThreshold << endl;

  
  // Compute statistics for uncorrected hookcurve
  HookcurveStatistics statsComputer(theModel, hookcurveModel, exprCalculater, options, *theChip, 
                                    controller.getAllProbesetIntensitiesArrays());
//  DataCollector::instance().insert(statsComputer.computeCurveStatistics(nsThreshold, "Uncorrected"));

// boost does not directly support overloaded functions, hence cast is needed for exp10
  UnaryIntensityFunction invLogFunction = (IntensityTypeFunction)exp10;
  UnaryIntensityFunction invGLogFunction = &MathUtil::gExp10;

// Print uncorrected Pms
  printVectorPairsToFile(exprCalculater->calculateProbeExpressionMeasures(*controller.getProbesetIntensityArray(kPmI), kNsPm,
                                                                          mem_fun_ref(&PmMmProbe::getSequence), invGLogFunction), "ExpressionMeasurePmUncorrected_Probe.dat", 
                         "#\t" + options.intensityFilename);


  controller.computeSequenceProfilesPmOnlyChip(kPmI, kPseudoMmI, nsThreshold, string("Primary")/*, false*/,
                                     options.probesetLimitForProfileComputation);

//    cout << "Printing Celfiles" << endl;
//  //    controller.writeIntensitiesToCelfile(options.intensityFilename, "CelRawIntensities.cel", kRawIntensitiesOnly);
//  //    controller.writeIntensitiesToCelfile(options.intensityFilename, "CelNsCorrectedIntensities.cel", kNonspecificProfilesOnly);
//
//  // @debug quits program after computation of profiles
////   if (options.universalParam == "computeProfilesOnly") {
////     cout << "Exiting" << endl;
////     return 0;
////   }
//  
//
//  printVectorPairsToFile(exprCalculater->calculateProbeExpressionMeasures(*controller.getProbesetIntensityArray(kPmI), kNsPm,
//                                                                          mem_fun_ref(&PmMmProbe::getSequence), invGLogFunction), 
//                         "ExpressionMeasurePmNscorrected_Probe.dat", 
//                         "#\t" + options.intensityFilename);
//  
  
  
// Recalculate hookcurve: first correct the provisional intensities.
  controller.updateCorrectedProbesPmOnly(kNonspecificProfilesOnly);
  cout << "\nPlotting new curve NsCorrected..." << endl;
  hookcurvePlot = controller.computeHookcurvePmOnly(kNonspecificProfilesOnly, hookcurveModel, gcToIndexMap);
  
  // @Debug: plots gc contant as function of intensity
//  gcPlot = controller.calculateGCContant(kSensitivityCorrectedPm, kSensitivityCorrectedPseudoMm, hookcurveModel);
//  printVectorPairsToFile(*gcPlot, "gcContant-Corrected");
//  printVectorPairsToFile(MathUtil::digitizeCurveNew(*gcPlot), "gcContantDigitized-Corrected");
  
  
  
  // Compute kink point
  nsThreshold = options.correctedKinkPoint;
  if (options.computeCorrectedKinkPoint) {
    /// @Debug: This line is to save the temporarilly calculated hookcurve.
    //nsThreshold = controller.computeKinkCoordinates("Hookplot-Correctedtmp.dat", MathUtil::digitizeCurveNew(*hookcurvePlot), theModel).first;
    nsThreshold = controller.computeKinkCoordinates("Hookplot-Corrected.dat", MathUtil::digitizeCurveNew(*hookcurvePlot), theModel).first;
  }
  cout << "Corrected kink point found at: " << nsThreshold << endl;
  
  
//  if (options.universalParam == "compute3PrimeBias") {
//    HookcurveStatistics statsComputer(theModel, hookcurveModel, exprCalculater, options, *theChip, controller.getAllProbesetIntensitiesArrays());
//    statsComputer.printProbeDegration("ProbeSumLogs.dat", kSensitivityCorrectedPm, kSensitivityCorrectedMm);
//    statsComputer.print3PrimeBiasAveragedHookcurve("3PrimeBias.dat", kSensitivityCorrectedPm, 
//                                                   kSensitivityCorrectedMm, nsThreshold);
//    cout << "Exiting" << endl;
//    return 0;
//  }


  controller.computeSequenceProfilesPmOnlyChip(kSensitivityCorrectedPm, kSensitivityCorrectedPseudoMm, nsThreshold, string("Corrected")/*,false*/, options.probesetLimitForProfileComputation);


  if (options.universalParam == "computeProfilesOnly") {
    cout << "Exiting" << endl;
    return 0;
  }

  
//  controller.computeSequenceProfiles
//     (nsThreshold, string("Corrected"), options.probesetLimitForProfileComputation);

   if (options.chipType != kChiptypeTiling) {
     controller.computeIntermediateSpecificProfilesPmOnly
       (nsThreshold, string("Corrected"), options.probesetLimitForProfileComputation);
   }
  
  
   controller.calculateNsChipPmOnly();
   
   // Apply the new seqence model
  controller.updateCorrectedProbesPmOnly(kNonspecificProfilesOnly);

  // Calculate Imax
  // @todo: Write  a better method to calculate the theoretical Imax.
  IntensityType iMaxTheor = controller.calculateImaxPmOnlyArray(/**controller.getProbesetIntensityArray(kSensitivityCorrectedPm)*/);
  IntensityType iMax = theModel->estimateImax(*controller.getProbesetIntensityArray(kSensitivityCorrectedPm), 15);
  cout << "Imax theor: " << iMaxTheor << "\tImax av(max 15 probes): " << iMax << endl; 
  theModel->setImax(iMax);
  DataCollector::instance().insert("saturationImax", theModel->getImax()); // Save Imax to the logger 

  
  
  cout << "Recalculate the curve correct with specific profiles." << endl;
  controller.updateCorrectedProbesPmOnly(kAllProfiles);
  hookcurvePlot = controller.computeHookcurvePmOnly(kAllProfiles, hookcurveModel, gcToIndexMap);
  nsThreshold = options.correctedKinkPoint;
  
  if (options.computeCorrectedKinkPoint) {
  	  nsThreshold 
        = controller.computeKinkCoordinates("Hookplot-Corrected.dat", MathUtil::digitizeCurveNew(*hookcurvePlot), theModel).first;
    }
   cout << "Final kink point: " << nsThreshold << endl;
   controller.updateCorrectedProbesPmOnly(kAllProfiles);
  
  
  // Correct the log intensities saved in the Probest instances by subtracting their nonspecific content.
  controller.subtractNBackgroundPmOnly(nsThreshold);
  
  
  // @debug: Plot NS Subtracted Intensities as Hookplot
  ProbesetIntensitiesArray p = *controller.getProbesetIntensityArray(kNsSubstractedPm);
  ProbesetIntensitiesArray m = *controller.createPseudoMmArray(*controller.getProbesetIntensityArray(kNsSubstractedPm), gcToIndexMap);
  IntensityMappingPtr fin = hookcurveModel->getAveragedHookcurvePlot((p+m)/2.0, p-m, options.hookcurveMovingAverageWindowSize);
  printVectorPairsToFile(*fin, "expressionHookPlot.dat"); 

  
  controller.calculateSpecificProfiles(kSensitivityCorrectedPm, kSensitivityCorrectedPseudoMm, kProbeSequenceLength, options.requiredProfileTypes, 
                                       options.probesetLimitForProfileComputation,
                                       nsThreshold);
  
  
  controller.calculateAndPrintExpressionMeasures(exprCalculater);
  
// Calculating and printing statistics:
  DataCollector::instance().insert(statsComputer.computeGeneralStatisticsPmOnly(nsThreshold));
  DataCollector::instance().insert(statsComputer.computeCurveStatisticsPmOnly(nsThreshold, "Corrected"));
  

//  if (options.flags.count("print-probe-statistics") > 0) {
    statsComputer.writeProbeDiagnosisPmOnly("probeStatistics", "\t");
//  }
  statsComputer.writeProbesetDiagnosisPmOnly("probesetStatistics", "\t");
//  statsComputer.print3PrimeBiasAveragedHookcurve("3PrimeBias.dat", kSensitivityCorrectedPm, kSensitivityCorrectedMm, nsThreshold);
//  
//  //   // @debug Write GGG content
//  //   std::vector<ProbesetExpressionPair> gggProbesetVector;
//  //   for (size_t probesetIndex = 0; probesetIndex < theChip->getProbesets().size(); ++probesetIndex) {
//  //     const Probeset& currentProbeset = theChip->getProbesets()[probesetIndex];
//  //     size_t ggg1Count = count_if(currentProbeset.beginProbe(), currentProbeset.endProbe(), 
//  //                                 ProbesequenceHasSubstring("GGG",0));
//  //     gggProbesetVector.push_back(ProbesetExpressionPair(currentProbeset.getProbesetId(), ggg1Count));
//  //   }
//  //   printVectorPairsToFile(gggProbesetVector, "ggg1Frequency.dat");
//
//  cout << endl;
//  
//  // Writing the logger to a file is the last thing we do (maybe we can put this in the destructor of DataCollector ?)
  DataCollector::instance().writeToFile("dataLogger.log");
  return 0;
}


