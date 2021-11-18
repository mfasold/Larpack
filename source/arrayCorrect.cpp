
/**
 * @file arrayCorrect.cpp -
 * @author $Author$ 
 * @date $Date$
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


int main (int argc, char *argv[])
{  
  // Read options from File/Commandline
  ProgramOptions options = parseCommandlineOptions(argc, argv); 
  OligoController controller(options);
  controller.logInitialOptions();
    
  cout << "Importing Experiment and applying background correction" << endl;
  PmMmProbePtrVector probes = controller.importProbeData();
  ChipPtr theChip = controller.createChip(probes);
  HookModelPtr theModel = controller.createHookModel();
  controller.initializeIntensityArrays();
  HookcurveAnalyzerPtr hookcurveModel = controller.createHookcurveModel();
  ExpressionMeasurePtr exprCalculater = ExpressionMeasurePtr(new ExpressionMeasure(*theChip, *theModel));

  cout << "There are "  << theChip->getProbes().size() 
       << " probes and " << theChip->getProbesets().size() << " probesets."  << endl;

  // debug: Write all probe data to file
  if (options.universalParam == "justDumpProbeInfo") {
    controller.dumpProbeInfoToFile("CompleteProbeInfo.dat");
    return 0;
  }

  
  // Debug: Use only probesets with certain size
  if (options.probesetSizeFilter > 0) {
    cout << "Waning: Debug filter applied.Using only probesets of size " 
         << options.probesetSizeFilter << endl;
    theChip->filterProbesets(ProbeCountEquals(options.probesetSizeFilter));
  }

  IntensityMappingPtr hookcurvePlot = controller.computeHookcurve(kRawIntensitiesOnly, hookcurveModel);
  IntensityType nsThreshold = options.uncorrectedKinkPoint;
  if (options.computeUncorrectedKinkPoint) {
    nsThreshold = 
      controller.computeKinkCoordinates("Hookplot-Primary.dat", 
                                        MathUtil::digitizeCurveNew(*hookcurvePlot), theModel).first;
  }
  cout << "Uncorrected kink point found at: " << nsThreshold << endl;

  // Compute statistics for uncorrected hookcurve
  HookcurveStatistics statsComputer(theModel, hookcurveModel, exprCalculater, options, *theChip, 
                                    controller.getAllProbesetIntensitiesArrays());
  DataCollector::instance().insert(statsComputer.computeCurveStatistics(nsThreshold, "Uncorrected"));

  // Boost does not directly support overloaded functions, hence cast is needed for exp10
  UnaryIntensityFunction invLogFunction = (IntensityTypeFunction)exp10;
  UnaryIntensityFunction invGLogFunction = &MathUtil::gExp10;

  // Compute general GGG stats
  if (1) {
    IntensityPredicate isLessThanNsThreshold = bind2nd(less<IntensityType>(), nsThreshold);
    controller.computeGGGStatistics(isLessThanNsThreshold, "gggAveragesNsPm.dat",
                                    kPmI, mem_fun_ref(&PmMmProbe::getSequence));
  }

  // Debug function to dump basic hook statistics
  if (options.universalParam == "computeUncorrectedHookStatistics") {
    statsComputer.writeProbesetDiagnosis(0.0, 0.0, "probesetStatistics");
    goto ExitProgram;
  }

  // Debug function to compute motifs statisics for GGG Paper
  if (options.universalParam == "computeMotifStats") {
    IntensityPredicate isLessThanNsThreshold = bind2nd(less<IntensityType>(), nsThreshold);
    controller.writeMotifIntensitiesDebug(isLessThanNsThreshold, "motifAveragesNsPm.dat",
                                          kPmI, mem_fun_ref(&PmMmProbe::getSequence));
    controller.writeMotifIntensitiesDebug(isLessThanNsThreshold, "motifAveragesNsMm.dat",
                                          kMmI, mem_fun_ref(&PmMmProbe::getSequenceMm));

    IntensityPredicate hasSpecificSumLogI = bind2nd(greater<IntensityType>(), nsThreshold + 0.5); 
    controller.writeMotifIntensitiesDebug(hasSpecificSumLogI, "motifAveragesSPm.dat",
                                          kPmI, mem_fun_ref(&PmMmProbe::getSequence));
    controller.writeMotifIntensitiesDebug(hasSpecificSumLogI, "motifAveragesSMm.dat",
                                          kMmI, mem_fun_ref(&PmMmProbe::getSequenceMm));

    IntensityPredicate isNotImpossiblyHigh = bind2nd(less<IntensityType>(), 20);
    controller.writeMotifIntensitiesDebug(isNotImpossiblyHigh, "motifAveragesAllPm.dat",
                                          kPmI, mem_fun_ref(&PmMmProbe::getSequence));
    controller.writeMotifIntensitiesDebug(isNotImpossiblyHigh, "motifAveragesAllMm.dat",
                                          kMmI, mem_fun_ref(&PmMmProbe::getSequenceMm));

    goto ExitProgram;
  }


  // Print uncorrected Pms
  printVectorPairsToFile(exprCalculater->calculateProbeExpressionMeasures
                         (*controller.getProbesetIntensityArray(kPmI), kNsPm,
                          mem_fun_ref(&PmMmProbe::getSequence), invGLogFunction), 
                         "ExpressionMeasurePmUncorrected_Probe.dat", "#\t" + options.intensityFilename);


  // Create Sequence profile
  controller.computeSequenceProfiles(nsThreshold, string("Primary"),
                                     options.probesetLimitForProfileComputation);
  
  printVectorPairsToFile(exprCalculater->calculateProbeExpressionMeasures
                         (*controller.getProbesetIntensityArray(kPmI), kNsPm,
                          mem_fun_ref(&PmMmProbe::getSequence), invGLogFunction), 
                         "ExpressionMeasurePmNscorrected_Probe.dat", "#\t" + options.intensityFilename);
  
  // Recalculate hookcurve: first correct the provisional intensities.
  cout << "\nPlotting new curve NsCorrected..." << endl;
  hookcurvePlot = controller.computeHookcurve(kNonspecificProfilesOnly, hookcurveModel);

  nsThreshold = options.correctedKinkPoint;
  if (options.computeCorrectedKinkPoint) {
    nsThreshold = controller.computeKinkCoordinates
      ("Hookplot-Corrected.dat", MathUtil::digitizeCurveNew(*hookcurvePlot), theModel).first;
  }
  cout << "Corrected kink point found at: " << nsThreshold << endl;

  if (options.universalParam == "compute3PrimeBias") {
    HookcurveStatistics statsComputer(theModel, hookcurveModel, exprCalculater, options, *theChip, controller.getAllProbesetIntensitiesArrays());
    statsComputer.printProbeDegration("ProbeSumLogs.dat", kSensitivityCorrectedPm, kSensitivityCorrectedMm);
    statsComputer.print3PrimeBiasAveragedHookcurve("3PrimeBias.dat", kSensitivityCorrectedPm, 
                                                   kSensitivityCorrectedMm, nsThreshold);
    goto ExitProgram;
  }

  // @debug Quit program after computation of profiles
  if (options.universalParam == "computeProfilesOnly") {
    goto ExitProgram;
  }


  controller.computeSequenceProfiles
    (nsThreshold, string("Corrected"), options.probesetLimitForProfileComputation);

  if (options.chipType != kChiptypeTiling) {
    controller.computeIntermediateSpecificProfiles
      (nsThreshold, string("Corrected"), options.probesetLimitForProfileComputation);
  }

  controller.calculateNsChip();
  // Correct the current intensities with the new profiles. (This is usually done in 
  // OligoController::computeHookcurve, but this time
  // we just need to update the intensities, and not calculate a further hookcurve.
  controller.updateCorrectedProbes(kNonspecificProfilesOnly);

  // Write Celfile with sequence corrected values if desired
  if (options.flags.count("update-celfiles")) {
    cout << "Printing Celfiles" << endl;
    string newCelfile = DataCollector::instance().getAs<string>("chipId") + "_SequenceCorrected.CEL";
    controller.writeIntensitiesToCelfile(options.intensityFilename, newCelfile, 
                                         exp10(*controller.getProbesetIntensityArray(kSensitivityCorrectedPm)),
                                         exp10(*controller.getProbesetIntensityArray(kSensitivityCorrectedMm)));
  }


  // @bug Should define these before tie (line x+10)
  // but the write statistics functions refers to these
  // (in case of no saturation they are not defined)
  IntensityType a, f;
  
  if (hookcurveModel->isSaturated(*hookcurvePlot, 30)) { // If the chip saturates we can correct the hookcurve with specific profiles.
    theModel->setImax(theModel->estimateImax(*controller.getProbesetIntensityArray(kPmI)));
    controller.calculateCompression();
  
    if (options.chipType == kChiptypeGenechip) {
      cout << "Recalculate the curve correct with specific profiles." << endl;
      hookcurvePlot = controller.computeHookcurve(kAllProfiles, hookcurveModel);

      nsThreshold = options.correctedKinkPoint;
      if (options.computeCorrectedKinkPoint) {
    	  nsThreshold 
          = controller.computeKinkCoordinates("Hookplot-Corrected.dat", MathUtil::digitizeCurveNew(*hookcurvePlot), theModel).first;
      }
      cout << "Final kink point: " << nsThreshold << endl;
    }
    else { // Tiling array: do not recompute hookcurve with specific profile
      cout << "Tiling array: do not recompute hookcurve with specific profile" << endl;
    }

    // Update once more since the intersection point has slightly changed.
    controller.updateCorrectedProbes(kAllProfiles, true);
    
    // The calculation of Imax by the Theoretic curve is optional 
    if (options.fitTheoreticHookcurve) {
      tie(f, a) = controller.computeImaxAndA(hookcurvePlot); // @note Imax is set in controller
      cout << "Calculated Theoretic Curve Parameters: a = " << log10(1/a) << ",  f = " << log10(1/f) << endl;
    }

    // Apply some scaling factor to Imax (for testing purpose, typically scaling is 1.0)
    if (options.ImaxScalingFactor != 1.0) {
      cout << "Appying a scaling factor (s=" << options.ImaxScalingFactor << ") to Imax." << endl;
      theModel->setImax(options.ImaxScalingFactor * theModel->getImax());
    }
  }
  else {
    cout << "Chip does not saturate. Skipping Saturation routines and further correction of the hookcurve." << endl;
    theModel->setImax(6.0); // else set Imax to some high value (high in log scale terms)
  }

  if (theModel->getImax() > 6) {
    cout << "We obtained an impossibly high Imax. Limiting to Imax=6.0." << endl;
    theModel->setImax(6.0); // else set Imax to some high value (high in log scale terms)
  }
 
  // hookcurveModel->setImax(4.8);

  

  DataCollector::instance().insert("saturationImax", theModel->getImax());
  DataCollector::instance().insert("saturationF", log10(1/f));
  DataCollector::instance().insert("saturationA", log10(1/a));
  
  // Correct the log intensities saved in the Probest instances by subtracting their nonspecific content.
  // (de-saturation also occurs here)
  controller.subtractNBackground(nsThreshold);
  
  controller.calculateSpecificProfiles(kSensitivityCorrectedPm, kSensitivityCorrectedMm, 
                                       kProbeSequenceLength, options.requiredProfileTypes, 
                                       options.probesetLimitForProfileComputation,
                                       nsThreshold);

  controller.calculateAndPrintExpressionMeasures(exprCalculater);
  
  // Calculating and printing statistics:
  DataCollector::instance().insert(statsComputer.computeGeneralStatistics(nsThreshold, a, f));
  DataCollector::instance().insert(statsComputer.computeCurveStatistics(nsThreshold, "Corrected"));

  cout << "Writing hook statistics to file " << flush;

  if (options.flags.count("print-probe-statistics")) {
    cout << " (incl. probe statistics)" << flush;
    statsComputer.writeProbeDiagnosis(a, f, "probeStatistics");
  }
  statsComputer.writeProbesetDiagnosis(a, f, "probesetStatistics");
  statsComputer.print3PrimeBiasAveragedHookcurve("3PrimeBias.dat", kSensitivityCorrectedPm, 
                                                 kSensitivityCorrectedMm, nsThreshold);
//   statsComputer.writeProbesetVariance("probesetVariance");
  cout << "Done." << endl;

  // Write Celfile with Expression Measures if desired
  if (options.flags.count("update-celfiles")) {
    if (exprCalculater->isEmPresent(kEmPmSingleIntegrationExp10_Probe) &&
        exprCalculater->isEmPresent(kEmMmSingleIntegrationExp10_Probe)) {
      cout << "Printing Celfile" << endl;
      string newCelfile = DataCollector::instance().getAs<string>("chipId") + "_ExpressionMeasure.CEL";
      controller.writeIntensitiesToCelfile(options.intensityFilename, newCelfile, 
                                           exprCalculater->getEmsArray(kEmPmSingleIntegrationExp10_Probe), 
                                           exprCalculater->getEmsArray(kEmMmSingleIntegrationExp10_Probe));
    }
  }

 ExitProgram:
  cout << endl;


  // That's the last thing we do (maybe we can put this in the destructor of DataCollector ?)
  DataCollector::instance().writeToFile("dataLogger.log");
  return 0;
}


