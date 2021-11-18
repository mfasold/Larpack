#include <string>
#include <vector>
#include "boost/tuple/tuple.hpp"
#include "MicroarrayExperiment.hpp"
#include "DataCollector.hpp"
#include "Probeset.hpp"
#include "HookcurveAnalyzer.hpp"
#include "AveragingProcedure.hpp"
#include "ProbesetComposition.hpp"
#include "DetectKinkPoint.hpp"
#include "ProbeFilter.hpp"
#include "OligoController.hpp"
#include "StlUtil.hpp"
#include "MathUtil.hpp"
#include "HookcurveStatistics.hpp"


extern "C" {
  
# include <Rinternals.h>
# include <Rdefines.h>

// define personal PROTECT macro that counts number of calls
#define _PROTECT(X) PROTECT(X); ++protectCount;


  using namespace larrpack;
  using namespace std;
  using namespace boost;

  inline size_t xyToIndex(const std::pair<int,int> xy, const size_t ncol) {
    return  xy.first + xy.second*ncol;
  }

  SEXP larpackOligo (SEXP celFile, SEXP seqFile, SEXP optionsList,
                     SEXP intensities, SEXP nrow, SEXP ncol) 
  {
    const char* celFileName         = CHAR(STRING_ELT(celFile,0));
    const char* seqFileName         = CHAR(STRING_ELT(seqFile,0));
    
    // Build up argc, argv data structures
    vector<string> argvList(0); 
    argvList.push_back("arrayCorrect"); 
    argvList.push_back("-i"); argvList.push_back(string(celFileName));
    argvList.push_back("-c"); argvList.push_back(string(seqFileName));

    // Add options from optionsList to argc, argv
    SEXP listElement = R_NilValue, listNames = getAttrib(optionsList, R_NamesSymbol);
    for (size_t i = 0; i < (size_t) length(optionsList); ++i) {
      string currentParamName = string(CHAR(STRING_ELT(listNames, i)));
      listElement = VECTOR_ELT(optionsList, i);
      string currentParamValue = string(CHAR(STRING_ELT(listElement,0)));
      argvList.push_back(string("--") + currentParamName);
      argvList.push_back(currentParamValue);
      // cout << currentParamName << "\t" << currentParamValue << endl;
    }

    // "Convert" to cstring
    char* argv[argvList.size()];
    for (size_t i = 0; i < argvList.size(); ++i) {
      argv[i] = const_cast<char*>(argvList[i].c_str()); // cast const away
    }

    // Read options from File/Commandline
    ProgramOptions options = parseCommandlineOptions(argvList.size(), argv);

    // Import expression values from R if desired
    IntensityMatrixPtr chipIntensities; // initialize as empty pointer
    BoolArrayPtr chipMaskedIntensities;
    if (string(celFileName) == "/useAffybatch") {
      size_t chipRows = (size_t)*INTEGER(nrow);
      size_t chipCols = (size_t)*INTEGER(ncol);
      chipIntensities = IntensityMatrixPtr(new IntensityMatrix(chipRows, chipCols));
      chipMaskedIntensities = BoolArrayPtr(new BoolArray(chipRows, chipCols, false));
      
      size_t intensityIndex = 0;
      for (size_t col = 0; col < chipCols; ++col) { // intensities are stored in column-major order
        for (size_t row = 0; row < chipRows; ++row) { 
          // intensityIndex = row + col*chipRows
          (*chipIntensities)[row][col] = REAL(intensities)[intensityIndex++];
        }
      }      
    }

    // Import Files
    OligoController controller(options);
    //    controller.logInitialOptions();
    PmMmProbePtrVector probes = controller.importProbeData(chipIntensities, chipMaskedIntensities);
    ChipPtr theChip = controller.createChip(probes);
    HookModelPtr theModel = controller.createHookModel();
    controller.initializeIntensityArrays();
    HookcurveAnalyzerPtr hookcurveModel = controller.createHookcurveModel();
    ExpressionMeasurePtr exprCalculater = ExpressionMeasurePtr(new ExpressionMeasure(*theChip, *theModel));

    // here begins arrayCorrect.cpp excerpt
    // -----------------------------------------------------------------------------------
    cout << "There are "  << theChip->getProbes().size() 
         << " probes and " << theChip->getProbesets().size() << " probesets."  << endl;
  
    // Debug: Use only probesets with certain size
    if (options.probesetSizeFilter > 0) {
      cout << "Waning: Debug filter applied.Using only probesets of size " 
           << options.probesetSizeFilter << endl;
      theChip->filterProbesets(ProbeCountEquals(options.probesetSizeFilter));
    }

    IntensityMappingPtr hookcurvePlot = controller.computeHookcurve(kRawIntensitiesOnly, hookcurveModel);
    IntensityMappingPtr hookcurvePlotRaw = hookcurvePlot;
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

    // boost does not directly support overloaded functions, hence cast is needed for exp10
    UnaryIntensityFunction invLogFunction = (IntensityTypeFunction)exp10;
    UnaryIntensityFunction invGLogFunction = &MathUtil::gExp10;

    // Print uncorrected Pms
    printVectorPairsToFile(exprCalculater->calculateProbeExpressionMeasures(*controller.getProbesetIntensityArray(kPmI), kNsPm,
                                                                            mem_fun_ref(&PmMmProbe::getSequence), invGLogFunction), "ExpressionMeasurePmUncorrected_Probe.dat", 
                           "#\t" + options.intensityFilename);

    // Create Sequence profile
    controller.computeSequenceProfiles(nsThreshold, string("Primary"),
                                       options.probesetLimitForProfileComputation);

    printVectorPairsToFile(exprCalculater->calculateProbeExpressionMeasures(*controller.getProbesetIntensityArray(kPmI), kNsPm,
                                                                            mem_fun_ref(&PmMmProbe::getSequence), invGLogFunction), 
                           "ExpressionMeasurePmNscorrected_Probe.dat", 
                           "#\t" + options.intensityFilename);
  
    // Recalculate hookcurve: first correct the provisional intensities.
    cout << "\nPlotting new curve NsCorrected..." << endl;
    hookcurvePlot = controller.computeHookcurve(kNonspecificProfilesOnly, hookcurveModel);
    IntensityMappingPtr hookcurvePlotCorrected = hookcurvePlot;
    vector<IntensityType>  correctedSumLogs    = controller.getAverageSumLogIArray(kSensitivityCorrectedPm, kSensitivityCorrectedMm);  

    nsThreshold = options.correctedKinkPoint;
    if (options.computeCorrectedKinkPoint) {
      nsThreshold = controller.computeKinkCoordinates("Hookplot-Corrected.dat", MathUtil::digitizeCurveNew(*hookcurvePlot), theModel).first;
    }
    cout << "Corrected kink point found at: " << nsThreshold << endl;


    controller.computeSequenceProfiles(nsThreshold, string("Corrected"),options.probesetLimitForProfileComputation);

    if (options.chipType != kChiptypeTiling) {
      controller.computeIntermediateSpecificProfiles(nsThreshold, string("Corrected"), options.probesetLimitForProfileComputation);
    }

    if (options.universalParam == "computeProfilesOnly") {
      cout << "Exiting" << endl;
      return 0;
    }

    controller.calculateNsChip();
    // Correct the current intensities with the new profiles. (This is usually done in OligoController::computeHookcurve, but this time
    // we just need to update the intensities, and not calculate a further hookcurve.
    controller.updateCorrectedProbes(kNonspecificProfilesOnly);
  
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
    }
    else {
      cout << "Chip does not saturate. Skipping Saturation routines and further correction of the hookcurve." << endl;
      theModel->setImax(6.0); // else set Imax to some high value (high in log scale terms)
    }
 
    DataCollector::instance().insert("saturationImax", theModel->getImax());
    DataCollector::instance().insert("saturationF", log10(1/f));
    DataCollector::instance().insert("saturationA", log10(1/a));
  
    IntensityPredicate isLessThanIntersectionSumLogI = bind2nd(less<IntensityType>(), nsThreshold);
    IntensityPredicate hasSpecificSumLogI = bind2nd(greater<IntensityType>(), nsThreshold);
  
    // Correct the log intensities saved in the Probest instances by subtracting their nonspecific content.
    controller.subtractNBackground(nsThreshold);
  

    controller.calculateSpecificProfiles(kSensitivityCorrectedPm, kSensitivityCorrectedMm, 
                                         kProbeSequenceLength, options.requiredProfileTypes, 
                                         options.probesetLimitForProfileComputation,
                                         nsThreshold);

    controller.calculateAndPrintExpressionMeasures(exprCalculater);
  
    // Calculating and printing statistics:
    DataCollector::instance().insert(statsComputer.computeGeneralStatistics(nsThreshold, a, f));
    DataCollector::instance().insert(statsComputer.computeCurveStatistics(nsThreshold, "Corrected"));

    // cout << "Writing hook statistics to file " << flush;
    // statsComputer.writeProbeDiagnosis(a, f, "probeStatistics", "\t");
    // statsComputer.writeProbesetDiagnosis(a, f, "probesetStatistics", "\t");
    // statsComputer.print3PrimeBiasAveragedHookcurve("3PrimeBias.dat", kSensitivityCorrectedPm, kSensitivityCorrectedMm, nsThreshold); 

    cout << endl;
  
    // That's the last thing we do (maybe we can put this in the destructor of DataCollector ?)
    DataCollector::instance().writeToFile("dataLogger.log");


    // here ends arrayCorrect.cpp excerpt
    // -----------------------------------------------------------------------------------
    // Export results to R
    SEXP  pSetNames, expr, nsFraction,
      hookcurveSumRaw, hookcurveDeltaRaw,
      hookcurveSumFinal, hookcurveDeltaFinal,
      interSection, dimnames, result_list,
      correctedSumLogsR,
      exprMeasuresPm, exprMeasuresFastDiff, exprMeasuresGcrmalikeDiff;
    size_t protectCount = 0;
    size_t exprMeasuresCount = 0;
    

    // lets assume that nProbeSets holds the number of probe sets,
    // eM is an array of expression measures an pSNames an array of 
    // probeset names
    vector<Probeset> probesets = theChip->getProbesets();
    size_t probesetCount = probesets.size();

    // allocate two vectors for the names of the probe sets and 
    // for the expression Measures
    _PROTECT(pSetNames    = allocVector (STRSXP, probesetCount));
    for (size_t i = 0; i < probesetCount; ++i) {
      SET_STRING_ELT (pSetNames, i, mkChar(probesets[i].getProbesetId().c_str()));       
    }

    if (exprCalculater->isEmPresent(kEmPmSingleIntegrationExp10) ) {
      _PROTECT(exprMeasuresPm = allocVector (REALSXP, probesetCount));
      ++exprMeasuresCount;
      for (size_t i = 0; i < probesetCount; ++i) {
        REAL (exprMeasuresPm)[i] = computeAverageIntensity(exprCalculater->getEms(kEmPmSingleIntegrationExp10, i));
      }
    }

    if (exprCalculater->isEmPresent(kEmDeltaSingleIntegrationExp10) ) {
      _PROTECT(exprMeasuresFastDiff = allocVector (REALSXP, probesetCount));
      ++exprMeasuresCount;
      for (size_t i = 0; i < probesetCount; ++i) {
        REAL (exprMeasuresFastDiff)[i] = computeAverageIntensity(exprCalculater->getEms(kEmDeltaSingleIntegrationExp10, i));
      }      
    }
    
    if (exprCalculater->isEmPresent(kEmDeltaGcrmaLikeIntegrationExp10) ) {
      _PROTECT(exprMeasuresGcrmalikeDiff = allocVector (REALSXP, probesetCount));
      ++exprMeasuresCount;
      for (size_t i = 0; i < probesetCount; ++i) {
        REAL (exprMeasuresGcrmalikeDiff)[i] = computeAverageIntensity(exprCalculater->getEms(kEmDeltaGcrmaLikeIntegrationExp10, i));
      }      
    }

    // Export probe-expression measures in the same format as affybatch stores exprs
    size_t chipRows = DataCollector::instance().getAs<size_t>("chipRowNumber"); // (size_t)*INTEGER(nrow);
    size_t chipCols = DataCollector::instance().getAs<size_t>("chipColumnNumber"); // (size_t)*INTEGER(ncol);
    _PROTECT(expr = allocVector (REALSXP, chipRows * chipCols));

    // Sensitivity-corrected intensities only:
    // const ProbesetIntensitiesArray pmExpression = exp10(*controller.getProbesetIntensityArray(kSensitivityCorrectedPm));
    // const ProbesetIntensitiesArray mmExpression = exp10(*controller.getProbesetIntensityArray(kSensitivityCorrectedMm));   
    // const ProbesetIntensitiesArray& pmExpression = exprCalculater->getEmsArray(kEmDeltaSingleIntegrationExp10_Probe);
    // const ProbesetIntensitiesArray& mmExpression = exprCalculater->getEmsArray(kEmDeltaGcrmaLikeIntegrationExp10_Probe);

    const ProbesetIntensitiesArray& pmExpression = exprCalculater->getEmsArray(kEmPmSingleIntegrationExp10_Probe);
    const ProbesetIntensitiesArray& mmExpression = exprCalculater->getEmsArray(kEmMmSingleIntegrationExp10_Probe);

    for (size_t probesetIndex = 0; probesetIndex < probesetCount; ++probesetIndex) {
      for (size_t probeIndex = 0; probeIndex < probesets[probesetIndex].getSize(); ++probeIndex) {
        PmMmProbePtr currentProbePtr = probesets[probesetIndex].getProbePtr(probeIndex);
        REAL(expr)[xyToIndex(currentProbePtr->getPositionPm(), chipCols)] = pmExpression[probesetIndex][probeIndex];
        REAL(expr)[xyToIndex(currentProbePtr->getPositionMm(), chipCols)] = mmExpression[probesetIndex][probeIndex];
      }
    }

    // Export uncorrected hookcurve
    _PROTECT(hookcurveSumRaw = allocVector (REALSXP, hookcurvePlot->size()));
    _PROTECT(hookcurveDeltaRaw = allocVector (REALSXP, hookcurvePlot->size()));
    for (size_t i = 0; i < hookcurvePlot->size(); ++i) {
      REAL (hookcurveSumRaw)[i] = (*hookcurvePlotRaw)[i].first;
      REAL (hookcurveDeltaRaw)[i] = (*hookcurvePlotRaw)[i].second;
    }

    // Export corrected hookcurve
    _PROTECT(hookcurveSumFinal = allocVector (REALSXP, hookcurvePlotCorrected->size()));
    _PROTECT(hookcurveDeltaFinal = allocVector (REALSXP, hookcurvePlotCorrected->size()));
    for (size_t i = 0; i < hookcurvePlotCorrected->size(); ++i) {
      REAL (hookcurveSumFinal)[i] = (*hookcurvePlotCorrected)[i].first;
      REAL (hookcurveDeltaFinal)[i] = (*hookcurvePlotCorrected)[i].second;
    }

    // Export corrected hookcurve
    _PROTECT(correctedSumLogsR = allocVector (REALSXP, correctedSumLogs.size()));
    for (size_t i = 0; i < correctedSumLogs.size(); ++i) {
      REAL (correctedSumLogsR)[i] = correctedSumLogs[i];      
    }

    // Allocate a 1 element real vector for the intersection point
    _PROTECT(interSection=allocVector (REALSXP, 1));
    REAL (interSection)[0] = nsThreshold; 

    // Allocate a 1 element real vector for the percentage of non-specific probes
    _PROTECT(nsFraction = allocVector (REALSXP, 2));
//     REAL (nsFraction)[0] = DataCollector::instance().getAs<IntensityType>("probeNsFractionUncorrected"); 
//     REAL (nsFraction)[1] = DataCollector::instance().getAs<IntensityType>("probesetNsFractionUncorrected"); 
    REAL (nsFraction)[0] = DataCollector::instance().getAs<IntensityType>("probeNsFractionCorrected"); 
    REAL (nsFraction)[1] = DataCollector::instance().getAs<IntensityType>("probesetNsFractionCorrected"); 

    // We use a list of key value pairs to return data
    _PROTECT(result_list=NEW_LIST(9 + exprMeasuresCount));

    // Now fill the result list with the stuff
    SET_VECTOR_ELT (result_list, 0, pSetNames);
    SET_VECTOR_ELT (result_list, 1, hookcurveSumRaw);
    SET_VECTOR_ELT (result_list, 2, hookcurveDeltaRaw);
    SET_VECTOR_ELT (result_list, 3, hookcurveSumFinal);
    SET_VECTOR_ELT (result_list, 4, hookcurveDeltaFinal);
    SET_VECTOR_ELT (result_list, 5, expr);
    SET_VECTOR_ELT (result_list, 6, nsFraction);
    SET_VECTOR_ELT (result_list, 7, interSection);    
    SET_VECTOR_ELT (result_list, 8, correctedSumLogsR);
    SET_VECTOR_ELT (result_list, 9, exprMeasuresPm);
    if (exprCalculater->isEmPresent(kEmDeltaSingleIntegrationExp10) ) {
      SET_VECTOR_ELT (result_list, 10, exprMeasuresFastDiff);
    }
    if (exprCalculater->isEmPresent(kEmDeltaGcrmaLikeIntegrationExp10) ) { // @bug: kEmDeltaSingleIntegrationExp10 must appear together with kEmDeltaSingleIntegrationExp10
      SET_VECTOR_ELT (result_list, 11, exprMeasuresGcrmalikeDiff);
    }

    // Name the elements of the list
    _PROTECT(dimnames = allocVector(VECSXP, 9 + exprMeasuresCount));
    SET_VECTOR_ELT (dimnames, 0, mkChar("ProbesetId"));
    SET_VECTOR_ELT (dimnames, 1, mkChar("RawHookcurveSumLogI"));
    SET_VECTOR_ELT (dimnames, 2, mkChar("RawHookcurveDeltaLogI"));
    SET_VECTOR_ELT (dimnames, 3, mkChar("CorrectedHookcurveSumLogI"));
    SET_VECTOR_ELT (dimnames, 4, mkChar("CorrectedHookcurveDeltaLogI"));
    SET_VECTOR_ELT (dimnames, 5, mkChar("ProbeExpressionMeasures"));
    SET_VECTOR_ELT (dimnames, 6, mkChar("FractionOfNonSpecificBinding"));
    SET_VECTOR_ELT (dimnames, 7, mkChar("nsThreshold"));
    SET_VECTOR_ELT (dimnames, 8, mkChar("CorrectedSumLogI"));
    SET_VECTOR_ELT (dimnames, 9, mkChar("ExpressionMeasurePm"));
    if (exprCalculater->isEmPresent(kEmDeltaSingleIntegrationExp10) ) {
      SET_VECTOR_ELT (dimnames, 10, mkChar("ExpressionMeasureDifference"));
    }
    if (exprCalculater->isEmPresent(kEmDeltaGcrmaLikeIntegrationExp10) ) {
      SET_VECTOR_ELT (dimnames, 11, mkChar("ExpressionMeasureGcrmalikeDifference"));
    }

    setAttrib(result_list, R_NamesSymbol, dimnames);

    // Clear singleton for next larpack run
    DataCollector::instance().clear();

    UNPROTECT (protectCount);
    return (result_list);
  }

  /**
   * Wrapper for methods to calculate the boundary between specific 
   * and non-specific binding (nsThreshold)
   *
   */
  SEXP larpackComputeNsThreshold (SEXP x, SEXP y, SEXP method, SEXP digitizeIntervals)
  {
    assert(length(x) == length(y));

    // Transform vectors x,y to internal graph data structure
    IntensityMapping graph(length(x));
    for (int i = 0; i < length(x); i++) {
	    graph[i].first = REAL(x)[i];
	    graph[i].second = REAL(y)[i];
    }

    // Compute NS threshold 
    DetectKinkPointPtr detectKinkPoint;
    switch(INTEGER(method)[0]) {
    case 0:
      detectKinkPoint.reset(new FitStraightLineAndParabola()); break;
    case 1:
      detectKinkPoint.reset(new FitTwoStraightLines(5, 0.6)); break;
    case 2:
      detectKinkPoint.reset(new FirstDerivativeAnalysis());break;
    case 3:
      detectKinkPoint.reset(new SecondDerivativeAnalysis()); break;
    default:
      error("Valid values for method are 0-3"); 
    }
    IntensityPair kinkPoint 
      = (*detectKinkPoint)(MathUtil::digitizeCurve(graph, INTEGER(digitizeIntervals)[0]));

    // Return pair elements as R list
    SEXP nsThreshold, nsDelta, dimnames, result_list;
    PROTECT (nsThreshold = allocVector (REALSXP, 1));
    REAL (nsThreshold)[0] = kinkPoint.first;

    PROTECT (nsDelta = allocVector (REALSXP, 1));
    REAL (nsDelta)[0] = kinkPoint.second;

    PROTECT (result_list=NEW_LIST(2));
    SET_VECTOR_ELT (result_list, 0, nsThreshold);
    SET_VECTOR_ELT (result_list, 1, nsDelta);

    PROTECT (dimnames = allocVector(VECSXP, 2));
    SET_VECTOR_ELT (dimnames, 0, mkChar("nsThreshold"));
    SET_VECTOR_ELT (dimnames, 1, mkChar("nsDelta"));
    setAttrib(result_list, R_NamesSymbol, dimnames);

		UNPROTECT(4);
    return(result_list);    
  }


  /**
   * Wrapper for method to a sequence profile
   *
   * !!!!! probesets must be sorted !!!!
   */
  SEXP larpackComputeSequenceProfile (SEXP sequences, SEXP probesets, SEXP intensities, SEXP modelRank)
  {
    assert(length(sequences) == length(probesets));
    assert(length(sequences) == length(intensities));
 
    string oneSequence(CHAR(STRING_ELT(sequences,1)));
    size_t probeSequenceLength = oneSequence.size();
    cout << probeSequenceLength << endl;

    // Transform input vectors to internal graph data structure  
    ProbesetIntensitiesArray intensitiesArray; 
    ProbesetSequencesArray   probesetSequences;

    string lastProbesetId( CHAR(STRING_ELT(probesets, 0) )); // ----------- TEST ---------------
    // cout <<  lastProbesetId << "\n";
    std::vector<std::string> currentProbesetSequences;
    std::vector<IntensityType> currentProbesetIntensities;

    for (int i = 0; i < length(sequences); ++i) {
      string currentProbesetId( CHAR(STRING_ELT(probesets, i)) );
      if (currentProbesetId != lastProbesetId) { // Store probeset data to larger array if new probeset begins
        probesetSequences.push_back(currentProbesetSequences);
        intensitiesArray.push_back(convertVectorToValarray(currentProbesetIntensities));

        currentProbesetSequences.clear();
        currentProbesetIntensities.clear();            
        lastProbesetId = currentProbesetId;
      }

      // Append probe data to probeset
      currentProbesetSequences.push_back(string( CHAR(STRING_ELT(sequences,i)) ));
      currentProbesetIntensities.push_back( REAL(intensities)[i] ); // Can be sped up by saving double *int = REAL(intensities)
    }

    // Append leftover probeset to list
    if (!currentProbesetSequences.empty()) {
      probesetSequences.push_back(currentProbesetSequences);
      intensitiesArray.push_back(convertVectorToValarray(currentProbesetIntensities));
    }

    // DEBUG
    // printValarray(cout, intensitiesArray[0]);

    // Create proper profile
    std::vector<std::string> gstackTuples; // initialize empty tuples
    size_t profileModelRank = INTEGER(modelRank)[0];
    size_t probesetLimitForProfileComputation = 5000; // 5000;
    SensitivityProfilePtr sensitivityProfileNSPm 
      = ProfileFactory::createNsProfile(probeSequenceLength, profileModelRank, gstackTuples);

    // Compute
    sensitivityProfileNSPm->computeProfile(probesetSequences, intensitiesArray, 
                                           probesetLimitForProfileComputation); 
  

    SimpleSensitivityProfilePtr simpleProfileNsPm 
      = boost::shared_static_cast<SimpleSensitivityProfile>(sensitivityProfileNSPm);

    // Return pair elements as R list
    SEXP profile, profileNames, tripletProfile, tripletProfileNames, correctedIntensities, dimnames, result_list;
    size_t protectCount = 0;

    // Store sequence profile
    size_t profileSize = simpleProfileNsPm->getProfileSize(profileModelRank);
    IntensityArray seqProfile = simpleProfileNsPm->getProfile();
    _PROTECT (profile = allocVector (REALSXP, profileSize));
    _PROTECT (profileNames = allocVector (STRSXP, profileSize)); 
    for (size_t i = 0; i < profileSize; ++i) { // can be optimized - copy in one run
      REAL (profile)[i] = seqProfile[i];    
      size_t seqPos = (i % simpleProfileNsPm->getValidSequenceCount());
      size_t baseIndex = size_t (i / simpleProfileNsPm->getValidSequenceCount());
      string baseString = sequenceutil::indexToSequence(baseIndex, profileModelRank);
      string tupleName = baseString + "." + boost::lexical_cast<std::string>(seqPos);
      SET_STRING_ELT (profileNames, i, mkChar(tupleName.c_str()));      
    }
    setAttrib(profile, R_NamesSymbol, profileNames);    

    // printValarray(cout, seqProfile); // DEBUG: Print array

    // For SNP: store Triplet profile
    size_t validPositions = (probeSequenceLength - 2);
    _PROTECT (tripletProfile = allocVector (REALSXP, 64*validPositions)); // TODO: remove constants
    _PROTECT (tripletProfileNames = allocVector (STRSXP, 64*validPositions)); 
    // _PROTECT (tripletProfileNames = allocVector (VECSXP, 64*23)); 
    if (profileModelRank == 2) { // now only works for model rank 2 -> make the function work for other ranks as well
      for (size_t baseIndex = 0; baseIndex < 64; ++baseIndex) {    
        for(size_t pos = 0; pos < validPositions; ++pos) {
          string baseString = sequenceutil::indexToSequence(baseIndex, 3);
          REAL (tripletProfile)[baseIndex*validPositions + pos] = simpleProfileNsPm->getSubTripleIncrement(baseString, pos);
          SET_STRING_ELT (tripletProfileNames, baseIndex*validPositions + pos, mkChar(baseString.c_str()));
          // SET_VECTOR_ELT (tripletProfileNames, baseIndex*24 + pos, mkChar("AAA")); // baseString.c_str()));
          // sensitivityProfileFile << pos + 1 << "\t" << baseString << "\t" << 
          //   simpleProfileNSPm->getSubTripleIncrement(baseString, pos) << endl;
        }
      }
    }
    setAttrib(tripletProfile, R_NamesSymbol, tripletProfileNames);    

    // Store corrected intensities
    _PROTECT (correctedIntensities = allocVector (REALSXP, length(intensities)));
    for (size_t i = 0; i < length(intensities); ++i) { 
      REAL (correctedIntensities)[i] = REAL(intensities)[i] // speed up: save pointer to REAL() 's 
        + sensitivityProfileNSPm->getSequenceIncrement(string( CHAR(STRING_ELT(sequences,i)) ));  
    }
        
    _PROTECT (result_list=NEW_LIST(3));
    SET_VECTOR_ELT (result_list, 0, profile);
    SET_VECTOR_ELT (result_list, 1, tripletProfile);
    SET_VECTOR_ELT (result_list, 2, correctedIntensities);

    _PROTECT (dimnames = allocVector(VECSXP, 3));
    SET_VECTOR_ELT (dimnames, 0, mkChar("affinities"));
    SET_VECTOR_ELT (dimnames, 1, mkChar("triplet.contributions"));
    SET_VECTOR_ELT (dimnames, 2, mkChar("intensities.corrected"));
    setAttrib(result_list, R_NamesSymbol, dimnames);

    UNPROTECT(protectCount);
    return(result_list);    
  }
}
