/**
 * @file OligoController.hpp Defines a controller to perform various steps of the OLIGO method
 * @author $Author: mario $ 
 * @author Mario Fasold
 * @date $Date: 2007-04-20 16:55:02 +0200 (Fri, 20 Apr 2007) $
 */
#ifndef _OLIGOCONTROLLER_
#define _OLIGOCONTROLLER_

#include <boost/scoped_ptr.hpp>
#include <Probeset.hpp>
#include <ExpressionMeasure.hpp>
#include <set>

#include "BackgroundSubtraction.hpp"
#include "DetectKinkPoint.hpp"
#include "HookcurveAnalyzer.hpp"
#include "ValarrayUtil.hpp"
#include "ProgramOptions.hpp"

//#include "Probe.hpp"

namespace larrpack {
  
  /**
   * Controller to perform various steps of the OLIGO method
   *
   *
   */
  class OligoController
  {
  public:
    OligoController(const ProgramOptions opt) : options(opt), mCorrectedIntensities(kIntensityModeCount) {}
    PmMmProbePtrVector importProbeData(IntensityMatrixPtr intensityArray = IntensityMatrixPtr(), 
                                       BoolArrayPtr intensityMaskArray = BoolArrayPtr());
    ChipPtr createChip(PmMmProbePtrVector& probes); // Set up probesets and Chip object
    
    void initializeIntensityArrays();
    
    ProbesetIntensitiesArrayPtr createPseudoMmArray(const ProbesetIntensitiesArray& pmArray, const std::vector< std::vector<size_t> >& gcToIndexMap);
    void initializeIntensityArraysPmOnly(const std::vector< std::vector<size_t> >& gcToIndexMap);
    
    //IntensityArray& getLogIs(IntensityMode t, size_t probesetIndex) const;
    
    void updateCorrectedProbes(const HookcurveCalculationType& calculationType,  bool print = false);
    
    void updateCorrectedProbesPmOnly(const HookcurveCalculationType& calculationType,  bool print = false);
    
    void calculateNsChip();
    
    void calculateNsChipPmOnly();
    
//    void exportIntensityArraysTofile(std::string f);

    HookModelPtr createHookModel();
    
    HookcurveAnalyzerPtr createHookcurveModel(); // This is temprerilly, should be removed or altered when Chip and Model work
    
    IntensityMappingPtr computeHookcurve(const HookcurveCalculationType& calculationType, 
                                         HookcurveAnalyzerPtr hookcurveModel);
    
    IntensityMappingPtr computeHookcurvePmOnly(const HookcurveCalculationType& calculationType, 
                                         HookcurveAnalyzerPtr hookcurveModel, const std::vector< std::vector<size_t> >& gcToIndexMap);
    
    IntensityMappingPtr calculateGCContant(IntensityMode pm, IntensityMode mm, HookcurveAnalyzerPtr hookcurveAnalyzer);
    
    IntensityPair computeKinkCoordinates(const std::string hookplotFilename, const IntensityMapping& hookcurvePlot,
                                         HookModelPtr hookModel);

    void calculateCompression();
    
    void computeIntermediateSpecificProfiles(const IntensityType nsThreshold,
                                             const std::string filenameSuffix, const size_t probesetLimit = 0);
    
    void computeIntermediateSpecificProfilesPmOnly(const IntensityType nsThreshold,
                                                              const std::string filenameSuffix, 
                                                              const size_t probesetLimit = 0);

    void computeSequenceProfiles(const IntensityType nsThreshold, const std::string filenameSuffix, 
                                 const size_t probesetLimit = 0);
    
    void computeSequenceProfilesPmOnlyChip(IntensityMode pmIntensityType, IntensityMode mmIntensityType, const IntensityType nsThreshold, const std::string filenameSuffix, 
    		/*bool renameBackgroundProbeIds,*/ const size_t probesetLimit = 0);
    
    boost::tuple<IntensityType, IntensityType>  computeImaxAndA(const IntensityMappingPtr hookcurvePlot);
    
    
    /// Calculates the specific sensitivity profiles that are needed.
    void calculateSpecificProfiles(IntensityMode pmMode, IntensityMode mmMode, const size_t& sequenceLength,
                                   const std::vector<ExpressionMeasureProfileType> profileTypes,
                                   const size_t& probesetLimit, const IntensityType& nsThreshold);

    /// Writes various options to the logger
    void logInitialOptions();
    
    void subtractNBackground(const IntensityType nsThreshold);
    
    void subtractNBackgroundPmOnly(const IntensityType nsThreshold);
    
    void inverseWashingOfIntensities();
    
    void calculateAndPrintExpressionMeasures(ExpressionMeasurePtr exprCalculater);

    void writeIntensitiesToCelfile(const std::string sourceFilename, 
                                   const std::string targetFilename,
                                   const ProbesetIntensitiesArray& intensitiesPm,
                                   const ProbesetIntensitiesArray& intensitiesMm);


    void writeDebuginfoToCelfile(const std::string sourceFilename, 
                                 const std::string targetFilename);
    void writeMotifIntensitiesDebug(const IntensityPredicate considerProbeset, 
                                    const std::string targetFilename,
                                    const IntensityMode intensityMode, 
                                    const PmMmProbeSequenceFunction& getProbeSequence);
    
    inline ProbesetIntensitiesArrayPtr getProbesetIntensityArray(IntensityMode mode) 
    { return mCorrectedIntensities[mode]; }
    
    IntensityType calculateImaxPmOnlyArray();

    inline ProbesetIntensitiesArray getSumLogIArray(IntensityMode mode1, IntensityMode mode2)
    { return (*mCorrectedIntensities[mode1] + *mCorrectedIntensities[mode2]) / 2.0; }

    inline ProbesetIntensitiesArray getDeltaLogIArray(IntensityMode mode1, IntensityMode mode2)
    { return (*mCorrectedIntensities[mode1] - *mCorrectedIntensities[mode2]); }


    std::vector<IntensityType> getAverageSumLogIArray(IntensityMode mode1, IntensityMode mode2);
    
    // @note: returns a reference! turn this to a pointer.
    const std::vector<ProbesetIntensitiesArrayPtr>& getAllProbesetIntensitiesArrays();// { return mCorrectedIntensities;}

  private:
    /// The program options
    ProgramOptions options;

    /// All probestes with sumLogI > nsThreshold + kSpecificThresholdIncrement bind predominantly specific
    static const IntensityType kSpecificThresholdIncrement;
	  
    /// Contains corrected intensities for different methods
    std::vector<ProbesetIntensitiesArrayPtr> mCorrectedIntensities;
    
    // Reference to the model
    HookcurveAnalyzerPtr mHookcurveModel;

    // Reference to the mathematical model
    HookModelPtr mHookModel;
    
    // Reference to the Oligo Chip
    ChipPtr mChip;
  };

  /// Sets default options and parses them from command line stings
  ProgramOptions parseCommandlineOptions(int argc, char *argv[]);

  /**
   * Returns a new sequence profile according to previously given set of paramters.
   * Uses the factory method design pattern.
   *
   */
  class ProfileFactory {
  public:
    static const size_t kMinimumSpecificProbesetsForNn = 4000;
    static const size_t kMinimumSpecificProbesetsForN  = 1000;
    static SensitivityProfilePtr createNsProfile(const size_t sequenceLength, const unsigned short modelRank, 
                                                 const std::vector<std::string> gstackTuples 
                                                 = std::vector<std::string>());
    static SensitivityProfilePtr createSpecificProfile(const size_t sequenceLength, 
                                                       const unsigned short modelRank, 
                                                       const std::vector<std::string> gstackTuples 
                                                       = std::vector<std::string>());
    static size_t getMinimumProbesetCountToEstimateModel(const size_t modelRank);
  };
}

#endif

