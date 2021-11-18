/**
 * @file HookcurveAnalyzer.hpp Provides methods for analysis of PmMmProbes and -probesets.
 * @author $Author$ 
 * @author Mario Fasold
 * @date $Date$
 */
#ifndef _ANALYZEDMICROARRAY_
#define _ANALYZEDMICROARRAY_

#include <string>
#include <vector>
#include <map>
#include <valarray>
#include <boost/function.hpp>
#include "boost/tuple/tuple.hpp"

#include "PmMmProbe.hpp"
#include "Probeset.hpp"
#include "SensitivityProfile.hpp"
#include "Chip.hpp"
#include "HookModel.hpp"
//#include "OligoController.hpp"
//#include "ExpressionMeasure.hpp"

#include <functional>


namespace larrpack {
  class HookcurveAnalyzer; // Forward definition to declare the next type

  /// Pointer to an hookcurve analyzer
  typedef boost::shared_ptr<HookcurveAnalyzer> HookcurveAnalyzerPtr;

  /// Const pointer to an hookcurve analyzer
  typedef boost::shared_ptr<const HookcurveAnalyzer> HookcurveAnalyzerConstPtr;

  /// A pair of probeset identifier and its expression value
  typedef std::pair<std::string, IntensityType> ProbesetExpressionPair;
  
  
  


/*  /// A function object to create the probeset composition
  typedef boost::function<ProbesetVector (PmMmProbePtrVector)> ProbesetCompositionFunction;*/
  // Not all compilers like the declaration above. Alternatively, try:
  // typedef boost::function1<ProbesetVector, PmMmProbePtrVector> ProbesetCompositionFunction;

/*  /// A function object to create the probeset composition
  typedef boost::function<bool (Probeset)> ProbesetPredicate;
  // Not all compilers like the declaration above. Alternatively, try:
  // typedef boost::function1<bool, Probeset> ProbesetPredicate;
*/  
  typedef boost::function<IntensityType (IntensityType)> UnaryIntensityFunction;
  
  typedef IntensityType (*IntensityTypeFunction)(IntensityType);
  
  
  enum DiffType {
    kFast,
    kGCRMA
  };

//  /// Types of intensities (log values)
//  enum IntensityMode {
//    kPmI = 0,
//    kMmI,    
//    kSensitivityCorrectedPm,
//    kSensitivityCorrectedMm,
//    kNsSubstractedFastDiff,
//    kNsSubstractedGcrmaDiff,
//    kNsSubstractedPm,
//    kNsSubstractedMm,
//    kIntensityModeCount // Contains the number of possible intensity types
//  };

  struct IntensityModeSelector {
    IntensityMode raw;
    IntensityMode corrected;
  };

  // Defines related modes for a PM or MM ananlysis
  const IntensityModeSelector intensityModeSelector[kProbeTypeCount] = {
    {kPmI, kSensitivityCorrectedPm},
    {kMmI, kSensitivityCorrectedMm}
  };


  /// Types of expression measures
  enum ExpressionMeasureProfileType{
    kEmProfilePm = 0, // Pm-only expression measure type
    kEmProfileMm,
    kEmProfileDiffFast,
    kEmProfileGcrmaFast,
    kExpressionMeasureProfileTypeCount
  };


  enum IntegrationType {
    kIntegrationSimple = 0, // Pm-only expression measure type
    kIntegrationDiff,
    kIntegrationDiffGcrma
  };

  struct ExpressionMeasureIntegrationParams {
    IntensityMode sourceIntensityMode;
    IntensityMode nsSubtractedIntensityMode;
    std::string gslParameterSet;
    // gsl_function integrationFunction;    
    IntegrationType integrationType;
    std::string normalizationParam;
  };

  const ExpressionMeasureIntegrationParams
  expressionMeasureIntegrationParams[kExpressionMeasureProfileTypeCount] = {
    {kPmI, kNsSubstractedPm, "gslPm", kIntegrationSimple, "normalizationSimple"}, // kEmProfilePm
    {kMmI, kNsSubstractedMm, "gslMm", kIntegrationSimple, "normalizationSimple"}, // kEmProfileMm
    {kNsSubstractedFastDiff, kNsSubstractedFastDiff, "gslSimpleDiff", kIntegrationDiff, "normalizationFastDiff"}, // kEmProfileDiffFast
    {kNsSubstractedGcrmaDiff, kNsSubstractedGcrmaDiff, "gslSimpleDiff", kIntegrationDiffGcrma,
     "normalizationFastDiff"} // kEmProfileGcrmaFast
  };
  
  
  struct DistributionParameters {
  	IntensityType mean;
  	IntensityType sigma;
  	IntensityType rightMargin;
  	IntensityType leftMargin;
   };
   
   struct BivariateDistributionParameters {
   	IntensityType meanPm;
   	//IntensityType meanMm;
   	IntensityType sigmaPm;
   	//IntensityType sigmaMm;
   	IntensityType rightMarginPm;
   	IntensityType leftMarginPm;
   	//IntensityType rightMarginMm;
   //IntensityType leftMarginMm;
   	IntensityType correlation;
   	IntensityType sigmaAverage;
   	IntensityType logB;
    };



  /**
   * Provides a function object to calculate how well the hookcurve
   * model fits a (digitized) graph. For certain parameters a,b,f, 
   * the function returns the sum of squares error.
   *
   */ 
  class ComputeHookcurveModelSumOfSquares : public std::unary_function<IntensityArray, IntensityType> 
  { 
  public:
    ComputeHookcurveModelSumOfSquares(IntensityMapping& intensityPlot);
    IntensityType operator() (IntensityArray& aAndF) const;
  //private:
   static IntensityType calculateRpmFromSumlogI(IntensityType a, IntensityType f, 
                                                IntensityType sumLogI);
   static IntensityType calculateRpmFromDeltalogI(IntensityType a, IntensityType f, 
                                                  IntensityType deltaLogI);

    IntensityMapping mIntensityPlot;
  };

  /**
   * Provides a function object to calculate how well two hookcurve
   * models fit a two (digitized) graphs. For certain parameters a1,b1,f1,
   * a2,b2,f2 the function returns the sum of two errors from 
   * ComputeHookcurveModelSumOfSquares.
   *
   */ 
  class ComputeSumOfTwoHookcurveModelFits : public std::unary_function<IntensityArray, IntensityType> 
  { 
  public:
    ComputeSumOfTwoHookcurveModelFits(IntensityMapping& intensityPlot1, IntensityMapping& intensityPlot2);
    IntensityType operator() (IntensityArray& params) const;
  private:
    ComputeHookcurveModelSumOfSquares mComputeFit1, mComputeFit2;
  };


  /**
   * @class HookcurveAnalyzer
   * @brief Aggregates processed data obtained from a single 
   * microarray experiment for hookcurve analysis.
   *
   * @todo Recheck if all member variables (as mAp,mBp, mNsChip) are required
   */
  class HookcurveAnalyzer
  {
  public:    
    HookcurveAnalyzer(const Chip& chip, // const
                      HookModel& hookModel/*,
                      IntensityArrayTable& correctedIntensities*/);
 /*   HookcurveAnalyzer(const ProbePtrVector& probes, 
                      const IntensitiesFunction& averageFunction);*/
    
    IntensityMappingPtr getAveragedHookcurvePlot(const ProbesetIntensitiesArray& sumLogIs, const ProbesetIntensitiesArray& deltaLogIs, const int movingAverageWindowSize = 5)  const;
    
    
    /// Calculates the specific portions of MM and PM and calculates specific sensitivity models
//    void calculateSpecificPortions(const IntensityPredicate& isLessThanIntersectionSumLogI,
//                                   const std::vector<ExpressionMeasureProfileType> profileTypes);
    
    
    ProbesetIntensitiesArrayPtr calculateSpecificPortionsMonovariate(const ProbesetIntensitiesArray&  rawIntensityArray, /*const ProbesetIntensitiesArray&  sensitivityCorrectedIntensityArray,*/ 
    		const SensitivityProfile& sensitivityProfile,
    		const DistributionParameters& distributionParameters, const PmMmProbeSequenceFunction& getProbeSequence = mem_fun_ref(&PmMmProbe::getSequence));
    
    ProbesetIntensitiesArrayPtr calculateSpecificPortionsBivariate(const ProbesetIntensitiesArray&  rawIntensityArrayPm, const ProbesetIntensitiesArray&  rawIntensityArrayMm,  
    		const SensitivityProfile& sensitivityProfilePm, const SensitivityProfile& sensitivityProfileMm,
    		const BivariateDistributionParameters& distributionParams,
    		double (*integrationFunction) (double, void*));
    

//    // Computes R for each probeset
//    std::vector<ProbesetExpressionPair> calculateProbesetR(const IntensityType nsThresholdSumLogI,
//                                                           const IntensityType a, const IntensityType f,
//                                                           const IntensityType offset) const;

     
    /// Checks if the chip saturates
    bool isSaturated(const IntensityMapping& graph, size_t intervals );
    
    
    /// Calculates the parameters for the compressionfactor function.
    void calculateCompressionFactorParameters(const ProbesetIntensitiesArray& probesetIntensitiesArray, const PmMmProbeSequenceFunction& getProbeSequence);
    
//    std::valarray<double> getSpecificityMarks(const ProbesetIntensitiesArray& sumLogIs,IntensityType threshold) const;

 
    /// Calculates the parameters a,b for the theoretic function
    boost::tuple<IntensityType, IntensityType, IntensityType> 
    estimateTheoreticFunctionParameters(const IntensityMapping& intensityPlot) const;

    static boost::tuple<IntensityType, IntensityType, IntensityType, IntensityType> 
    estimateTwoTheoreticFunctionParameters(IntensityMappingPtr hookcurvePlot1,
                                           IntensityMappingPtr hookcurvePlot2);

    
  private:
    
    IntensityArray integrationHelperMonovariate(const Probeset& currentProbeset, /*const size_t& probesetIndex, */
    		const IntensityArray& rawIntensities, /*const IntensityArray& sensitivityCorrectedIntensities,*/
    		const DistributionParameters& distribParams,
            const SensitivityProfile& sensitivityProfile,
            const PmMmProbeSequenceFunction& getProbeSequence);
    
    
    IntensityArray integrationHelperBivariate(
    		const Probeset& currentProbeset, 
            const IntensityArray& rawIntensitiesPm, const IntensityArray& rawIntensitiesMm,
    		const SensitivityProfile& sensitivityProfilePm, const SensitivityProfile& sensitivityProfileMm,
    		const BivariateDistributionParameters& distribParams,
    		const IntensityType& denominatorIntegral, 
    		double (*integrationFunction) (double, void*));

    
 
  private:  
    /// All probes used for this analysis
    const PmMmProbePtrVector& mProbes;

    /// Composition of probes to probesets with the
    /// constraint: Each probeset contains probes from a contigous segment of mProbes and for 
    /// each probe of \f$ probeset_i \f$ has a smaller index in mProbes than any other probe 
    /// from probeset \f$ i+1..n \f$.
    const ProbesetVector& mProbesets;

        
    /// Ap and Bp are the parameters of the compression function (as in Preibisch und page 48)
    double mAp;
    double mBp;

    
    
//     std::vector< boost::shared_ptr<ProbesetIntensitiesArray> > mIntensitiesArray;
    
    const Chip& mChip;
    // const
    HookModel& mHookModel; // @todo Make const when all setParam functions sourced out XSfrom HookcurveAnalyzer
    /*const IntensityType& mImax;*/
  };

}
#endif
