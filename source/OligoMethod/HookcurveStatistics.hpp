/**
 * @file HookcurveStatistics.hpp Defines functions to compute statistics on microarray and hookcurve data
 * @author $Author: mario $ 
 * @author Mario Fasold
 * @date $Date: 2008-04-04 17:54:25 +0200 (Fri, 04 Apr 2008) $
 */
#ifndef _HOOKCURVESTATISTICS_
#define _HOOKCURVESTATISTICS_

#include "OligoController.hpp"
#include "HookcurveAnalyzer.hpp"
#include "ExpressionMeasure.hpp"
#include "HookModel.hpp"


namespace larrpack {
  /// Maps statistics names to intensity values
  typedef std::map<std::string, IntensityType> IntensityTypeStatistics;

  /**
   * @class HookcurveStatistics
   * @brief Computes and prints statistics on microarray/hookcurve data.
   * 
   *
   *
   * @todo Convert HookcurveAnalyzerPtr to HookcurveAnalyzerConstPtr! 
   */
  class HookcurveStatistics
  {
  public:
    HookcurveStatistics(HookModelPtr hookModel, HookcurveAnalyzerPtr hookcurveModel,
    					ExpressionMeasurePtr exprMeasure,
                        const ProgramOptions options, const Chip& chip,
                        const  std::vector< boost::shared_ptr<ProbesetIntensitiesArray> >& intensitiesArray);

    IntensityTypeStatistics computeGeneralStatistics(IntensityType nsThreshold, IntensityType a, IntensityType f);
    IntensityTypeStatistics computeGeneralStatisticsPmOnly(IntensityType nsThreshold);
    
    IntensityTypeStatistics computeCurveStatistics(IntensityType nsThreshold, std::string suffix = "") const;
    IntensityTypeStatistics computeCurveStatisticsPmOnly(IntensityType nsThreshold, std::string suffix = "") const;
    
    /// Calculates the distribution for nonspecific Pm or Mms. If a filename is given, the result is plotted to the named data file.
//    boost::shared_ptr< std::vector<IntensityType> >//     std::vector<IntensityType>
//    calculateNonspecificDistribution(const IntensityPredicate filter, 
//                                     const IntensityMode intensityType,
//                                     std::string filename = std::string(""));


    void writeProbeDiagnosis(const IntensityType a, const IntensityType f, const std::string& fileName, const std::string& delimiter = "\t");
    void writeProbeDiagnosisPmOnly(const std::string& fileName, const std::string& delimiter = "\t");
    void writeProbesetDiagnosis(const IntensityType a, const IntensityType f, const std::string& fileName, const std::string& delimiter = "\t");
    void writeProbesetDiagnosisPmOnly(const std::string& fileName, const std::string& delimiter = "\t");


    void writeProbesetVariance(const std::string fileName, const std::string d = "\t") const;

    // @debug temporary function
    void printProbeDegration(std::string filename, IntensityMode pmType, IntensityMode mmType) const;

    // @debug temporary function
    void print3PrimeBiasAveragedHookcurve(std::string filename,IntensityMode pmType, IntensityMode mmType,
                                          IntensityType nsThreshold) const;

  private:
    IntensityMappingPtr getAveragedHookcurveInProbesetOrder(IntensityMode pmType, IntensityMode mmType, 
                                                            IntensityType defaultIntensity = 0.0) const;
  
    HookModelPtr mHookModel;

    // Reference to the model
    HookcurveAnalyzerPtr mHookcurveModel;
    
    //Reference to the Expression Measure Calculater
    ExpressionMeasurePtr mExpressionMeasureCalculater;
    
    // A reference to the intensities arrays 
    const  std::vector< boost::shared_ptr<ProbesetIntensitiesArray> >& mIntensitiesArray;

    /// The program options
    ProgramOptions options;
    
    const Chip& mChip;
  };

}

#endif

