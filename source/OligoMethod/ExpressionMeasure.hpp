/**
 * @file ExpressionMeasure.hpp Provides methods for calculating the final Expression Measures
 * @author $Author: jan $ 
 * @author Jan Bruecker
 * @date $Date: 2008-09-08 15:04:08 +0200 (Mon, 08 Sep 2008) $
 */


#ifndef EXPRESSIONMEASURE_HPP_
#define EXPRESSIONMEASURE_HPP_

//#include "HookcurveAnalyzer.hpp"

//#include "Probe.hpp"
#include <boost/scoped_ptr.hpp>
#include "PmMmProbe.hpp"
#include "Probe.hpp"
#include "HookcurveAnalyzer.hpp"

namespace larrpack
{

  class ExpressionMeasure;

  /// Pointer to an hookcurve analyzer
  typedef boost::shared_ptr<ExpressionMeasure> ExpressionMeasurePtr;
  
  //  /// A pair of probeset identifier and its expression value
  //  typedef std::pair<std::string, IntensityType> ProbesetExpressionPair;

  enum ExpressionMeasureType {
    kEmPmConstantNs = 0, /// Em = Expression Measure
    kEmMmConstantNs,
    //   kEmDeltaConstantNs,
    //   kEmPmSingleIntegration,
    //   kEmMmSingleIntegration,
    //   kEmDeltaSingleIntegration,
    //   kEmDeltaGcrmaLikeIntegration,  /// Gcrma-like computation
    kEmPmConstantNsExp10, 
    kEmMmConstantNsExp10,
    kEmDeltaConstantNsExp10,
    kEmPmSingleIntegrationExp10,
    kEmMmSingleIntegrationExp10,
    kEmDeltaSingleIntegrationExp10,
    kEmDeltaGcrmaLikeIntegrationExp10,  /// Gcrma-like computation
    kEmPmSingleIntegrationExp10_Probe,
    kEmMmSingleIntegrationExp10_Probe,
    kEmDeltaSingleIntegrationExp10_Probe,
    kEmDeltaGcrmaLikeIntegrationExp10_Probe,  /// Gcrma-like computation
    kExpressionMeasureTypeCount // Contains the number of possible expression measures
  };

  /// Parameters to calculate various *specific profiles* 
  struct ExpressionMeasureProfileParams {
    IntensityMode intensityMode;
    ProfileType sProfile;
    std::string filenameSuffix;
    // std::string profileSuffix;
    PmMmProbeSequenceFunction getSequenceFunction;
  };

 
  const ExpressionMeasureProfileParams 
  expressionMeasureProfileParams[kExpressionMeasureProfileTypeCount] = {
    {kPmI, kSPm, "SPm", mem_fun_ref(&PmMmProbe::getSequence)}, // kEmProfilePm
    {kMmI, kSMm, "SMm", mem_fun_ref(&PmMmProbe::getSequenceMm)}, // kEmProfileMm
    {kNsSubstractedFastDiff, kSFastDiff, "SFastDiff", mem_fun_ref(&PmMmProbe::getSequence)}, // kEmProfileDiffFast
    {kNsSubstractedGcrmaDiff, kSGcrmaDiff, "SGCRMADiff", mem_fun_ref(&PmMmProbe::getSequence)} // kEmProfileGcrmaFast
  };
 
 

  class ExpressionMeasure
  {
  public:
    ExpressionMeasure(const Chip& chip, // const
                      HookModel& hookModel);
	
    //	virtual ~ExpressionMeasure();
	
    /// Initializes the needed Expression Measure
    void initializeExpressionMeasure(ExpressionMeasureType t);
    
    /// Calculates the expression values by integrating over the Ns distributions of the PM probes.
    std::vector<ProbesetExpressionPair> 
    calculateProbesetExpressionMeasure(const ProbesetIntensitiesArray& nsSubtractedIntensities, 
                                       const ProfileType profileType,
                                       const PmMmProbeSequenceFunction& getProbeSequence,
                                       UnaryIntensityFunction& invLogarithm,
                                       ExpressionMeasureType expMeasureType = kExpressionMeasureTypeCount);// If the default is active the reults won't be written to the array

	
    std::vector<ProbesetExpressionPair> 
    calculateProbeExpressionMeasures(const ProbesetIntensitiesArray& nsSubtractedIntensities,
                                     const ProfileType profileType,
                                     const PmMmProbeSequenceFunction& getProbeSequence,
                                     UnaryIntensityFunction& invLogarithm,
                                     ExpressionMeasureType expMeasureType = kExpressionMeasureTypeCount);
		
    /// Returns the expression measure of a given type for all probesets
    inline IntensityArray& getEms(ExpressionMeasureType t, size_t probesetIndex) const
    { return (*mExpressionMeasures[t])[probesetIndex]; }

    /// Returns all expression measures of a given type for all probesets
    inline const ProbesetIntensitiesArray& getEmsArray(ExpressionMeasureType t) const
    { return (*mExpressionMeasures[t]); }
    
    inline bool isEmPresent(ExpressionMeasureType t) const
    {return mExpressionMeasuresCalculated[t];}
	
	
  private:
    /// Contains expression measures for different methods
    std::vector<ProbesetIntensitiesArrayPtr> mExpressionMeasures;
    
    /// Basically the same vector as above, bur holds for each expression measure a bool if it is calculated.
    std::vector<bool> mExpressionMeasuresCalculated;
    
    const Chip& mChip;
    HookModel& mHookModel;
  };

}

#endif /*EXPRESSIONMEASURE_HPP_*/
