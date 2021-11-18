/**
 * @file HookModel.hpp Declares the model used to describe microarray
 * intensities 
 * @author $Author: jan $ 
 * @date $Date: 2008-06-17 12:53:14+0200 (Tue, 17 Jun 2008) $
 */
#ifndef _HOOKMODEL_
#define _HOOKMODEL_

#include <map>
//#include <functional>
#include "SensitivityProfile.hpp"
#include "MathUtil.hpp"


namespace larrpack {

  class HookModel; // Forward declaration

  /// Pointer to an hookcurve model
  typedef boost::shared_ptr<HookModel> HookModelPtr;

  /// Different ways to calculate the hookcurve
  enum HookcurveCalculationType {
    kRawIntensitiesOnly,
    kNonspecificProfilesOnly,
    kAllProfiles
  };
  
  /// Type of a probe: either perfect-match or mismatch
  enum ProbeType {
    kPm,
    kMm,
    //kDiff // Not a real probetype, but treated as one when we calculate Expression Measures by difference
    kProbeTypeCount
  };

  /// Types of sensitivity profiles
  enum ProfileType {
    kNsPm = 0, // Profile type for non-specific (absent target) region, perfect-match values
    kNsMm,
    kSPm, 
    kSMm,
    kSFastDiff,
    kSGcrmaDiff,
    kProfileTypeCount // Contains the number of possible profile types
  };

  struct RelatedCalculationParams {
    ProbeType probeType;
    ProfileType nsProfile;
    ProfileType sProfile;
    std::string nsChip;
    PmMmProbeSequenceFunction getSequence;
  };

  // Defines all the related parameters used for a PM or MM ananlysis
  const RelatedCalculationParams probeTypeParams[kProbeTypeCount] = {
    {kPm, kNsPm, kSPm, "NsPmChip", mem_fun_ref(&PmMmProbe::getSequence)},
    {kMm, kNsMm, kSMm, "NsMmChip", mem_fun_ref(&PmMmProbe::getSequenceMm)}
  };


  /**
   *
   *
   *
   */
  class HookModel 
  {
  public:
    //     // \name The main estimator steps of the algorithm
    //     //@{ 
    //     IntensityPair detectNsDomain(const IntensityMappingPtr& graph) = 0; /*Hook dependent method*/ 

    //     /// Calculates the ns profiles on the ns domain
    //     void computeNonspecificSequenceProfiles(const Chip& chip, IntensityArray& intensities,
    //                                             size_t probesetLimit) = 0;

    //     boost::tuple<IntensityType, IntensityType>  computeImaxAndA(const IntensityMappingPtr& graph) = 0;

    //     /// After destuaration and subtracting Ns contant the specific profiles can be calculated
    //     void computeSpecificSequenceProfiles(const Chip& chip, IntensityArray& intensities) = 0; 
    //     //@}
    HookModel();
	  
    virtual ~HookModel() {}
  
    // \name Model-solving functions for the specific, sequence corrected contant
    //@{ 
    /// Gives back the unsuturated, unlogged intensity value.
    IntensityType getDesaturatedIntensity(const IntensityType logSaturatedIntensity, 
                                          const IntensityType cutOffMargin = 0.90) const;

    
    inline IntensityArray getDesaturatedIntensities(IntensityArray logSaturatedIntensities, 
                                              const IntensityType cutOffMargin = 0.90)
    { 
      IntensityArray desaturatedIntensities = logSaturatedIntensities;
    	for (size_t probeIndex = 0; probeIndex < logSaturatedIntensities.size(); ++probeIndex) {
    		desaturatedIntensities[probeIndex] = getDesaturatedIntensity(logSaturatedIntensities[probeIndex]);
    	}
    	//logSaturatedIntensities.apply( std::mem_fun_ref(&HookModel::getDesaturatedIntensity) );  
       return desaturatedIntensities; 
    }

    /// Computes the fraction of non-specific binding accoring to hook paper for large beta: 10^(-  (sigma - nsThreshold))
    inline IntensityType getNsFraction(IntensityType sumLogI) const 
    {
      return MathUtil::trimToRange( exp10(-(sumLogI - mThreshold.first )) , 0.0, 1.0);
    }
    
    // @note Rewrite as correctIntensities
    IntensityArray calculateIncrements(const Probeset& probeset, const IntensityArray& intensities,
                                       ProbeType probeType, IntensityType nonspecificFraction)  const;

    void desaturateIntensities(/*Parameters are pointers to the intensity arrays and Imax*/);
    void subtractNsBackground();
    
    void subtractSpecificSequenceProfiles();
    //@}
    
    // \name Parameter estimator functions (maybe source out)
    //@{ 
    
    IntensityType estimateImax(const std::vector< std::valarray<double> >& probesetIntensitiesArray, 
                               const size_t topX = 10);

  //protected:
    /// Returns the compression factor that is multiplied to the specific sensitivity profiles.
    IntensityType getCompressionFactor(IntensityType intensity) const;
  public:
  
    //@}

    
//    /// Gives back the insaturated, unlogged intensity values of a valarray
//    std::valarray<double> getDesaturatedIntensities(IntensityArray logSaturatedIntensities) const;
//    
//    std::vector< valarray<double > > getDesaturatedIntensities(const std::vector< valarray<double > >& rawProbesetIntensitiesarray) const;
    
    // \name Getter and setter methods
    //@{
    void setImax(IntensityType iMax);   
    const IntensityType& getImax() const;
    
    /// Get the coordinates of the recent intersection point.
    IntensityPair getNsThreshold() const;
    /// Set the new intersectionpoint coordinates.
    void setNsThreshold(IntensityPair sThreshold);
    
    std::valarray<double> markPresentProbes(const ProbesetIntensitiesArray& intensityArray, IntensityType threshold) const;

    void setParam(const std::string id, const IntensityType value) { mParams[id] = value; }
    IntensityType getParam(const std::string id) const { return (*mParams.find(id)).second; }
      //return mParams[id]; } //@todo make const
    
    /// Returns a sequence profile of a particular type
    inline boost::shared_ptr<SensitivityProfile> getProfile(const ProfileType profileType) const
    { return mSensitivityProfiles[profileType]; }

    /// Sets the sequence profile of a particular type
    void setProfile(const ProfileType profileType, const boost::shared_ptr<SensitivityProfile>& profile) 
    { mSensitivityProfiles[profileType] = profile; }
    
    //@}
  protected:  	  
    void initializeSensitivityProfiles();
	    
    /// The maximum intensity value (as log value)
    IntensityType mImax;

    /// Current kink point coordinates
    IntensityPair mThreshold;

    /// The mean non-specific intensity (=cross-hybridization fraction)
    IntensityType mNsIntensity;

    // Stores minor model parameters
    std::map<std::string, IntensityType> mParams;

    /// Contains all profiles to calculate sequence specific increments (for NsPm,NsMm,SPm,...)
    /// @note We decided to use a vector since std::map does not provide constant-time
    /// access using operator []. As a drawback, the vector contains a pointer for
    /// every possible profile type, though not every needn't be in use.
    std::vector<SensitivityProfilePtr> mSensitivityProfiles;
    // std::map<ProfileType, SensitivityProfilePtr> mSensitivityProfiles;
  };


  class PmMmHookModel : private HookModel {

  };
    
}
#endif

