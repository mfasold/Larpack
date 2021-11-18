/**
 * @file SubsetSensitivityProfile.hpp Provides a class for sequence profile of a subset of nucleotides
 * @author $Author$ 
 * @author Mario Fasold
 * @date $Date: 2008-04-17 10:12:26 -0400 (Thu, 17 Apr 2008) $
 */
#ifndef _SUBSETSENSITIVITYPROFILE_
#define _SUBSETSENSITIVITYPROFILE_

#include <map>
#include <boost/shared_ptr.hpp>

#include "SensitivityProfile.hpp"

namespace larrpack {

  /**
   * @class SubsetSensitivityProfile
   * @brief Represents a sequence dependent intensity increment for some base tuples of each probe.
   *
   * In our model, each probe's measured intensity depends upon its sequence.
   * We denote that influence on the intensity "sequence dependent increment"
   * \f$ \delta y \f$ (See Binder, OligoMethodDescription).
   *
   * @todo Decide which functions better be private.
   */
  class SubsetSensitivityProfile : public SimpleSensitivityProfile
  {
  public:
    SubsetSensitivityProfile(std::vector<std::string> relevantBasetuples, size_t sequenceLength, 
                             unsigned short modelRank = 1);

    SubsetSensitivityProfile(const SubsetSensitivityProfile& s, unsigned short modelRank) : SimpleSensitivityProfile(s, modelRank), mRelevantBasetuples(s.mRelevantBasetuples)
    {}


    virtual IntensityType getSequenceIncrement(const std::string& sequence) const;

    virtual void computeProfile(const ProbesetSequencesArray& probesetSequencesArray,
                                const ProbesetIntensitiesArray& probesetIntensitiesArray,
                                const size_t probesetLimit);

    // @todo Rethink: does cloning with different rank make sense for these? 
    // Return a copy of this instance
    virtual SensitivityProfilePtr cloneWithNewRank(const size_t modelRank) const
    {
      //            SensitivityProfilePtr newProfile = SensitivityProfilePtr(new SubsetSensitivityProfile(*this, mModelRank));
      SensitivityProfilePtr newProfile = 
        SensitivityProfilePtr(new SubsetSensitivityProfile(mRelevantBasetuples, mSequenceLength,
                                                           mModelRank));
      boost::shared_static_cast<SubsetSensitivityProfile>(newProfile)->mProfile = mProfile;
      
      return newProfile; // @note Temporary shared pointers might be forbidden
      
//       SensitivityProfilePtr newProfile = SensitivityProfilePtr(new SubsetSensitivityProfile(*this, modelRank));
//       return newProfile; // @note Temporary shared pointers might be forbidden
    }

    // @note This is a hack to get the SequenceProfile functions to work easily here
    virtual size_t getProfileSize(unsigned short modelRank// , size_t sequenceLength, size_t nucleotideCount
                                  )
    { return mSequenceLength * mRelevantBasetuples.size(); }

    /// Returns the number of possible (Poly)nucleiotides within the subset
    virtual size_t getBasetupleCount() const 
    { return  mRelevantBasetuples.size();}

    // Returns the index of seq in the relevant base tuples if seq occurs, 
    // and an index > mRelevantBasetuples.size otherwise
    size_t sequenceToSubsetIndex(std::string seq) const
    {
      return find(mRelevantBasetuples.begin(), mRelevantBasetuples.end(), seq) - mRelevantBasetuples.begin();
    }

    virtual std::string indexToSequence(size_t baseIndex) const
    { return mRelevantBasetuples[baseIndex]; }

  private:
    /// A list of base tuples the profile shall be calculated for
    std::vector<std::string> mRelevantBasetuples;
  };

  /**
   * @class AdditivePofile
   * @brief A sequence profile that constitutes of multiple other profiles.
   *
   * This profile is combines two profiles. For example,  a simple profile 
   * (rank=1) might be used to estimate base effects while another profile 
   * computes GGG-based effects. To compute these profiles, the probessets
   * are split into two sets: (A) those containing probes that cover the 
   * relevant basetuples and (B) those that don't. The small first profile
   * is computed from set B, while the second profile is computed from set 
   * A.
   *
   * @todo Rethink handling of modelRank within class and constructor.
   */
  class AdditivePofile : public SensitivityProfile
  {
  public:
    AdditivePofile(std::vector<std::string> relevantBasetuples,
                   size_t sequenceLength, 
                   unsigned short modelRank);
    
    AdditivePofile(std::vector<std::string> relevantBasetuples, std::vector<SensitivityProfilePtr> profiles);

    virtual void computeProfile(const ProbesetSequencesArray& probesetSequencesArray,
                                const ProbesetIntensitiesArray& probesetIntensitiesArray,
                                const size_t probesetLimit);
  
    virtual IntensityType getSequenceIncrement(const std::string& sequence) const;
    virtual void exportToDatafile(const std::string filename) const;
    virtual void importFromDatafile(const std::string filename) {} // @todo Implement

    virtual SensitivityProfilePtr cloneWithNewRank(const size_t modelRank) const;

    virtual IntensityType getSigma(size_t baseIndex, size_t position) const;  
    virtual size_t getValidSequenceCount() const;
    virtual void zeroMiddlebase(); 

  private:
    std::vector<SensitivityProfilePtr> mProfiles;

    /// A list of base tuples 
    std::vector<std::string> mRelevantBasetuples;
  };
}
#endif
