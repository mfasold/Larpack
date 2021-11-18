/**
 * @file SensitivityProfile.hpp Contains classes to calculate sequence specific intensity increments
 * @author $Author$ 
 * @author Mario Fasold
 * @date $Date$
 */
#ifndef _SENSITIVITYPROFILE_
#define _SENSITIVITYPROFILE_

#include <valarray>

#include "Probe.hpp"
#include "PmMmProbe.hpp"
#include "Probeset.hpp"
// #include "boost/tuple/tuple.hpp"
#include "SequenceUtil.hpp"

namespace larrpack {
  class SensitivityProfile; // Forward declaration

  /// A pointer to a sensitivity profile
  typedef boost::shared_ptr<SensitivityProfile> SensitivityProfilePtr;

  class SimpleSensitivityProfile; // Forward declaration
  typedef boost::shared_ptr<SimpleSensitivityProfile> SimpleSensitivityProfilePtr;

  /// A function that returns a list of all intensity values of a probeset
  typedef boost::function<IntensityArray (Probeset)> ProbesetIntensitiesFunction;
  // Not all compilers like the declaration above. Alternatively, try:
  // typedef boost::function1<IntensityArray, Probeset> ProbesetIntensitiesFunction;

  /**
   * @class SensitivityProfile
   * @brief Interface for essential profile functions.
   * 
   */
  class SensitivityProfile
  {   
  public:
    virtual void computeProfile(const ProbesetSequencesArray& probsetSequencesArray,
                                const ProbesetIntensitiesArray& probesetIntensitiesArray,
                                const size_t probesetLimit) = 0;

    // Alternate interface for compute-profile (backwards compatibility)
    void computeProfile(const ProbesetVector& probesets,
                        const ProbesetIntensitiesArray& probesetIntensitiesArray,
                        const PmMmProbeSequenceFunction& getProbeSequence,
                        const size_t probesetLimit = 0);

    virtual SensitivityProfilePtr cloneWithNewRank(const size_t modelRank) const = 0; 

    virtual IntensityType getSequenceIncrement(const std::string& sequence) const = 0;

    virtual void exportToDatafile(const std::string filename) const = 0;
    virtual void importFromDatafile(const std::string filename) = 0;

    // @todo Potential design flaw: different profiles might use different indices
    // These functions are only used in HookcurveAnalyzer::calculateCompressionFactorParameters
    // and should be replaced/removed.
    virtual IntensityType getSigma(size_t baseIndex, size_t position) const = 0;
    virtual size_t getValidSequenceCount() const = 0;
    virtual void zeroMiddlebase() = 0;
    virtual ~SensitivityProfile() {};
  };

  /**
   * @class SimpleSensitivityProfile
   * @brief Calculates a sequence dependent intensity increment for each probe.
   *
   * In our model, each probe's measured intensity depends upon its sequence.
   * We denote that influence on the intensity "sequence dependent increment"
   * \f$ \delta y \f$ (See Binder, OligoMethodDescription).
   *
   * @note Encapsulation of the profile is reasonable since we eventually
   *       want to support different models (single nt, nearest-neighbor, ...)
   *
   * @todo Decide which functions better be private.
   */
  class SimpleSensitivityProfile : public SensitivityProfile
  {
  public:
    SimpleSensitivityProfile(size_t sequenceLength, unsigned short modelRank = 1);
    SimpleSensitivityProfile(size_t sequenceLength, unsigned short modelRank, size_t profileSize);
    SimpleSensitivityProfile(const SimpleSensitivityProfile& s, unsigned short modelRank);

    virtual IntensityType getSequenceIncrement(const std::string& sequence) const;

    virtual void computeProfile(const ProbesetSequencesArray& probsetSequencesArray,
                                const ProbesetIntensitiesArray& probesetIntensitiesArray,
                                const size_t probesetLimit);

    // void computePartlyInterpolatedProfile(const ProbesetVector& probesets, 
    //                                       const ProbesetIntensitiesArray& probesetIntensitiesArray,
    //                                       const PmMmProbeSequenceFunction& getProbeSequence,
    //                                       const std::vector<std::string> subsequences,
    //                                       const size_t sourceRank = 2,
    //                                       const size_t probesetLimit = 0,
    //                                       const size_t probesetSelectLimit = 5000);

    unsigned short getModelRank() const;
    static IntensityArray interpolateFromProfile(const SimpleSensitivityProfile& s, unsigned short modelRank);

    IntensityType getSigma(size_t baseIndex, size_t position) const;
    
    virtual void zeroMiddlebase();

    virtual void exportToDatafile(const std::string filename) const;
    virtual void importFromDatafile(const std::string filename);

    virtual IntensityType 
    computeSumOfSquares(const ProbesetVector& probesets,
                        const ProbesetIntensitiesArray& probesetIntensitiesArray,
                        const PmMmProbeSequenceFunction& getProbeSequence);

    static void
    exportConditionalSumOfSquares(const ProbesetVector& probesets,
                                  const ProbesetIntensitiesArray& probesetIntensitiesArray,
                                  const PmMmProbeSequenceFunction& getProbeSequence,
                                  const size_t baseCount, const std::string filename,
                                  const SensitivityProfilePtr profile);

    SensitivityProfilePtr cloneWithNewRank(const size_t modelRank) const;
    
    // Some small helper functions
    static size_t getProfileSize(unsigned short modelRank, size_t sequenceLength, size_t nucleotideCount);

    virtual size_t getProfileSize(unsigned short modelRank)
    { return getProfileSize(modelRank, mSequenceLength, sequenceutil::kNucleotideCount); }

    /// Returns the sequence length
    size_t getSequenceLength() const 
    { return mSequenceLength; } 

    /// Returns the number of possible (Poly)nucleiotides with the current model rank
    virtual size_t getBasetupleCount() const 
    { return sequenceutil::getBasetupleCount(mModelRank); }

    /// Returns number of base position taken into account with the current model rank
    virtual size_t getValidSequenceCount() const
    { return sequenceutil::getValidSequenceCount(mSequenceLength, mModelRank); }

    /// Gets index in profile for given base tuple and sequence index
    // @todo Use it everywhere in Sequenceprofile
    inline size_t getProfileIndex(std::string tuple, size_t sequenceIndex) const
    { return sequenceutil::sequenceToIndex(tuple)*getValidSequenceCount() + sequenceIndex; }

    static inline size_t getProfileIndex(size_t tupleIndex, size_t sequenceIndex, 
                                         size_t sequenceLength, size_t modelRank)
    { return tupleIndex*sequenceutil::getValidSequenceCount(sequenceLength, modelRank) + sequenceIndex; }

    static inline size_t getProfileIndex(std::string tuple, size_t sequenceIndex, 
                                         size_t sequenceLength, size_t modelRank)
    { return getProfileIndex(sequenceutil::sequenceToIndex(tuple),sequenceIndex, sequenceLength, modelRank); }

    /// Gets profile index to sequence 
    /// @todo replace all occurences
    virtual std::string indexToSequence(size_t baseIndex) const
    { return sequenceutil::indexToSequence(baseIndex, mModelRank); }

    IntensityType getSubTripleIncrement(const std::string& subTuple, size_t pos) const;
  protected: 
    IntensityArray getProfile() const  
    { return mProfile; }


  protected:
    /// The length of the probe sequences (->determines the size of the profile)
    const size_t mSequenceLength;

    /// Single-base, double-base or triple-base model
    unsigned short mModelRank;
  
    /// The profile 
    IntensityArray mProfile;
  private:   
    void probeIterationHelper(math::Matrix<double>& delta, std::valarray<double>& theta, 
                              const ProbesetSequencesArray& probsetSequencesArray, 
                              const size_t& probesetIndex, 
                              const ProbesetIntensitiesArray& probesetIntensitiesArray,
                              const size_t& modelNucleotides, const size_t& validSequencePositions, 
                              size_t profileSize);

  };


  /**
   * @class IncrementCalculator
   * @brief Calculates the sequence specific increment \f$ \delta S \f$.
   *
   * Given a sequence profile and the function to get the probe 
   * sequence, this class can be used as a function object in
   * algorithms.
   * 
   */
  class IncrementCalculator
  {
  public:
    /**
     * Constructor
     *
     * @param profile The sequence profile which is used
     * @param getProbeSequenceFunction The function to get a probe's sequence
     *
     */
    IncrementCalculator(const SensitivityProfilePtr& profile, 
                        const PmMmProbeSequenceFunction& getProbeSequenceFunction) 
      : mSensitivityProfile(profile), getProbeSequence(getProbeSequenceFunction) {}

    /**
     * Calculates the sequence specific increment for a probe.
     *
     * @return \f[ \delta S \f].
     */
    IntensityType operator() (const PmMmProbePtr& probe)
    {
      return mSensitivityProfile->getSequenceIncrement(getProbeSequence(*probe));
    }
  private:
    /// The sequence profile which is used
    const SensitivityProfilePtr mSensitivityProfile;

    /// The function to get a probe's sequence
    const PmMmProbeSequenceFunction getProbeSequence;
  };
}
#endif
