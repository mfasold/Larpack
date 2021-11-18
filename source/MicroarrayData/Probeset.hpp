/**
 * @file Probeset.hpp Defines a datastructure for aggregated probes (PmMmProbes).
 * @author $Author$ 
 * @author Mario Fasold
 * @date $Date$
 */
#ifndef _PROBESET_
#define _PROBESET_

#include <vector>
#include <valarray>

#include "PmMmProbe.hpp"
#include "SnpProbe.hpp"
#include "Matrix.hpp"

namespace larrpack {
  /// A list of valarrays
  typedef std::vector<IntensityArray> ProbesetIntensitiesArray;

  /// A pointer to a list of valarrays
  typedef boost::shared_ptr<ProbesetIntensitiesArray> ProbesetIntensitiesArrayPtr;

  // A list of sequences 
  // @todo Find appropriate place for this
  typedef std::vector< std::vector<std::string> > ProbesetSequencesArray;
  typedef boost::shared_ptr<ProbesetSequencesArray> ProbesetSequencesArrayPtr;

  /// Function to robestely compute average intensity (used in many classes)
  extern IntensitiesFunction computeAverageIntensity;

  class Probeset; // forward declaration to define the next types

  typedef std::vector<Probeset> ProbesetVector;

  
  typedef std::vector< boost::shared_ptr<ProbesetIntensitiesArray> > IntensityArrayTable;
  
  typedef boost::shared_ptr< std::vector< boost::shared_ptr<ProbesetIntensitiesArray> > > IntensityArrayTablePtr;

  /// Types of intensities (log values)
  enum IntensityMode {
    kPmI = 0,
    kMmI,    
    kPseudoMmI,
    kSensitivityCorrectedPm,
    kSensitivityCorrectedMm,
    kSensitivityCorrectedPseudoMm,
    kNsSubstractedFastDiff,
    kNsSubstractedGcrmaDiff,
    kNsSubstractedPm,
    kNsSubstractedMm,
    kWashedPm,
    kIntensityModeCount // Contains the number of possible intensity types
  };
  
  /**
   * @class Probeset
   * @brief Aggregates individual PmMmProbes and provides methods and collects data
   * relevant for the further analysis of the set.
   *
   * Currently, probesets are regarded a a set of PmMmProbes, since the analysis
   * is based on the mismatch values. The fundamental purpose of the probeset datastructure
   * is not only aggegating the probes but also providing analysis specific functions
   * as getAverageSumLogI.
   *
   * @note Probesets should only be contigous slices of the probe list.
   * 
   */
  class Probeset 
  {
  public:
    // Probeset();
    Probeset(const PmMmProbePtrVector& probes);
    Probeset(PmMmProbePtrVector::const_iterator probesetBegin,
             PmMmProbePtrVector::const_iterator probesetEnd,
             std::string probesetId = "");

    /// Gets the number of elements within the set
    size_t getSize() const;

    // Getter Functions
    PmMmProbe getProbe(const size_t probeIndex) const;
    PmMmProbePtr getProbePtr(const size_t probeIndex) const;

    // Get iterators
    // @todo Move to cpp if reasonable
    PmMmProbePtrVector::const_iterator beginProbe() const { return mProbes.begin(); }
    PmMmProbePtrVector::const_iterator endProbe() const { return mProbes.end(); }
  
    // Returns the probeset id
    std::string getProbesetId() const;
    
    // @debug
    void setProbesetId(std::string id);

    IntensityArray getProbeIntensities(const IntensityMode pmOrMm) const;
    //std::vector<std::string> getProbeSequences(const IntensityMode pmOrMm) const;
    std::vector<std::string> getProbeSequences(const PmMmProbeSequenceFunction& getProbeSequence) const;



    static ProbesetIntensitiesArrayPtr initializeIntensitiesArray(const ProbesetVector& probesets);
    static ProbesetIntensitiesArrayPtr initializeIntensitiesArray(const ProbesetVector& probesets, 
                                                                  const IntensityMode pmOrMm);

  private:
    /// Pointers to probes in this set
    PmMmProbePtrVector mProbes;

    std::string mProbesetId;
  };

  /**
   * @class ProbeCountEquals
   * @brief Defines an unary operator for probesets such that it returns true
   * if and only if the set contains exactly x probes.
   * 
   * @todo Source out to cpp and document.
   */
  class ProbeCountEquals : public std::unary_function<Probeset, bool> 
  {
    size_t mProbeCount;
  public:
    explicit ProbeCountEquals(const size_t probeCount) : mProbeCount(probeCount) {}
    bool operator() (const Probeset& probeset) const 
    { return probeset.getSize() == mProbeCount;  }    
  };



  // ------------------------------------------------------------------------------------------
  // @todo Source the following out to SnpProbeset.hpp
  // ------------------------------------------------------------------------------------------
  class SnpProbeset; // forward declaration to define the next types

  typedef std::vector<SnpProbeset> SnpProbesetVector;
  /**
   * @class SnpProbeset
   * @brief Aggregates individual SnpProbes and provides methods and collects data
   * relevant for the further analysis of the set.
   *
   * Basic idea is to seperate between Allele (Perfect-match) and Cross (Mismatch)
   * Probes. 
   *
   * @note Probesets should only be contigous slices of the probe list.
   * @todo Releate with Probeset Class. Problem: different types of probe vector. 
   *       Maybe create templated probeset!
   * 
   * @todo SnpProbesets should be a data-structur that actually contains two probesets: Cross and Allele
   */
  class SnpProbeset 
  {
  public:
    // Probeset();
    // SnpProbeset(const SnpProbePtrVector& probes);
    SnpProbeset(SnpProbePtrVector::const_iterator probesetBegin,
                SnpProbePtrVector::const_iterator probesetEnd,
                std::string probesetId,
                const char call);

    /// Gets the number of elements within the set
    size_t getSize() const;

    // Getter Functions
    // PmMmProbe getProbe(const size_t probeIndex) const;
    SnpProbePtr getProbePtr(const size_t probeIndex) const;

    // Get iterators
    SnpProbePtrVector::const_iterator beginProbe() const { return mProbes.begin(); }
    SnpProbePtrVector::const_iterator endProbe() const { return mProbes.end(); }
  
    // Returns the probeset id
    std::string getProbesetId() const;

    char getCall() const { return mCall; }
    
    // IntensityArray getProbeIntensities(const IntensityMode pmOrMm) const;
    // //std::vector<std::string> getProbeSequences(const IntensityMode pmOrMm) const;
    //std::vector<std::string> getProbeSequences(const PmMmProbeSequenceFunction& getProbeSequence) const;
    std::vector<std::string> getProbeSequences() const;

    static ProbesetIntensitiesArrayPtr initializeIntensitiesArray(const SnpProbesetVector& probesets);
    // static ProbesetIntensitiesArrayPtr initializeIntensitiesArray(const ProbesetVector& probesets, 
    //                                                               const IntensityMode pmOrMm);

    // --------------------------------------- REMOVE!  ---------------------------------------------
    // The following functions swaps Allele and Cross Probes
    // @todo REMOVE! This is very bad programming practice, function used only temporary for convenience. All 
    //       Probeset functions should be availble separately for Cross and Allele.
    void swapCross() { 
      SnpProbePtrVector tmp = mProbes; // Copies vector of smart pointers
      mProbes = mProbesCross; 
      mProbesCross = tmp; 
    }
    // --------------------------------------- REMOVE!  ---------------------------------------------
  private:
    /// Pointers to probes in this set
    SnpProbePtrVector mProbes, mProbesCross;

    std::string mProbesetId;
    
    /// Genpotype Call
    char mCall;  
  };

}
#endif
