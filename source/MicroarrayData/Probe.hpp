/**
 * @file Probe.hpp Defines fundamental data types for a probe centric view.
 * @author $Author$ 
 * @author Mario Fasold
 * @date $Date$
 */
#ifndef _PROBE_
#define _PROBE_

#include <string>
#include <vector>
#include <valarray>
#include <boost/function.hpp>

#include <boost/shared_ptr.hpp>

namespace larrpack {
  /// Numerical type of the intensities
  typedef double IntensityType;
  
  /// A pair of intensity values
  typedef std::pair<IntensityType,IntensityType> IntensityPair;

  /// A set of intensity pairs
  typedef std::vector<IntensityPair> IntensityMapping;

  /// A pointer to a set of intensity pairs
  typedef boost::shared_ptr<IntensityMapping> IntensityMappingPtr;  

  /// An array of intensity values
  typedef std::valarray<IntensityType> IntensityArray;

  /// The length of a probe 
  const size_t kProbeSequenceLength = 25;

  /// A function object (handler) that computes on a list of intensity values
  typedef boost::function<IntensityType (IntensityArray)> IntensitiesFunction;
  // Not all compilers like the declaration above. Alternatively, try:
  // typedef boost::function1<IntensityType, IntensityArray > IntensitiesFunction;

  /// A function object predicate on intensity values
  typedef boost::function<bool (IntensityType)> IntensityPredicate;
  // Not all compilers like the declaration above. Alternatively, try:
  // typedef boost::function1<bool, IntensityType> IntensityPredicate;

  class Probe; // forward declaration to define the next types

  /// Pointer to a probe
  typedef boost::shared_ptr<Probe> ProbePtr;

  /// List of probes
  typedef std::vector<ProbePtr> ProbePtrVector;
  
  

  

  /**
   * @class Probe
   * @brief This class represents a single probe within a array experiment.
   * 
   * @todo Change pmPosition type to std::pair<size_t,size_t>.
   */
  class Probe 
  {
  public:
    Probe(const std::string& chromosome, const char strand, const std::string& sequence, 
          const int location, const IntensityType pm, 
          const std::pair<int,int>& pmPosition,
          const std::string& probesetId = "");

    std::string toString();


    // @name Getter and setter routines
    // @{
    signed int getLocation() const;
    IntensityType getPm() const;
    std::string getChromosome() const;
    std::string getSequence() const;
    std::string getProbesetId() const;
    char getStrand() const;
    
    // @debug:
    void setProbesetId(std::string id);
    
    std::pair<int,int> getPositionPm() const;
    
    // @}

  protected:
    /// The chromosome the probe is located on.(Generally, this could also be an arbitraty sequence
    /// identifier in case mutiple cell types are observed on one array.)
    std::string mChromosome;
  
    /// The strand
    char mStrand;

    /// The DNA/RNA sequence
    std::string mSequence;

    /// The genomic location on the chromosome
    signed int mLocation;

    /// The perfect match value
    IntensityType mPm;

    /// The position of the perfect match spot on the array
    std::pair<int,int> mPmPosition;

    /// Identifier of the corresponding probeset (in genechips)
    /// @note We need this to calculate the probeset composition
    ///       for genechips.
    std::string mProbesetId;
  };
}
#endif
