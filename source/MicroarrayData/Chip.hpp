/**
 * @file Chip.hpp Provides a data-structure to hold all chip-specific data
 * @author $Author: jan $ 
 * @date $Date: 2008-06-17 12:53:14 +0200 (Tue, 17 Jun 2008) $
 */
#ifndef CHIP_HPP_
#define CHIP_HPP_

#include <boost/function.hpp>
#include <boost/scoped_ptr.hpp>

#include "PmMmProbe.hpp"
#include "Probeset.hpp"

namespace larrpack
{
  class Chip; // Forward declaration

  /// Pointer to an hookcurve analyzer
  typedef boost::shared_ptr<Chip> ChipPtr;

  typedef boost::function<ProbesetVector (PmMmProbePtrVector)> ProbesetCompositionFunction;

  /// A function object to create the probeset composition
  typedef boost::function<bool (Probeset)> ProbesetPredicate;
  // Not all compilers like the declaration above. Alternatively, try:
  // typedef boost::function1<bool, Probeset> ProbesetPredicate;

  

  /// Defines different microarray chip types
  enum ChipType {
    kChiptypeGenechip = 0,
    kChiptypeTiling,
    kChiptypeSnp,
    kChiptypeExon,
    kChiptypeGenomicfile,
    kChiptypeExpressionArray
  };

  /**
   * @class Chip
   * @brief Aggregates microarray chip related data such as probes and probesets 
   */
  class Chip
  {
  public:
    // Ctors
    Chip(const PmMmProbePtrVector& probes, ProbesetVector probesets);
    Chip(const ProbePtrVector& probes);
    Chip(const ProbePtrVector& probes, ProbesetVector probesets);

    virtual ~Chip() {}
	
    /// Filters the probesets by a supplied filter function.
    void filterProbesets(ProbesetPredicate filterProbeset);
	
    // \name Getter methods
    //@{ 
    /// Get the probesets
    const  std::vector<Probeset>& getProbesets() const;
    /// Get the probes
    const PmMmProbePtrVector& getProbes() const;
    
    const size_t getProbesetCount() const;
    //@}
  private:
    PmMmProbePtrVector mProbes;
    ProbesetVector mProbesets;
  };

}

#endif /*CHIP_HPP_*/
