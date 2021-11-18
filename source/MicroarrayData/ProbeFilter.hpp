/**
 * @file ProbeFilter.hpp Defines filters for (PmMm-)Probes.
 * @author $Author$ 
 * @author Mario Fasold
 * @date $Date$
 */
#ifndef _PROBEFILTER_
#define _PROBEFILTER_

#include <functional>
#include "Probe.hpp"
#include "SnpProbe.hpp"

namespace larrpack {

  /**
   * @class ProbeintensityIsLessThan
   * @brief Defines an unary operator that filters probes with an intensity
   * less than a given theshold.
   * 
   * @note Dependencies: This class tests for PmMmProbes as well.
   */
  class ProbeintensityIsLessThan : public std::unary_function<ProbePtr, bool> 
  {
  private:
    IntensityType threshold;
  public:
    explicit ProbeintensityIsLessThan(const IntensityType& intensity) : threshold(intensity) { }
    bool operator() (const ProbePtr& probe) const; 
      
  };
  
  /**
   * @class IsNotInProbesetFilter
   * @brief Defines a unary operator that returns true iff a probe's 
   * probeset_id is in a given list of probesets ids (validProbesetids)
   * 
   */
  class IsInProbesetFilter : public std::unary_function<ProbePtr, bool>
  {
    private:
      const std::vector<std::string> validProbesetIds;
    public:
      explicit IsInProbesetFilter(const std::vector<std::string> validIds) : validProbesetIds(validIds) {
      }
      bool operator() (const ProbePtr& probe) const;
  };

  /**
   * @class ProbeintensityIsLessThan
   * @brief Defines an unary operator that filters probes the sequence of which
   * contains a substring at a given position.
   * 
   * @note Dependencies: This class tests for PmMmProbes as well.
   */
  class ProbesequenceHasSubstring : public std::unary_function<ProbePtr, bool> 
  {
  private:
    std::string mSubstring;
    size_t mIndex;
  public:
    explicit ProbesequenceHasSubstring(const std::string substring, 
                                       const size_t index = std::string::npos) 
      : mSubstring(substring), mIndex(index) { }
    bool operator() (const ProbePtr& probe) const; 
      
  };

  /**
   * @class ProbeChromosomeEquals
   * @brief Defines an unary operator for probes that 
   * tests for a certain chromosome.
   * 
   * @note Dependencies: This class tests for PmMmProbes as well.
   */
  class ProbeChromosomeEquals : public std::unary_function<ProbePtr, std::string> 
  {
  private:
    std::string mChromosome;
  public:
    explicit ProbeChromosomeEquals(const std::string chromosome) : mChromosome(chromosome) { }
    bool operator() (const ProbePtr& probe) const;       
  };


  /**
   * @class ProbesetIdIsLexicallySmaller
   * @brief A binary operator to compare the probeset id of two probes.
   * 
   */
  class ProbesetIdIsSmaller : public std::binary_function<ProbePtr,ProbePtr, bool> 
  {
  public:
    explicit ProbesetIdIsSmaller() { }
    bool operator() (const ProbePtr& probeLhs, const ProbePtr& probeRhs) const;       
  };

  /**
   * @class ProbesetIdIsLexicallySmaller
   * @brief A binary operator to compare the probeset id, allele and 
   *        interrogation of two SNP probes.
   * 
   */
  class ProbesetIdAndAlleleAndInterrogationIsSmaller : public std::binary_function<SnpProbePtr, SnpProbePtr, bool> 
  {
  public:
    explicit ProbesetIdAndAlleleAndInterrogationIsSmaller() { }
    bool operator() (const SnpProbePtr& probeLhs, const SnpProbePtr& probeRhs) const;       
  };

}
#endif
