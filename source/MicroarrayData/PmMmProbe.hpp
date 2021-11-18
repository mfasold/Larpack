/**
 * @file PmMmProbe.hpp Defines data type for a pm-mm-probe centric view.
 * @author $Author$ 
 * @author Mario Fasold
 * @date $Date$
 */
#ifndef _PMMMPROBE_
#define _PMMMPROBE_

#include <string>
#include <vector>
#include <valarray>

#include <boost/shared_ptr.hpp>

#include "SequenceUtil.hpp"

#include "Probe.hpp"

namespace larrpack {  
  class PmMmProbe; // Forward definition to declare the next types

  /// Pointer to a PmMmProbe
  typedef boost::shared_ptr<PmMmProbe> PmMmProbePtr;

  /// List of probes
  typedef std::vector<PmMmProbePtr> PmMmProbePtrVector;

  /// A function that returns the seqence of a PmMmProbe
  typedef boost::function<std::string (PmMmProbe)> PmMmProbeSequenceFunction;
  // Not all compilers like the declaration above. Alternatively, try:
  // typedef boost::function1<std::string, PmMmProbe> PmMmProbeSequenceFunction;


  /**
   * @class PmMmProbe
   * @brief Represents a probe pair containing perfect match
   *        as well as mismatch intensity values.
   * 
   * PmPmProbe is used for affymetrix chips which yield an intensity 
   * value for the (actual) perfect match probe and another intensity
   * for a so-called mismatch probe. That is, one base in the middle
   * of the probe sequence is not complementary.
   */
  class PmMmProbe : public Probe
  {
  public:
    PmMmProbe(const std::string& chromosome, const char strand, const std::string& sequence, 
          const int location, const IntensityType pm, const IntensityType mm, 
              const std::pair<int,int>& pmPosition, const std::pair<int,int>& mmPosition,
              const std::string& probesetId = "");

    // @name Getter and setter routines
    // @{
    IntensityType getMm() const;
    std::pair<int,int> getPositionMm() const;
    std::string getSequenceMm() const;
    // @}

    IntensityType getDeltaLogI() const;
    IntensityType getSumLogI() const;
    
    double getSensitivityPm(std::valarray<double>& sensitivityArray) const;
    double getSensitivityMm(std::valarray<double>& sensitivityArray) const;

  private:
    /// The mismatch value
    IntensityType mMm;

    /// The position of the mismatch spot on the array
    std::pair<int,int> mMmPosition;
  };

  PmMmProbePtrVector extractPmMmProbesFromProbelist(const ProbePtrVector& probes);
}
#endif
