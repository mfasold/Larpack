/**
 * @file SpikedinProbe.hpp Contains a class for spiked-in probes.
 * 
 * @author $Author$
 * @date $Date$
 */
#ifndef _SPIKEDINPROBE__
#define _SPIKEDINPROBE_

#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "Probe.hpp"

namespace larrpack {  
  class SpikedinProbe; // forward definition to declare the next types

  /// Pointer to a spikedin probe
  typedef boost::shared_ptr<SpikedinProbe> SpikedinProbePtr;

  /// List of probe pointers
  typedef std::vector<SpikedinProbePtr> SpikedinProbePtrVector;

  /**
   * @class SpikedinProbe
   * @brief A spiked in probe is a special probe with a fixed and
   *        well defined concentration.
   */
  class SpikedinProbe : public Probe
  {
  public:
  protected:
  };
}
#endif
