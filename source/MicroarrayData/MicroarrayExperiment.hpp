/**
 * @file MicroarrayExperiment.hpp Contains MicroarrayExperiment class and defines common 
 *       datatypes.
 * @author $Author$ 
 * @date $Date$
 */
#ifndef _MICROARRAYEXPERIMENT_
#define _MICROARRAYEXPERIMENT_

#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>

#include "Probe.hpp"
#include "PmMmProbe.hpp"
#include "SnpProbe.hpp"
#include "Matrix.hpp"

namespace larrpack {
  /// Matrix of Intensity values
  typedef math::Matrix<IntensityType> IntensityMatrix;

  /// Pointer to an array of intensities
  /// @note: since the arrays get very, very large, we want to avoid
  /// copy operations (copy constructor of vector is used by default)
  /// @note: We use shared pointers - it's the easiest way of checking
  /// that all heap objects get deleted sometime!
  typedef boost::shared_ptr<IntensityMatrix> IntensityMatrixPtr;

  /// Array of bools
  typedef math::Matrix<bool> BoolArray;

  /// Pointer to an array of bools
  typedef boost::shared_ptr<BoolArray> BoolArrayPtr;

  /// Array of strings
  typedef math::Matrix<std::string> StringArray;

  /// Pointer to an array of strings
  typedef boost::shared_ptr<StringArray> StringArrayPtr;

  /**
   * @class MicroarrayExperiment
   * @brief Provides input functions for an microarray experiment.
   * 
   * There are various types of microarray data. This class provides 
   * methods to read all those formats into the internal, probe centered 
   * data type.
   */
  class MicroarrayExperiment
  {
  public:
    void readFromGenomicFile(const std::string& filename, const bool readProbesetComposition = false);
    void readFromProbesequenceFile(const std::string& filename, const IntensityMatrix& intensities);
    void readFromExonSequenceFile(const std::string& filename, const IntensityMatrix& intensities);
    void readPmOnlySequenceFile(const std::string& filename, const IntensityMatrix& intensities);
    void readFromSnpSequenceFile(const std::string& filename, const IntensityMatrix& intensities);
    void readFromBpmapfile(const std::string& filename, const IntensityMatrix& intensities, 
                           const StringArrayPtr& probesetIds = StringArrayPtr());
    void readFromPmapfile(const std::string& filename, const IntensityMatrix& intensities, 
                          const StringArrayPtr& probesetIds = StringArrayPtr());


    static void readArrayFromCelfile(const std::string& filename, IntensityMatrixPtr& intensityArrayPtr, 
                                     BoolArrayPtr& isMaskedArrayPtr);

    static void updateCelfileIntensities(const std::string sourceFilename, const std::string targetFilename, 
                                    IntensityMatrixPtr& intensityArrayPtr);

    static StringArrayPtr readProbesetIdsFromCdffile(const std::string& filename);

    ProbePtrVector getProbes();

  protected:
    /// All Probes on an microarray
    /// Datatype is an vector of smart pointers
    /// REASON: Probes are the fundamental object in the microarray analysis. Other objects which 
    /// perfom data analysis contain references to them. Smart pointers prevent running into
    /// memory issues. We use shared_ptr, which automaticly frees memory when no smart pointers 
    /// points to the object anymore.
    ProbePtrVector mProbes;
  };

}
#endif
