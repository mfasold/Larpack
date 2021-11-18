/**
 * @file BackgroundSubtraction.hpp Declares class and methodes for optical background correction (precorrection).
 * @author $Author$
 * @author Mario Fasold 
 * @date $Date$
 */
#ifndef _BACKGROUNDSUBTRACTION_
#define _BACKGROUNDSUBTRACTION_

#include <vector>

#include "Probe.hpp"
#include "OpticalBackgroundCorrection.hpp"
#include "MicroarrayExperiment.hpp"

namespace larrpack {

/**
 * @class BackgroundSubtraction
 * @brief Provides routines to calculate background correction on microarrays.
 * 
 * @todo Use typedefs IntensityArray etc.
 */
class BackgroundSubtraction : public OpticalBackgroundCorrection
{

private:
  /// Percentile of the values that is defined as noise.
  static const double kNoiseFraction;

  /**
   * @todo Find out why i cannot separate that brief description
   * to two lines
   */
  /// Small factor beeing added to assure that weight factor never be devided by zero.
  static const double kSmooth;

  /// The square root of the number of rectangular zones, the arry 
  /// is devided into (= number of rows/colums).
  int mZoneCount;

  /// Multiplicative factor applied to the background value to be subtracted from
  /// each intensity
  double mSubtractionRatio;

  /// Percentage of probes that is considered to be mainly optical backgrousnd
  double mOpticalBackgroundPercentage;

public:
  BackgroundSubtraction(const int zoneCount, const double subtractionRatio = 1.0, 
                        const double opticalBackgroundPercentage = 0.02);

  // Performs background correction
  virtual void computeCorrection(IntensityMatrix& intensities, 
                                 BoolArray& isMasked);
private:
  static int getZone(int x, int y, int maxX, int maxY, int zoneCount);
  static std::pair<double,double> getZoneCenter(int zoneX, int zoneY, int maxX, int maxY, int zoneCount);
  static IntensityType calculateWeightedProbeBackground(const std::pair<int, int>& spotPosition, 
                                                         const std::vector<IntensityType>& zoneValues, 
                                                         const std::vector< std::pair<double, double> >& zoneCenters);
  IntensityType calculateCorrectedIntensity(const IntensityType intensity, 
                                             const std::pair<int, int>& spotPosition, 
                                             const std::vector<IntensityType>& zoneBackgroundMeans, 
                                             const std::vector<IntensityType>& zoneBackgroundStandardDeviations,
                                             const std::vector< std::pair<double, double> >& zoneCenters);
};

}

#endif
