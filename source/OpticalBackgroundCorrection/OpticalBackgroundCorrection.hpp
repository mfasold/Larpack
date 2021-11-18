/**
 * @file OpticalBackgroundCorrection.hpp Defines an interface for optical background correction methods.
 * @author $Author$ 
 * @author Mario Fasold
 * @date $Date$
 */
#ifndef _OPTICALBACKGROUNDCORRECION_
#define _OPTICALBACKGROUNDCORRECION_

#include <vector>
#include "Probe.hpp"
#include "MicroarrayExperiment.hpp"

namespace larrpack {

/**
 * @class OpticalBackgroundCorrection
 * @brief An abstract interface used by various methods correcting optical 
 * background of microarrays.
 * 
 */
class OpticalBackgroundCorrection 
{
public:
  // Virtual classes must have a virtual destructor
  virtual ~OpticalBackgroundCorrection() {}

  /**
   * Abstract function to compute the background correction.
   *
   * @param intensities 2-dimensional array of intensity values
   * @param isMasked 2-dimensional array, stating for each intensity value whether or not it should be used
   */
  virtual void computeCorrection(IntensityMatrix& intensities, 
                                 BoolArray& isMasked) = 0;

};

}
#endif
