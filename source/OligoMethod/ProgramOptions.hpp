/**
 * @file ProgramOptions.hpp Reads options from the command-line
 * @author $Author: mario $ 
 * @author Mario Fasold
 * @date $Date: 2007-04-20 16:55:02 +0200 (Fri, 20 Apr 2007) $
 */
#ifndef _PROGRAMOPTIONS_
#define _PROGRAMOPTIONS_

#include <boost/scoped_ptr.hpp>
#include <Probeset.hpp>
#include <ExpressionMeasure.hpp>
#include <set>


#include "DetectKinkPoint.hpp"
#include "OpticalBackgroundCorrection.hpp"
#include "BackgroundSubtraction.hpp"




namespace larrpack {
  /// All major options needed by the controller
  struct ProgramOptions {
    // Contrains all (boolean) option flags
    // @todo Insert all boolean options
    std::set<std::string> flags;

    DetectKinkPointPtr detectKinkPoint;
    boost::shared_ptr<OpticalBackgroundCorrection> backgroundCorrection;
    ChipType chipType;
    size_t profileModelRank;
    std::string intensityFilename;
    std::string chipDiscriptionFilename;
    std::string genotypeCallFilename;
    size_t hookcurveMovingAverageWindowSize;

    std::vector<std::string> gstackTuples;

    /// Maximum probes to use for the computation of sequence profiles
    size_t probesetLimitForProfileComputation;

    /// Shall the theoretic hookcurve be fitted to compute f (saturation)
    bool fitTheoreticHookcurve; 

    /// Parameters selecting kink point computation
    bool computeCorrectedKinkPoint;
    bool computeUncorrectedKinkPoint;
    IntensityType correctedKinkPoint;
    IntensityType uncorrectedKinkPoint;

    // Skaling factor for I_max (used for testing, typically equals 1)
    double ImaxScalingFactor;

    /// Paramteres for optical background correction
//     size_t backgroundZones;
//     double backgroundRatio;
    
    /// Parameter for expression measure selection
    std::string expressionMeasureIdentifier;
    std::vector<ExpressionMeasureProfileType> requiredProfileTypes; // Profiles that need to be computed

    /// Probes with 1/2(PM+MM) distant to the breakpoint by the cutOff
    /// are considered to calculate the specific sensitivity profiles.
    IntensityType specificityCutoff;

    /// Parameters for probeset shuffling
    size_t optimalProbesPerSet;
    size_t maxProbesPerSet;
    int maxProbeDistance;
    size_t shufflingMovingAverageWindowSize;
    double specificityShufflingTreshold;  

    /// Parameters used for debug/testing
    bool isDebugMode;
    std::string maskfile;
    size_t probesetSizeFilter;
    std::string universalParam;
    bool readProfilesFromFile;
  };  

  /// Sets default options and parses them from command line stings
  ProgramOptions parseCommandlineOptions(int argc, char *argv[]);
}

#endif

