/**
 * @file BackgroundSubtraction.cpp Defines methodes for background subtraction.
 * @author $Author$
 * @author Mario Fasold
 * @date $Date$
 */ 

#include <vector>
#include <algorithm>
#include <cmath>

// ______DEBUG
#include <iostream>
// ----------

#include "MathUtil.hpp"
#include "DataCollector.hpp"
#include "BackgroundSubtraction.hpp"
  
using namespace std;
using namespace larrpack;


const double BackgroundSubtraction::kNoiseFraction = 0.5; 
const double BackgroundSubtraction::kSmooth = 0.00001;


/**
 * Constructor.
 *
 * @param zoneCount The square root of the number of rectangular zones, the arry 
 *        is devided into (= number of rows/colums)
 * @param subtractionRatio Multiplicative factor applied to the background value 
 *        to be subtracted from each intensity.
 * @param opticalBackgroundPercentage Percentage of probes that is considered to 
 *        be mainly optical background
 */
BackgroundSubtraction::BackgroundSubtraction(const int zoneCount, const double subtractionRatio, 
                                             const double opticalBackgroundPercentage) 
  : mZoneCount(zoneCount), mSubtractionRatio(subtractionRatio), mOpticalBackgroundPercentage(opticalBackgroundPercentage)
{} 

/**
 * Calculates the zone which a spot on the microarray is located in.
 * It trys to devide in the area in \f$zoneCount^2\f$ equal rectangles. 
 * If that is not 
 * possible, the remainder elements are added to the last zone. To be accurate,
 * for all zones (i,j) hold
 *   \f[start_i=\lfloor i*W \rfloor,start_j=\lfloor i*H \rfloor \f]
 * and
 *   \f[end_i=\begin{cases}
 *     maxX& \text{if $i=zoneCount$},\\
 *     \lfloor (i+1)*W \rfloor -1& \text{otherwise}.
 *   \end{cases}\f]
 * and
 *   \f[end_j=\begin{cases}
 *     maxY& \text{if $j=zoneCount$},\\
 *     \lfloor (j+1)*H \rfloor -1& \text{otherwise}.
 *   \end{cases}\f]
 *
 * where \f$W=(maxX+1)/zoneCount\f$ is the width of zone and
 * \f$H=(maxY+1)/zoneCount\f$ is the height of zone.
 *
 *
 *
 * @param x X-coordinate of the spot on the chip.
 * @param y Y-coordinate of the spot on the chip.
 * @param maxX Maximum x-coordinate on the chip. (The smallest is 0).
 * @param maxY Maximum y-coordinate on the chip. (The smallest is 0).
 * @param zoneCount The square root of the number of rectangular zones, the arry 
 *        is devided into (= number of rows/colums).
 *
 * @return Index of that zone.
 *
 */
inline int BackgroundSubtraction::getZone(int x, int y, int maxX, int maxY, int zoneCount) 
{
  int elementsPerZoneX = (maxX + 1) / zoneCount; // minX = 0
  int elementsPerZoneY = (maxY + 1) / zoneCount; // minY = 0
  return min(y/elementsPerZoneY, zoneCount - 1)*zoneCount + min(x/elementsPerZoneX, zoneCount - 1);
}

/**
 * Calculates center coordinates of a zone. 
 *
 * @see BackgroundSubtraction::getZone
 *
 * @param zoneX The zone's column (\f$0 <= ZoneX < zoneCount\f$).
 * @param zoneY The zone's row. (\f$0 <= ZoneY < zoneCount\f$).
 * @param maxX Maximum x-coordinate on the chip. (The smallest is 0).
 * @param maxY Maximum y-coordinate on the chip. (The smallest is 0).
 * @param zoneCount The square root of the number of rectangular zones, the arry 
 *        is devided into (= number of rows/colums).
 *
 * @return Center coordinates of a zone.
 */
pair<double,double> BackgroundSubtraction::getZoneCenter(int zoneX, int zoneY, int maxX, int maxY, int zoneCount) 
{
  int elementsPerZoneX = (maxX + 1) / zoneCount; // minX = 0
  int elementsPerZoneY = (maxY + 1) / zoneCount; // minY = 0

  // Calculate last point within each zone, with
  // special case for the last zone
  int endZoneX = (zoneX + 1)*elementsPerZoneX - 1;
  if (zoneX == zoneCount - 1) {
    endZoneX = maxX;
  }

  int endZoneY = (zoneY + 1)*elementsPerZoneY - 1;
  if (zoneY == zoneCount - 1) {
    endZoneY = maxY;
  }

  // return (startpoint + endpoint) / 2
  return pair<double,double>(double(zoneX*elementsPerZoneX + endZoneX) / 2.0,
                             double(zoneY*elementsPerZoneY + endZoneY) / 2.0);
}

/**
 * Calculates background value \f$ b(x,y) \f$ as weighted mean of the zone centers. See Preibisch, 
 * 2006 S.38 for further explanation of the formula.
 *
 * @param spotPosition The x/y coordinates of the spot on the chip.
 * @param zoneValues A List containing the backdround values (either mean \f$ bZ_k \f$ or standard
 *        deviation \f$ nZ_k \f$) for each zone.
 * @param zoneCenters A List containing center coordinates for each zone.
 *
 * @return Background value b(x,y).
 *
 * @todo Test for zoneValue.size = zoneCenters.size
 * @todo Make function inline?
 * @todo Check if double precision is needed.
 */
IntensityType BackgroundSubtraction::calculateWeightedProbeBackground(const std::pair<int, int>& spotPosition, 
                                                                   const std::vector<IntensityType>& zoneValues, 
                                                                   const std::vector< std::pair<double, double> >& zoneCenters)
{
  double weightTotal = 0;
  double weightedBackgroundTotal = 0;
  for (size_t zone = 0; zone < zoneValues.size(); ++zone) {
    // see Preibisch, S.38
    double weight_factor = 1/(sqrt(pow((double) spotPosition.first - zoneCenters[zone].first, 2) + 
                                   pow((double) spotPosition.second - zoneCenters[zone].second, 2)) + kSmooth);
    weightedBackgroundTotal += weight_factor*(double) zoneValues[zone];
    weightTotal += weight_factor;
  }
  return (IntensityType) weightedBackgroundTotal/weightTotal;
}


/**
 * Calculates a single corrected background intensity according to the affymetrix
 * formulas.
 *
 * @param intensity The uncorrected intensity.
 * @param spotPosition The x/y coordinates of the spot on the chip.
 * @param zoneBackgroundMeans A List containing the backdround mean \f$ bZ_k \f$ for each zone.
 * @param zoneBackgroundStandardDeviations A List containing the backdround standard
 *        deviation \f$ nZ_k \f$ for each zone.
 * @param zoneCenters A List containing center coordinates for each zone.
 *
 * @return Corrected intensity.
 */
IntensityType BackgroundSubtraction::calculateCorrectedIntensity(const IntensityType intensity, const std::pair<int, int>& spotPosition, 
                                   const std::vector<IntensityType>& zoneBackgroundMeans, 
                                   const std::vector<IntensityType>& zoneBackgroundStandardDeviations,
                                   const std::vector< std::pair<double, double> >& zoneCenters)
{
  // Calculate background to the coordinate
  IntensityType backgroundValue = calculateWeightedProbeBackground(spotPosition,zoneBackgroundMeans, zoneCenters);
  // Calculate noise to the coordinate
  IntensityType backgroundNoise = calculateWeightedProbeBackground(spotPosition, zoneBackgroundStandardDeviations, zoneCenters);

  // Calculate corrected intensity using those values
  IntensityType correctedIntensity = max(intensity, 0.5);
  return max(correctedIntensity - mSubtractionRatio*backgroundValue, kNoiseFraction * backgroundNoise);   
}


/**
 * Performs a background subtraction as suggested by Affymetrix Inc. Refer to 
 * Affymetrix technical whitepaper for exact description of the algorithm.
 *
 * @param intensities 2-dimensional array of intensity values
 * @param isMasked 2-dimensional array, stating for each intensity value whether or not it should be used
 * 
 * @see http://www.affymetrix.com/support/technical/whitepapers/sadd_whitepaper.pdf , Page 4
 * 
 * @todo Check row-major order of all functions. 
 * @todo Check Matrix rowNum and colNum for non quadradic matrices
 */
void BackgroundSubtraction::computeCorrection(IntensityMatrix& intensities, 
                                              BoolArray& isMasked
                                              )
{
  // Check basic preconditions
  assert(intensities.rowNum() == isMasked.rowNum());
  assert(intensities.columnNum() == isMasked.columnNum());
  assert(mZoneCount > 0);
  
  // Calculate maximal x and y locations a genechip (minX and minY are 0 by default)
  int maxX = intensities.rowNum() - 1;
  int maxY = intensities.columnNum() - 1;
  
  // Initialize array of size mZoneCount*mZoneCount of variable-sized lists! The 
  // list shall contain ALL PM/MM values (no distinction between those) within that zone.
  vector< vector<IntensityType> > zoneValues(mZoneCount*mZoneCount); 

  // [Note: another possibilty would be to save the current lowest 2% of each zone. this would
  // be less memory demanding but not faster (than n*log n) and much more complicated]

  // For each spot on the array
  for (int x = 0; x <= maxX; ++x) {
    for (int y = 0; y <= maxY; ++y) {
      if (!isMasked[x][y]) {
        //  Determine the zone the value lies in
        //  and it to that zone's list
        zoneValues[getZone(x, y, maxX, maxY, mZoneCount)].push_back(intensities[x][y]);
      }
    }
  }

  // Sort the arrays of values for each zone
  for (vector< vector<IntensityType> >::iterator zi = zoneValues.begin(); zi != zoneValues.end(); ++zi) {
    sort((*zi).begin(), (*zi).end()); // @note partial_sort should be enough
  }

  //   Initialize arrays for background mean and -variance (of size mZoneCount*mZoneCount)
  vector<IntensityType> zoneBackgroundMeans;
  vector<IntensityType> zoneBackgroundStandardDeviations;

  //   For each zone
  for (vector< vector<IntensityType> >::iterator zi = zoneValues.begin(); zi != zoneValues.end(); ++zi) {
    //     Calculate the mean intensity of the lowest 2% probes
    int background_fraction = int(mOpticalBackgroundPercentage*double(zi->size())) + 1 ;    
    // @note The first k values are copied into a new vector to pass it to the calculateMean 
    // function -> this is not efficient
    zoneBackgroundMeans.push_back(MathUtil::calculateMean(vector<double>(zi->begin(), 
                                                                         zi->begin() + background_fraction)));
    
    //     Calculate the variance intensity of the lowest 2% probes
    zoneBackgroundStandardDeviations.push_back(MathUtil::calculateStandardDeviation(vector<double>(zi->begin(), 
                                                                                                   zi->begin() + background_fraction)));
  }  

  //   Initialize array of pairs.
  //   And fill that array with center position x/y of each zone
  //   for better performance in following calculations
  vector< pair<double,double> > zoneCenters;
  for (int y = 0; y < mZoneCount; ++y) {
      for (int x = 0; x < mZoneCount; ++x) {
        zoneCenters.push_back(getZoneCenter(x,y, maxX, maxY, mZoneCount));
      }
  }

  //   For each spot on the array
  for (int x = 0; x <= maxX; ++x) {
    for (int y = 0; y <= maxY; ++y) {      
      if (!isMasked[x][y]) {
        //     Update value with corrected intensity
        intensities[x][y] = calculateCorrectedIntensity(intensities[x][y],
                                                        pair<int, int> (x, y),
                                                        zoneBackgroundMeans, zoneBackgroundStandardDeviations, zoneCenters);
       
//     int zone = getZone((*pi)->getPositionPm().first, (*pi)->getPositionPm().second, maxX, maxY, mZoneCount);
//     cout << (*pi)->getPm() << "\t(" << (*pi)->getPositionPm().first << "," << (*pi)->getPositionPm().second << ")\t";
//     cout << zone <<  "\t" << pmCorrected << endl;
//     cout << (*pi)->getPm() << "\t" << (*pi)->getMm() << "\t";
//     cout << pmCorrected << "\t" << mmCorrected << endl;
      }
    }
  }

  // Log some background subtraction parameters
  double bgMeansAverage = accumulate(zoneBackgroundMeans.begin(), zoneBackgroundMeans.end(), 0.0)
    / zoneBackgroundMeans.size();
  DataCollector::instance().insert("backgroundSubtractionBgMean", log10(bgMeansAverage));
}


