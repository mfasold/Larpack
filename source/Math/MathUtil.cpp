/**
 * @file MathUtil.cpp Defines various mathematical functions.
 * @author $Author$ 
 * @author Mario Fasold
 * @date $Date$
 *
 */
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cassert>

#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include <gsl/gsl_statistics.h>

#include "MathUtil.hpp"
#include "Matrix.hpp"
#include "MathVector.hpp"
#include "ValarrayUtil.hpp"
#include "StlUtil.hpp"


using namespace std;
using namespace larrpack;
using namespace MathUtil;

/**
 *     
 * Calculates the moving average \f$ m_i \f$ for all values \f$ v_i, i=0..n-1 \f$ of an input list. 
 * The result vector has the same size as the input vector. It contains the usual moving average 
 * for those values we can compute it for. For the others - the first and the last values of the 
 * list - constant values are used. Precisely, the gerneralized moving average is  
 *   \f[m_i=\begin{cases}
 *     \frac{1}{s} \sum_{k=i-\lfloor s/2 \rfloor}^{i+\lfloor s/2 \rfloor} v_k \
 *       & \text{if $\lfloor s/2 \rfloor <= i < n - \lfloor s/2 \rfloor$} \\
 *     m_{\lfloor s/2 \rfloor} & \text{if $i < \lfloor s/2 \rfloor$} \\
 *     m_{n - \lfloor s/2 \rfloor - 1} & \text{if $n - \lfloor s/2 \rfloor <= i$}
 *    \end{cases}\f]
 * where \f$ s = movingAverageWindowSize \f$.
 * 
 * @note Why not use a normal moving average? It is not defined for the first and the last indexes.
 *       By using constant values for those indexes, we avoid monkeying with index corrections.
 * @note Return by value (using the copy constructor) needs more time
 * @note Use uneven movingAverageWindowSizes for precise average calculation.
 *
 * @todo Try compiling that slice thing with other compilers/systems.
 * 
 * @param values List of values.
 * @param movingAverageWindowSize Number of values to average over.
 * @return Generalized moving average of values.
 */
std::valarray<FloatType> MathUtil::calculateGeneralizedMovingAverage(const std::valarray<FloatType>& values, 
                                                                  size_t movingAverageWindowSize) 
{
  assert(values.size() > movingAverageWindowSize);

  // Initialize new (temporal) vector
  valarray<FloatType> movingAverages(0.0, values.size());
  
  // Create a view on those indexes, the moving average can be calculated for
  int firstAveragedIndex = movingAverageWindowSize/2 + 1 - 1; // integer division
  int validMovinAveragesCount = values.size() - movingAverageWindowSize + 1; // TEST
  slice_array<FloatType> validMovingAverages = movingAverages[slice(firstAveragedIndex, 
                                                                 validMovinAveragesCount, 
                                                                 1)];
                                                                 
  // Calculate the sum of values in the moving average window.
  // For all offsets in the moving average window
  for (size_t offset = 0; offset < movingAverageWindowSize; ++offset) {
    // Add this offset's values to the sum of all offsets
    // @note With my current knowledge, it seems not possible to create a const slice_array 
    // of a const valarray. Therefore, i have to cast the const away. However, i never
    // change the valarray!
    const slice_array<FloatType> currentValuesView = const_cast<std::valarray<FloatType>& >(values)[slice(offset, validMovinAveragesCount, 1)];
    validMovingAverages += currentValuesView;
  }

  // Devide all values by window size
  movingAverages /= (FloatType)movingAverageWindowSize;

  // Set constant values for first and last indices
  movingAverages[slice(0, firstAveragedIndex, 1)] = movingAverages[firstAveragedIndex];
  int lastAveragedIndex = firstAveragedIndex + validMovinAveragesCount - 1;
  movingAverages[slice(lastAveragedIndex + 1, 
                       movingAverages.size() - lastAveragedIndex - 1, 1)] = movingAverages[lastAveragedIndex];

  return movingAverages;
}

FloatPairListPtr
MathUtil::calculateMovingAverage(const FloatPairList& mapping, size_t movingAverageWindowSize)
{
  // Calculate first Point which the moving average can be calculated for
  size_t firstAveragedDatapoint = movingAverageWindowSize/2 + 1 - 1; // integer division

  // Check if averaging list large enough 
  if (mapping.size() <= movingAverageWindowSize) {
    // @note Is returning an empty vector good? Maybe return a single average
    // (using an accumulate and a select2nd functor)
    FloatPairListPtr emptyPlot(new FloatPairList(0));
    return emptyPlot;
  }
 
  // Initialize averagedPlot with zeros
  FloatPairListPtr averagedPlot(new FloatPairList(mapping.size() - movingAverageWindowSize + 1));
  
  // Calculate the sum of values in the moving average window:
  // For all frames in moving average window
  for (size_t averageFrame = 0; averageFrame < movingAverageWindowSize; ++averageFrame) {
    // Add this frame's values to the sum of all frames
    for (size_t i = 0; i < averagedPlot->size(); ++i) {
      (*averagedPlot)[i].second += mapping[i + averageFrame].second;
    }
  }

  for (size_t i = 0; i < averagedPlot->size(); ++i) {
    // Devide all values by window size
    (*averagedPlot)[i].second /= movingAverageWindowSize;

    // and copy first coordinate from the intensities
    (*averagedPlot)[i].first = mapping[i + firstAveragedDatapoint].first;
  }    

  return averagedPlot;
}



/**
 * Calculates the mean of a list of values.
 *
 * @note Uses vector's array values directly, which might not be clean.
 * 
 * @param values List of values.
 * @return Mean of the values.
 */
FloatType MathUtil::calculateMean(const vector<FloatType>& values) 
{
  // use gsl function to calculate mean
  return gsl_stats_mean(&values[0], 1, values.size());
}


/**
 * Calculates the sample standard deviation of a list of values.
 *
 * @note Uses vector's array values directly, which might not be clean.
 * 
 * @param values List of values.
 * @return  Standard deviation of the values.
 */
FloatType MathUtil::calculateStandardDeviation(const vector<FloatType>& values) 
{
  // use gsl function to calculate SAMPLE standard deviation ( \sqrt(1/{N-1} * \sum_1^N ( x - \mu))  )
  return gsl_stats_sd(&values[0], 1, values.size());
  //return gsl_stats_absdev(&values[0], 1, values.size());
}


/**
 * Calculates the mean of a list of values.
 *
 * @todo Solve MSC oddity (note it works for vectors).
 *
 * @param values List of values.
 * @return Mean of the values.
 */
FloatType MathUtil::calculateMean(const valarray<FloatType>& values) 
{
  // use gsl function to calculate mean
#ifdef _MSC_VER
  const FloatType* tmp = (const double*)&values;	// VC++, istead of &values[0]
  return gsl_stats_mean(tmp, 1, values.size());
#else
  return gsl_stats_mean(&values[0], 1, values.size());
#endif           // _MSC_VER 
}

/**
 * Calculates the sample standard deviation of a list of values.
 *
 * @todo Solve MSC oddity (note it works for vectors).
 *
 * @param values List of values.
 * @return  Standard deviation of the values.
 */
FloatType MathUtil::calculateStandardDeviation(const valarray<FloatType>& values) 
{
  // use gsl function to calculate SAMPLE standard deviation ( \sqrt(1/{N-1} * \sum_1^N ( x - \mu))  )
#ifdef _MSC_VER
  const FloatType* tmp = (const double*)&values;	// VC++, istead of &values[0]
  return gsl_stats_sd(tmp, 1, values.size());
#else
  return gsl_stats_sd(&values[0], 1, values.size());
#endif           // _MSC_VER 

}


FloatType MathUtil::calculateMedian(const valarray<FloatType>& values)
{
	valarray<FloatType> sortedValues(values);
	size_t valueCount = sortedValues.size();
	sort(&sortedValues[0], &sortedValues[valueCount]);
	FloatType median;
	if (valueCount % 2 == 1) {
		median = sortedValues[(valueCount - 1)/ 2];
	    }
	else {
		median = (sortedValues[(valueCount/2) - 1] + sortedValues[(valueCount/2)]) / 2.0;
	}
	return median;
//	valarray<FloatType> absDeviations = abs(values - median);
//	sort(&absDeviations[0], &absDeviations[valueCount]);
//	FloatType mad;
//	if (valueCount % 2 == 1) {
//		mad = absDeviations[(valueCount - 1)/ 2];
//	    }
//	else {
//		mad = (absDeviations[(valueCount/2) - 1] + absDeviations[(valueCount/2)]) / 2.0;
//	}
//	return mad;
}


/**
 * Calculates the MAD (median absolute deviation) from a valarray of floats
 * 
 * @todo This could be done more nicely using the median class of AveragingProcedure.hpp ?
 * 
 * @param values valarray of floats
 * 
 * @return the MAD of the valarrays data
 * 
 */
FloatType MathUtil::calculateMAD(const valarray<FloatType>& values)
{
//	valarray<FloatType> sortedValues(values);
//	size_t valueCount = sortedValues.size();
//	sort(&sortedValues[0], &sortedValues[valueCount]);
//	FloatType median;
//	if (valueCount % 2 == 1) {
//		median = sortedValues[(valueCount - 1)/ 2];
//	    }
//	else {
//		median = (sortedValues[(valueCount/2) - 1] + sortedValues[(valueCount/2)]) / 2.0;
//	}
	size_t valueCount = values.size();
	FloatType median = calculateMedian(values);
	valarray<FloatType> absDeviations = abs(values - median);
	sort(&absDeviations[0], &absDeviations[valueCount]);
	FloatType mad;
	if (valueCount % 2 == 1) {
		mad = absDeviations[(valueCount - 1)/ 2];
	    }
	else {
		mad = (absDeviations[(valueCount/2) - 1] + absDeviations[(valueCount/2)]) / 2.0;
	}
	return mad;
}


/**
 * Calculates the Covariance
 * 
 * @note Uses the same vector to array transformation as before.
 * @param valuesX List of values of random variable X
 * @param valuesY List of values of random variable Y
 */
FloatType MathUtil::calculateCovariance(const vector<FloatType>& valuesX, const vector<FloatType>& valuesY)
{
  assert( valuesX.size() == valuesY.size() );
  return gsl_stats_covariance(&valuesX[0], 1, &valuesY[0], 1, valuesX.size());
}

/**
 * Calculates the Covariance when the means are known
 * 
 * @note Uses the same vector to array transformation as before.
 * @param valuesX List of values of random variable X
 * @param valuesY List of values of random variable Y
 * @param muX mean of random variable X
 * @param muY mean of random variable X
 */
FloatType MathUtil::calculateCovariance(const vector<FloatType>& valuesX, const vector<FloatType>& valuesY, FloatType muX, FloatType muY)
{
  assert( valuesX.size() == valuesY.size() );
  return gsl_stats_covariance_m(&valuesX[0], 1, &valuesY[0], 1, valuesX.size(), muX, muY);
}


/**
 * Calculates the Correlation Coefficient (Pearson's)
 * 
 * @note Uses the same vector to array transformation as before.
 * @param valuesX List of values of random variable X
 * @param valuesY List of values of random variable Y
 */
FloatType MathUtil::calculateCorrelation(const vector<FloatType>& valuesX, const vector<FloatType>& valuesY)
{
  assert( valuesX.size() == valuesY.size() );
  FloatType cov  = gsl_stats_covariance(&valuesX[0], 1, &valuesY[0], 1, valuesX.size());
  FloatType stdX = gsl_stats_sd(&valuesX[0], 1, valuesX.size());
  FloatType stdY = gsl_stats_sd(&valuesY[0], 1, valuesY.size());
  return cov / (stdX * stdY);
}

/**
 * Calculates the Correlation Coefficient (Pearson's) when the standard deviations and means are known
 * 
 * @note Uses the same vector to array transformation as before.
 * @param valuesX List of values of random variable X
 * @param valuesY List of values of random variable Y
 * @param muX mean of random variable X
 * @param muY mean of random variable X
 * @param stdX standard deviation of random variable X
 * @param stdY standard deviation of random variable Y
 */
FloatType MathUtil::calculateCorrelation(const vector<FloatType>& valuesX, const vector<FloatType>& valuesY, 
                                      FloatType muX, FloatType muY, FloatType stdX, FloatType stdY)
{
  assert( valuesX.size() == valuesY.size() );
  FloatType cov  = gsl_stats_covariance_m(&valuesX[0], 1, &valuesY[0], 1, valuesX.size(), muX, muY);
  //FloatType stdX = gsl_stats_sd(&valuesX[0], 1, valuesX.size());
  //FloatType stdY = gsl_stats_sd(&valuesY[0], 1, valuesY.size());
  return cov / (stdX * stdY);
}





void print(const math::Matrix<FloatType>& m) {
//    for(size_t i = 0; i < m.dim1(); ++i)  {
//       for(size_t j = 0; j < m.dim2(); ++j)
   for(size_t i = 0; i < m.rowNum(); ++i)  {
      for(size_t j = 0; j < m.columnNum(); ++j)
        cout << '\t' << m(i,j);
      cout << endl;
   }
   cout << endl;
}

/// Define the default gLog parameter
FloatType MathUtil::gGLogParameter = 4.0;

/**
 * Calculates the generalized logarithmic function of a given value 
 * using log10,
 * \f$ log_{10}(\frac {z + \sqrt{z^2 + c^2}}{2}) \f$.
 *
 * @note For \f$ c = 2 \f$ gLog for negative values behaves 
 * like \f$ -log(-z) \f$.
 *
 * @warning The global gGLogParameter has to be initialized properly.
 * 
 * @param z Value the glog is caluclated for.
 * @param c gLog parameter. 
 * @return The glog10 Value
 */
FloatType MathUtil::gLog10(FloatType z) {
  return log10((z + sqrt(z*z + MathUtil::gGLogParameter*MathUtil::gGLogParameter))/2.0);
}

/**
 * Calculates the inverse glog function value of a given value (uses log10)
 *
 * @warning The global gGLogParameter has to be initialized properly.
 * 
 * @param log10Value Value the glog value.
 * @return The inverse glog
 */
FloatType MathUtil::gExp10(FloatType log10Value){
  return exp10(log10Value) - MathUtil::gGLogParameter*MathUtil::gGLogParameter * exp10(-log10Value)/4.0;
}


// FloatType MathUtil::exp10(FloatType log10Value){
//   return exp(log10Value) / exp(10);
// }
// FloatType exp10(FloatType log10Value);


/**
 * Computes 
 * \f$ e^{\displaystyle -\frac{1}{2} (\frac{x - \mu}{2*\mu / \sigma})^2} \f$.
 * Note that \f$ \mu \neq 0 \f$ must hold.
 * 
 *
 * @param x X point.
 * @param mu \f$ \mu \f$ - the center of the gauss curve.
 * @param sigma \f$ \sigma \f$ - the bigger this value, the more spiky the curve 
 *        is. For \f$ \sigma = 1 \f$, the curve matches a parabola \f$ y=-x^2 \f$.
 *
 * @note Multiplying by \f$ \displaystyle \frac{1}{\sqrt{2*\pi}} \f$  is not 
 * necessary, since weighting eliminates this constant factor.
 */
FloatType MathUtil::getGaussianWeight(const FloatType x, const FloatType mu, const FloatType sigma)
{
  assert( mu != 0.0 );
  return exp(-0.5 * pow((x - mu) / ((mu*2.0)/sigma), 2.0));
}


/**
 * Robust version of the digitizeCurve. First tests if in each of the given intervals at least one data point can be found, 
 * otherwise it resets the number of intervals to a size where enough data points are available. If intervalCound is 0
 * it always
 * @param intervalCount the number of desired intervals. if it is set to 0 an "optimized" number is calculated
 * i.e. it calculates the biggest number of equally spaced intervals where at least one data point is in each interval.
 */
std::vector<FloatPair> MathUtil::digitizeCurveRobust(const std::vector<FloatPair>& graph, size_t intervalCount, FloatType sigma)
{ 
  //FloatType intervalSize;
  // Calculate biggest distance between two points
  FloatType maxDistance = 0.0; // Set to initial value
  for (size_t i = 0; i < graph.size()-1; ++i) {
    FloatType d = abs(graph[i+1].first - graph[i].first); // Calculates distance of two neighbored values
    if (d > maxDistance) {
      maxDistance = d; 
    }
  }
  size_t intervalCountTmp = (size_t) ((graph.back().first - graph.front().first) / maxDistance); // Calculates how many intervals of length maxDistance are needed
  intervalCountTmp -= 2; // Use little less intervals due to the previous int casting...
  if (intervalCount > intervalCountTmp) {
    if (intervalCountTmp < 20) intervalCountTmp = 20; // We should at least have something like 15 data points.
    cout << "set digitize-intervals from " << intervalCount << flush;
    intervalCount = intervalCountTmp;
    cout << " to " << intervalCount << " (maximal number that is valid for this graph. 15 is minimum)" << endl;
    }
  return digitizeCurve(graph, intervalCount, sigma);
  }

/**
 * Computes a new graph with only sampleCount points. Each of the new points
 * is the gaussian weighted average of a interval of points in the original
 * graph.
 *
 * @pre The graph must be sorted w.r.t. to the first x-axis.
 *
 * @note The value of an interval with no graph points will be NaN.
 * 
 * @param graph Sorted set of points in the original graph.
 * @param intervalCount Number of points the new graph shall contain. 
 * @param sigma Describes behaviour of the gaussian weight (@see MathUtil::getGaussianWeight)
 *
 */
std::vector<FloatPair> MathUtil::digitizeCurve(const std::vector<FloatPair>& graph, size_t intervalCount, FloatType sigma)
{
  FloatType intervalSize;
  
  std::vector<FloatPair> digitizedCurve;

  // Compute the size of an interval
  intervalSize = (graph.back().first - graph.front().first) / intervalCount;

  // Initialize graph point iterator
  std::vector<FloatPair>::const_iterator currentPoint = graph.begin();
  
  // For each interval
  for (size_t currentInterval = 0; currentInterval < intervalCount; ++currentInterval) {
    // Update current interval start and end
    FloatType intervalBegin = graph.front().first + currentInterval * intervalSize;
    FloatType intervalEnd = graph.front().first +  (currentInterval + 1) * intervalSize;

    // Compute central x-point
    FloatType centralX = (intervalBegin + intervalEnd) /2;    

    // Add gaussian weighted Y-point and total weight for 
    // each point in that interval
    FloatType weightedYTotal = 0.0;
    FloatType weightTotal = 0.0;
    size_t pointsInInterval = 0;
    while ((currentPoint != graph.end()) && 
           (currentPoint->first <= intervalEnd)) {
      weightedYTotal += getGaussianWeight(currentPoint->first, centralX, sigma) * currentPoint->second;
      weightTotal += getGaussianWeight(currentPoint->first, centralX, sigma);
      ++currentPoint;
      ++pointsInInterval;
    }

    if (pointsInInterval >= 2) { // Only add a digitized point if at least 5 points were in the interval (to unrobust otherwise)
    // Add new point to result graph
    digitizedCurve.push_back(FloatPair(centralX, weightedYTotal/weightTotal));
    }
  }

  return digitizedCurve;
}



/**
 * A second digitize function, that digitizes rather the raw data than the already smoothed hookcurve
 * Apperently this version shows good results when used for the theoretical hookcurve fit.
 * 
 *
 * @pre The graph must be sorted w.r.t. to the first x-axis.
 *
 * @param graph Sorted set of points in the original graph.
 * @param stepSize Distance between two adjacent digitized points on the X axis 
 * @param intervalSize Size of the area that is integrated into the center point
 *
 */
FloatPairList MathUtil::digitizeCurveNew(const std::vector<FloatPair>& graph, FloatType stepSize, FloatType intervalSize)
{
//  FloatType intervalSize;
  FloatType start = graph.front().first;
  FloatType end   = graph.back().first;
  size_t intervalCount = static_cast<size_t>((end - start) / stepSize);
  //cout << "Bereich " << start << " bis " << end << ", Intervals: " << intervalCount << endl;
  FloatType endOfGraph = graph.back().first;
  
  FloatPairList averagedPlot;

  // Compute the size of an interval
  size_t previousSizeOfQueue = 0;
  FloatType previousQueueSum = 0.0;

  //Initialize graph point iterator
  std::vector<FloatPair>::const_iterator currentPoint1 = graph.begin();
  std::vector<FloatPair>::const_iterator currentPoint2 = graph.begin();
  
  FloatPair outOfBounce = FloatPair(10.0,0);
  
  // For each interval
  for (size_t currentInterval = 0; currentInterval < intervalCount; ++currentInterval) {
    // Update current interval start and end 
    FloatType intervalBegin = start + currentInterval * stepSize - intervalSize/2.0;
    FloatType intervalEnd   = start + currentInterval * stepSize + intervalSize/2.0;
    
    size_t    newPointsCount = 0;   // number of new points entering the queue
    FloatType newPointsSum   = 0.0; // sum of these points
    size_t oldPointsCount = 0;      // number of points to pushed out of the queue
    FloatType oldPointsSum = 0;  
    while (currentPoint1->first < intervalBegin) {
      ++oldPointsCount;
      oldPointsSum += currentPoint1->second;
      ++currentPoint1;
    }
    while ((currentPoint2 != graph.end()) && (currentPoint2->first < intervalEnd)) {
        ++newPointsCount;
        newPointsSum += currentPoint2->second;
        ++currentPoint2;
    }

    // Compute central x-point
    FloatType centralX;
    if (currentPoint2 == graph.end()) centralX = (intervalBegin + (currentPoint2-1)->first)/2; // So the interval gets smaller when we approach the end of the curve
    else centralX = (intervalBegin + intervalEnd) /2;    
    previousQueueSum    -= oldPointsSum;
    previousQueueSum    += newPointsSum; // Set the new sum
    previousSizeOfQueue -= oldPointsCount;
    previousSizeOfQueue += newPointsCount;
    if (previousSizeOfQueue < 3) {
      //cout << "Too little elements, we skip that point at " << centralX << "  " << previousSizeOfQueue << endl;
      continue;
    }
//     averagedPlot->push_back(FloatPair(centralX, previousQueueSum/previousSizeOfQueue));
    averagedPlot.push_back(FloatPair(centralX, previousQueueSum/previousSizeOfQueue));
  }
  return averagedPlot;
}

// Replace this by some sort of pair accumlate
// @todo accumulate(currentStartPoint, currentEndPoint, make_pair(0.0, 0.0), ); // 
FloatType sliceSum(const std::vector<FloatPair>& graph, std::vector<FloatPair>::const_iterator start, std::vector<FloatPair>::const_iterator end)
{
  FloatType slSum = 0.0;
  for (std::vector<FloatPair>::const_iterator itr = start; itr != end; ++itr) {
    slSum += itr->second;
  }
  return slSum;
}

/**
 * Another function to calculate the moving average. Here the moving window shrinks
 * when we approach the right part of the curve.
 * 
 * @note Interestingly, the results in the left part do not always seem to agree to 100% with the other average method.
 * @param graph 
 * @param windowSize Size of the sliding window
 * 
 */
FloatPairListPtr MathUtil::movingAv2(const std::vector<FloatPair>& graph, size_t windowSize, size_t minimalWindowSize)
{
  size_t firstAveragedDatapoint = windowSize/2 + 1 - 1;  
  size_t winSize = windowSize; // need to copy it, since it won't stay constant during this method (getting smaller in the end)
  FloatPairListPtr averagedPlot(new FloatPairList());
  FloatType sum = 0.0;

  //Initialize graph point iterator
  std::vector<FloatPair>::const_iterator currentStartPoint = graph.begin(); // start of the interval
  std::vector<FloatPair>::const_iterator currentMidPoint   = graph.begin() + firstAveragedDatapoint; // middle where the averaged value is place
  std::vector<FloatPair>::const_iterator currentEndPoint   = graph.begin() + winSize; // end of the interval.
  
  // Calculate the sum of the points inside the first interval.
  FloatType tmpSum = sliceSum(graph, currentStartPoint, currentEndPoint);

  // Store first average
  averagedPlot->push_back(FloatPair(currentMidPoint->first, tmpSum/winSize)); 

  // Loop until window arrives at end point
  while(currentEndPoint != graph.end()) {
    // The average is not caluclated for every new slice, instead the sum over the array is 
    // updated by substracting and adding the elements that differ from the previous array.
    tmpSum    -= currentStartPoint->second;
    tmpSum    += currentEndPoint->second; // Note that endPoint is always lastElement + 1

    // Advance window + 1
    ++currentMidPoint; ++currentStartPoint; ++currentEndPoint; 

    // Place the average at the mid point
    averagedPlot->push_back(FloatPair(currentMidPoint->first, tmpSum/winSize));  
  }
 

  // For the remaining elements, each iteration shrink window size by 2 
  // and move the middle element one to the right (so windows is centered around a single middle element)
  while (winSize > minimalWindowSize) {
    winSize = winSize - 2;

    // Move two further
    tmpSum    -= currentStartPoint->second;
    ++currentStartPoint;
    tmpSum    -= currentStartPoint->second;
    ++currentStartPoint;

    ++currentMidPoint;

    averagedPlot->push_back(FloatPair(currentMidPoint->first, tmpSum/winSize));  
  }

  return averagedPlot;
}




/**
 * Using a simple heuristic, the function decides wheather or not 
 * a graph (a hookcurve) descends in its rightermost part. For this,
 * the graph is digitized into very few points. It is assumed to 
 * descend if the last point is smaller than the one before.
 *
 *
 * @param graph Sorted set of points in the original graph.
 */
bool MathUtil::graphDescends(const std::vector<FloatPair>& graph)
{
  size_t pointCount = 5;
  std::vector<FloatPair> digitizedCurve = MathUtil::digitizeCurve(graph, pointCount );
  return digitizedCurve[pointCount-1].second < digitizedCurve[pointCount-2].second;
}

/**
 * Solves a system of equations Ax=b using singular value decompsition.
 *
 * @note Uses vector's array values directly, which might not be clean.
 *
 * 
 * @param delta Matrix A of equation system.
 * @param theta Vector b of equation system.
 * @param zeroingCutoff All singular values below this constant are set to zero
 *
 * @return Solution of equation (x).
 *
 * @todo Make parameters const
 * @todo Assert correct sizes
 * @todo Document Zeroing.
 */
std::valarray<FloatType> MathUtil::solveEquationsSvd(math::Matrix<FloatType>& delta, std::valarray<FloatType>& theta,
                                                  const FloatType zeroingCutoff)
{
//math::Matrix<FloatType> delta(5, 4, 0.0);
//gsl_matrix_view deltaView = gsl_matrix_view_array(&(delta.flat()[0]), 4, 5);
//gsl_matrix* u = gsl_matrix_alloc(4, 5);
//gsl_matrix_set_all(u,0.0);
//delta[0][0] = 0
//print(delta);
//   FloatType delta_[2][2] = {{3.0, 2.0},
//                          {2.0, 5.0}};

//   math::Matrix<FloatType> delta_(2,2);

//   delta_[0][0] = 3.0;
//   delta_[0][1] = 2.0;
//   delta_[1][0] = 4.0;
//   delta_[1][1] = 5.0;

//   print(delta_);
//   valarray<FloatType> tmp2 = delta_.flat();
//   for (size_t i = 0; i < tmp2.size(); ++i) {
//     cout << tmp2[i] << " ";
//   }
//   cout << endl;

//   FloatType theta_[2] = {7.0, 14.0};
  //cout<<"This is delta\n";
  //print(delta);
  //cout<<"\n\n\nAnd this is theta\n";
  //for (int i =0; i<theta.size(); ++i){
   // if (i%25==0){
   // cout<<theta[i]<<"\t";
   // }
   // }

  //print(delta);

  size_t equationCount = delta.rowNum();
  size_t variableCount = delta.columnNum(); // RECHECK

  // Create solution vector
  // std::valarray<FloatType> x(variableCount);
  // std::vector<FloatType> x(variableCount);

  // Create gsl views on the matrices
//   gsl_matrix_view deltaView = gsl_matrix_view_array(&delta_(0,0), equationCount, variableCount); // RECHECK
  gsl_matrix_view deltaView = gsl_matrix_view_array(&(delta.flat()[0]), equationCount, variableCount); // RECHECK
  gsl_vector_view thetaView = gsl_vector_view_array(&theta[0], equationCount);
  //  gsl_vector_view xView = gsl_vector_view_array(&x[0], variableCount);
  gsl_vector* x = gsl_vector_alloc (variableCount); 

  // gsl_vector_int_const_view v = gsl_vector_int_const_view_array(&stl_v[0], stl_v.size());

  gsl_matrix* u = gsl_matrix_alloc(equationCount, variableCount); // (rows, colums)
  gsl_matrix* v = gsl_matrix_alloc(equationCount, variableCount); 
  gsl_vector* s = gsl_vector_alloc(variableCount);
  gsl_vector* tmp = gsl_vector_alloc(variableCount);    



  // Solve linear equations with SVD
  gsl_matrix_memcpy(u, &deltaView.matrix); // SVD overwrites original matrix with u
  gsl_linalg_SV_decomp(u, v, s, tmp);

  //zeroing: for (almost) singular matrices some values of v get very small which often leads to large
  //solution vectors. to avoid this, these small elements are set to zero which makes the gsl-svd routine
  //to ignore these entries. often this leads to better results (in our case it obviously does) (see numerical recipies for details)
  // matrix s is saved as a vector since it only holds diagonal elements.
  for (size_t i = 0; i < equationCount; ++i){
    if (gsl_vector_get(s,i) < zeroingCutoff){
      gsl_vector_set(s,i,0.0);
    }
  }

  gsl_linalg_SV_solve(u, v, s, &thetaView.vector, x);    

  // Free solution vector memory but save it before
  //std::valarray<FloatType> solution(gsl_vector_ptr(x,0), variableCount);
  std::valarray<FloatType> solution(variableCount);
  for (size_t i = 0; i < variableCount; ++i) {
    solution[i] = gsl_vector_get(x,i);
  }
  gsl_vector_free(x);

  // Free memory (SVD)
  gsl_matrix_free(u); gsl_matrix_free(v);
  gsl_vector_free(s); gsl_vector_free(tmp);
  

  // cout << gsl_vector_get(x,0) << " " << gsl_vector_get(x,1) << endl;
  return solution;
}


/**
 * Solves a system of equations Ax=b using LU
 *
 * @note Uses vector's array values directly, which might not be clean.
 *
 * @todo Make parameters const
 * 
 * @param delta Matrix A of equation system.
 * @param theta Vector b of equation system.
 * @return Solution of equation (x).
 *
 * @todo Assert correct sizes
 */
std::valarray<FloatType> MathUtil::solveEquationsLU(math::Matrix<FloatType>& delta, std::valarray<FloatType>& theta)
{
  size_t equationCount = delta.rowNum();
  size_t variableCount = delta.columnNum(); // RECHECK

  // Create gsl views on the matrices
//   gsl_matrix_view deltaView = gsl_matrix_view_array(&delta_(0,0), equationCount, variableCount); // RECHECK
  gsl_matrix_view deltaView = gsl_matrix_view_array(&(delta.flat()[0]), equationCount, variableCount); // RECHECK
  gsl_vector_view thetaView = gsl_vector_view_array(&theta[0], equationCount);
  //  gsl_vector_view xView = gsl_vector_view_array(&x[0], variableCount);
  gsl_vector* x = gsl_vector_alloc (variableCount); 

  // gsl_vector_int_const_view v = gsl_vector_int_const_view_array(&stl_v[0], stl_v.size());

  // Allocate LU specific memory
  gsl_matrix* lu = gsl_matrix_alloc (equationCount, variableCount); 
  gsl_permutation* p = gsl_permutation_alloc (variableCount); 
  int signum;

  // Solve linear equations with LU
  gsl_matrix_memcpy(lu, &deltaView.matrix); // LU decomp overwrites original matrix with LU
  gsl_linalg_LU_decomp(lu, p, &signum);
  gsl_linalg_LU_solve(lu, p, &thetaView.vector, x);



  // Free solution vector memory but save it before
  //std::valarray<FloatType> solution(gsl_vector_ptr(x,0), variableCount);
  std::valarray<FloatType> solution(variableCount);
  for (size_t i = 0; i < variableCount; ++i) {
    solution[i] = gsl_vector_get(x,i);
  }
  gsl_vector_free(x);

  // Free memory (LU)
  gsl_permutation_free(p);
  gsl_matrix_free(lu); 

  return solution;
}

/**
 * Fits a straight line a*x+b into the (sub)graph using least squares fitting.
 *
 * @param beginPoint First point of the graph
 * @param endPoint End point of the graph.
 *
 * @return An array [a,b] containing the parameters of the line ax+b
 */
Polynomial<FloatType> MathUtil::fitStraightLine(const GraphIterator& beginPoint, const GraphIterator& endPoint)
{
  // Initialize system of linear equations
  math::Matrix<FloatType> delta(2,2, 0.0);
  valarray<FloatType> theta(0.0, 2);

  // Fill the system for linear fit
  for (GraphIterator currentPoint = beginPoint; currentPoint != endPoint; ++currentPoint) {
    FloatType x = currentPoint->first;
    FloatType y = currentPoint->second;
    delta[0][0] += x*x;
    delta[0][1] += x;
    delta[1][0] += x;
    delta[1][1] += 1;
    theta[0] += x*y;
    theta[1] += y;
  }

  // Solve the system and return ascent
  return Polynomial<FloatType>(MathUtil::solveEquationsLU(delta, theta));
}

void printStraightLine(std::valarray<FloatType> line, FloatType startPoint, FloatType endPoint, FloatType step,
                       std::ostream& stream)
{
  for(FloatType x = startPoint; x < endPoint; x += step) {
    stream << x << "\t" << (line[0]*x + line[1]) << endl;
  }
}

/**
 * Calculates the point before the leftmost steepest ascend of a 
 * graph.
 *
 * @todo describe algorithm precisely
 *
 * @note: This is another possibility to calculate the breakpoint
 * between specific and non-specific binding in the hookcurve.
 * 
 * @param graphPoints The (x,y) pairs in the graph sorted ascending in x.
 * @param Number of intervals the graph is devided into each iteration.
 *
 */
// FloatType MathUtil::getPointBeforeSteepestAscent(const std::vector<FloatPair>& graphPoints, 
//                                                  size_t intervalCount)
// {
//   // Debug for "plot "testplot_shuffle.dat" using 1:2, "breakpointDebug.dat" using 1:2"
//   ofstream debugFile;   // DEBUG
// 	debugFile.open("breakpointDebug.dat");   // DEBUG


//   assert(graphPoints.size() > 0);
//   assert(intervalCount > 0);

//   const int minimalPoints = 7; // Stop iterating if less points in an interval left
//   const FloatType rubustnessFactor = 0.7;

//   // Iterate as long enough points left
//   GraphIterator beginPoint = graphPoints.begin(),
//     rightPoint = graphPoints.end() - 1;
//   while (rightPoint - beginPoint > minimalPoints) {
//     cout << "beginPoint: " << beginPoint->first << endl;
//     cout << "rightPoint: " << rightPoint->first << endl;    

//     // Get interval borders and ascents
//     vector< pair<GraphIterator, FloatType> > intervalsAndAscends = MathUtil::getIntervalsAndAscents(beginPoint, rightPoint, intervalCount);

//     // Get find leftmost robust maximum Index 
//     // @todo Either use stl functions or merge with above
//     size_t maxIndex = intervalsAndAscends.size() - 1;
//     for (int i = intervalsAndAscends.size() - 2; i >= 0; --i) {
//       if (atan(intervalsAndAscends[i].second) > rubustnessFactor*atan(intervalsAndAscends[maxIndex].second)) {
//         maxIndex = i;
//       }
//     } 


//     // Stop if maximum is left of maxIndex
// //     if (maxIndex == 0) {
// //       cout << "Aborting breakpoint search" << endl;
// //       break;
// //     }
      
//     // Update left and right
//     if (maxIndex > 1) { // Left point does not change iff the steepest ascent is within the first two intervals
//       beginPoint = intervalsAndAscends[maxIndex-2].first; // interval[0].first RIGHTMOST point interval 0
//     } 
//     rightPoint = intervalsAndAscends[maxIndex].first;


//     // Print ascent to file
//     printStraightLine(MathUtil::fitStraightLine(beginPoint, rightPoint + 1), beginPoint->first, rightPoint->first, 0.001, debugFile);
//   for (FloatType y = 0; y < 0.2; y += 0.005) {
//     debugFile << beginPoint->first << "\t" << y << endl;
//   }

    
// //     debugFile << beginPoint->first << "\t" << beginPoint->second << endl;
// //     debugFile << rightPoint->first << "\t" << rightPoint->second << endl;
//   }
//   debugFile.close();

//   cout << "beginPoint->first: " << beginPoint->first << endl;
  
//   return beginPoint->first;
// }


/**
 * Trims a value v to the range \f$ minValue \leq v \leq maxValue \f$.
 *
 * @param value The value
 * @param minValue Smallest desired value for v.
 * @param maxValue The biggest desired value for v.
 * 
 * @return Value trimmed to desired range.
 */
FloatType MathUtil::trimToRange(FloatType value, FloatType minValue, FloatType maxValue)
{
  return max(min(value, maxValue), minValue);
}

/**
 * Returns a value for NaN (Not a Number), if supported by the current platform.
 * Otherwise returns 0.
 *
 * @return NaN 
 */
FloatType MathUtil::getNaN() 
{
  if (numeric_limits<FloatType>::has_signaling_NaN) {
    return numeric_limits<FloatType>::signaling_NaN();
  }
//   else if(numeric_limits<FloatType>::has_quiet_NaN) {
//     d=numeric_limits<FloatType>::quiet_NaN();
//   }

// btw. for later:
// if(numeric_limits<float>::has_infinity)
//  f=numeric_limits<float>::infinity();
  assert( !numeric_limits<FloatType>::has_signaling_NaN );
  return 0;
}


/**
 * Returns a value for NaN (Not a Number), if supported by the current platform.
 * Otherwise returns 0.
 *
 * @return NaN 
 */
size_t MathUtil::getNaNSizeT() 
{
  if (numeric_limits<size_t>::has_signaling_NaN) {
    return numeric_limits<size_t>::signaling_NaN();
  }
//   else if(numeric_limits<FloatType>::has_quiet_NaN) {
//     d=numeric_limits<FloatType>::quiet_NaN();
//   }

// btw. for later:
// if(numeric_limits<float>::has_infinity)
//  f=numeric_limits<float>::infinity();
  assert( !numeric_limits<size_t>::has_signaling_NaN );
  return 0;
}


/**
 * Computes the minimal value of an multivariate fitness function f using
 * the method of gradient descent.
 *
 * @param initialPoint Starting point for the algorithm.
 * @param f The function to optimize.
 */
std::valarray<FloatType> MathUtil::computeGradientDescent(const std::valarray<FloatType>& initialPoint, 
                                                          const MathOptimizationFunction& f)
{
  // Initialize constants 
  FloatType differentialDistance = 0.000000001;
  size_t dimension = initialPoint.size();

  // Initialize variables for aborting conditions
  const size_t maxIterations = 5000;
  size_t iterationCount = 0;

//   std::valarray<FloatType>optimum(getNaN(), dimension);
  std::valarray<FloatType>optimum(0.0, dimension);
  FloatType optimumDistanceAbort = 0.00001;
  FloatType mingradientAbort = 0.0000000001;

  // Compute starting values
  std::valarray<FloatType> currentPoint = initialPoint;
  FloatType currentValue = f(currentPoint);
  std::valarray<FloatType> numericGradient(mingradientAbort + 1.0, dimension); // Initialize with something > maxgradientAbort
  
  // Iterate until abortion coditions fulfilled:
  while ((iterationCount < maxIterations) 
         && (isnan(optimum.max()) || std::valarray<FloatType>(currentPoint - optimum).max() > optimumDistanceAbort) // if optimum defined, abort iff close
         /* && (abs(numericGradient.max()) > mingradientAbort) */ ) {

    // Compute (numerical) gradient for each dimension
    for (size_t i = 0; i < dimension; ++i) {
      std::valarray<FloatType> differentialPoint = currentPoint;
      differentialPoint[i] += differentialDistance;
      numericGradient[i] = (f(differentialPoint) - currentValue) / differentialDistance;
    }

//     for (size_t i = 0; i < 100000; ++i) FloatType t = sqrt(i);
//     cout << "iterationCount: " << iterationCount << endl;
//     cout << "currentPoint:                     "<< currentPoint << endl;
//     cout << "numericGradient: " << numericGradient << endl;
    
    // Compute new point
    FloatType gradientNorm = std::valarray<FloatType>(numericGradient * numericGradient).sum();
    std::valarray<FloatType> decayFactor = 0.1*abs(currentPoint)/(1 + 0.001*iterationCount);
//     cout << "decayFactor: " << decayFactor << endl;
    currentPoint -= decayFactor * numericGradient/sqrt(gradientNorm);
    currentValue = f(currentPoint);

    ++iterationCount;
  }

  // Debug: Plot the optimizing function
//   ofstream gnuplotFile;
//   std::valarray<FloatType> testPoint(2);
// 	gnuplotFile.open("landscape.dat");
//   for (FloatType a = 0.001; a <= 1; a+= 0.0001 + a/100) {
//     for (FloatType ff = 0.00001; ff <= 0.01; ff+= 0.00001 + ff/100) {
//       testPoint[0] = a; testPoint[1] = ff;
//       gnuplotFile << a << "\t" << ff << "\t" << (f(testPoint)) << endl;
//     }
//   }
  
  ofstream gnuplotFile;
  std::valarray<FloatType> testPoint(2);
  gnuplotFile.open("landscape.dat");
  for (FloatType a = -1; a <= 0; a+= 0.01){// + a/100) {
    for (FloatType ff = -5; ff <= -1; ff+= 0.01){// + ff/100) {
      testPoint[0] = exp10(a); testPoint[1] = exp10(ff);
      gnuplotFile << a << "\t" << ff << "\t" << (log10(f(testPoint))) << endl;
    }
  }
  
  gnuplotFile.close();

  return currentPoint;
}


/**
 * Trims a value v to the range \f$ minValue \leq v \leq maxValue \f$.
 *
 * @param value The value
 * @param minValue Smallest desired value for v.
 * @param maxValue The biggest desired value for v.
 * 
 * @return Value trimmed to desired range.
 */
/*FloatType MathUtil::trimToRange(FloatType value, FloatType minValue, FloatType maxValue)
{
  return max(min(value, maxValue), minValue);
}*/


/**
 * Returns the function value of the normal distribution
 *   \f[ N = \frac{1}{\sqrt{2 \pi} \sigma} \exp{- \frac{(x - \mu)^2}{2 \sigma^2}} \f]
 *
 * @param x X
 * @param mu \f$ \mu \f$
 * @param sigma \f$ \sigma \f$
 */
FloatType MathUtil::normalDistribution(FloatType x, FloatType mu, FloatType sigma)
{
  return (1.0/(sqrt(2 * M_PI) * sigma)) * exp( - (pow((x - mu) , 2) / (2 * pow(sigma , 2)) ));
}


/**
 * Returns the function value of the simplified bivariate normal distribution
 *   \f[ N = \frac{1}{2 \pi \sigma^2 \sqrt{1 - \rho^2}} \exp{-\frac{(x - \mu)^2}{2 \sigma^2}} \f]
 *
 * @param x X
 * @param mu \f$ \mu \f$
 * @param sigma \f$ \sigma \f$
 * @param rho \f$ \rho \f$
 */

FloatType MathUtil::simpleBivariateDistribution(FloatType x, FloatType mu, FloatType sigma, FloatType rho)
{
	FloatType sSquare = pow(sigma, 2);
	FloatType ret = (1.0/(2 * M_PI * sSquare * sqrt(1 - pow(rho, 2)))) // This is to see the value in debugging modus
    * exp( - pow((x - mu), 2) / (2 * sSquare));
  return ret;
}


/// Helper functions for the gsl integration:
double MathUtil::simpleDistributionIntegral(double x, void *p) {
  const struct SimpleIntegralParameters* params 
    = (struct SimpleIntegralParameters *)p;
  return MathUtil::normalDistribution(x, params->mean, params->std);
}

/// Helper function for PM or MM integration
double MathUtil::simpleIntegralFunction(double x, void *p) {  
  struct SimpleIntegralParameters* params 
    = (struct SimpleIntegralParameters *)p;
  //   return MathUtil::gLog10( ((1.0/sqrt(2 * M_PI * pow(std , 2))) * exp( - (pow((x - mu) , 2) / (2 * pow(std , 2)) ))) * (intensity - exp10(x) * deltaNS) );
  double distrib = MathUtil::normalDistribution(x, params->mean, params->std);
  double expd = exp10(x + params->sensitivity);
  double glogged = MathUtil::gLog10((params->desaturatedIntensity - expd ));
  return distrib * glogged;
}

/// Helper function for the gsl bivariate integration:
double MathUtil::simpleBivariateDistributionIntegral(double x, void *p) {
  struct SimpleBivariateIntegralDistributionParameters* params 
    = (struct  SimpleBivariateIntegralDistributionParameters*)p; 
  return MathUtil::simpleBivariateDistribution(x, params->mean, params->std, params->correlationCoefficient);
}

// Function for simple bivariate integration
double MathUtil::simpleBivariateIntegralFunction(double x, void *p) {  
  struct SimpleBivariateIntegralParameters* params 
    = (struct SimpleBivariateIntegralParameters *)p;
  double mu          = (params->mean);
  double std         = (params->std);
  double intensityPm = (params->desaturatedIntensityPm);
  double intensityMm = (params->desaturatedIntensityMm);
  double yNSPm       = (params->sensitivityPm); //exp10(params->sensitivity);
  double yNSMm       = (params->sensitivityMm);
  double rho         = (params->correlationCoefficient);
  double logB        = (params->logB);
  //   return MathUtil::gLog10( ((1.0/sqrt(2 * M_PI * pow(std , 2))) * exp( - (pow((x - mu) , 2) / (2 * pow(std , 2)) ))) * (intensity - exp10(x) * deltaNS) );
  //   if (simpleBivariateDistribution(x, mu, std, rho) * MathUtil::gLog10(((intensityPm - exp10(yNSPm + x)) 
  //       - (intensityMm - exp10(yNSMm + mu*(1-rho) - logB + rho * x))), c) 
  //       == MathUtil::getNaN()) {cout << "NaN: "<<endl;}
  
  //   cout << mu << "\t" << std << "\t" << intensityPm << "\t" << intensityMm << "\t" << yNSPm << "\t" << yNSMm << "\t" << c << "\t"      
  //       << rho << "\t" << logB << endl;
  //   cout << simpleBivariateDistribution(x, mu, std, rho) <<endl;
  //   cout << MathUtil::gLog10((
  //                            (intensityPm - exp10(yNSPm + x))
  //       - (intensityMm - exp10(yNSMm + mu*(1-rho) - logB + rho * x)))
  //       , c) << endl;
  double distrib = MathUtil::simpleBivariateDistribution(x, mu, std, rho);
  double deltaLPm = (intensityPm - exp10(yNSPm + x));
  double deltaLMm = (intensityMm - exp10(yNSMm + mu*(1-rho) - logB + rho * x));
  double glogged = MathUtil::gLog10(( deltaLPm - deltaLMm));
  return distrib * glogged;
}



// Function for simple bivariate integration
double MathUtil::gcrmaLikeBivariateIntegralFunction(double x, void *p) {  
  struct SimpleBivariateIntegralParameters* params 
    = (struct SimpleBivariateIntegralParameters *)p;
  double mu          = (params->mean);
  double std         = (params->std);
  double intensityPm = (params->desaturatedIntensityPm);
  double intensityMm = (params->desaturatedIntensityMm);
  double yNSPm       = params->sensitivityPm; //exp10(params->sensitivity);
  double yNSMm       = params->sensitivityMm;
  double rho         = params->correlationCoefficient;
  double logB        = params->logB;
  FloatType distrib = MathUtil::simpleBivariateDistribution(x, mu, std, rho);
  FloatType exponentenKram = exp10(yNSPm - yNSMm + (x - mu) * (1 - rho) + logB);
  FloatType glogged = MathUtil::gLog10((intensityPm - intensityMm * exponentenKram));
  return distrib * glogged;
  //MathUtil::gLog10(intensity - exp10(x) * yNS);    
}
