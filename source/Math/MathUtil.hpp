/**
 * @file MathUtil.hpp Provides various mathematical functions.
 * @author $Author$ 
 * @author Mario Fasold
 * @date $Date$
 *
 */
#ifndef _MATHUTIL_
#define _MATHUTIL_

#include <vector>
#include <valarray>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <sstream>

#include "Matrix.hpp"
#include "Polynomial.hpp"

// Define functions for Microsoft VC++ compatibility
// @todo: can we use general HAVE_DECL_EXP10 directive?
#ifdef _MSC_VER
// MS VC++ does not provide ::exp10 (conrary to gcc stdlib),so
// we have to provide our own version.
inline double exp10(double d) {
  return exp(d * log(10.0));
}

// isnan function isn't provided either
#define isnan(x) ((x) != (x))
#endif 

/**
 * @namespace MathUtil
 * @brief Provides utility functions for mathematical calculations.
 * 
 * @todo Komplett auf valarray umstellen
 */
namespace MathUtil {
  typedef double FloatType;
  typedef std::pair<FloatType, FloatType> FloatPair;
  typedef std::vector<FloatPair> FloatPairList;
  typedef boost::shared_ptr<FloatPairList> FloatPairListPtr;
  typedef std::vector<FloatPair>::const_iterator GraphIterator;  
  
  struct SimpleIntegralParameters {
    double mean;
    double std;
    double desaturatedIntensity;
    double sensitivity;    
  };

  struct SimpleBivariateIntegralParameters {
    double mean;
    double std;
    double desaturatedIntensityPm;
    double desaturatedIntensityMm;
    double sensitivityPm;    
    double sensitivityMm;    
    double correlationCoefficient;
    double logB;
  };

  struct SimpleBivariateIntegralDistributionParameters {
    double mean;
    double std;
    double correlationCoefficient;
  };

  // Calculates a moving average where size(result) = size(input)
  std::valarray<FloatType> calculateGeneralizedMovingAverage(const std::valarray<FloatType> &values, 
                                                             size_t movingAverageWindowSize);
    
  FloatPairListPtr  calculateMovingAverage(const FloatPairList& values, size_t movingAverageWindowSize);
  
  FloatPairListPtr movingAv2(const std::vector<FloatPair>& graph, size_t windowSize, size_t minimalWindowSize = 4);

  FloatType calculateMean(const std::vector<FloatType>& values);
  FloatType calculateStandardDeviation(const std::vector<FloatType>& values) ;

  FloatType calculateMean(const std::valarray<FloatType>& values);
  FloatType calculateStandardDeviation(const std::valarray<FloatType>& values) ;
  
  FloatType calculateMedian(const std::valarray<FloatType>& values);
  // Calculates the MAD (median absolute deviation) of some float data
  FloatType calculateMAD(const std::valarray<FloatType>& values);

  FloatType calculateCovariance(const std::vector<FloatType>& valuesX, 
                                const std::vector<FloatType>& valuesY);

  /// Calculate covaricance with known means
  FloatType calculateCovariance(const std::vector<FloatType>& valuesX, 
                                const std::vector<FloatType>& valuesY, 
                                FloatType muX, FloatType muY);
  FloatType calculateCorrelation(const std::vector<FloatType>& valuesX, 
                                 const std::vector<FloatType>& valuesY);
  /// Calculate Correlation with known means and stds
  FloatType calculateCorrelation(const std::vector<FloatType>& valuesX, 
                                 const std::vector<FloatType>& valuesY, 
                                 FloatType muX, FloatType muY, FloatType stdX, FloatType stdY);

  /// Parameter c used for gLog computation
  extern FloatType gGLogParameter;
  
  FloatType gLog10(FloatType z);
  FloatType gExp10(FloatType log10Value);
//   FloatType exp10(FloatType log10Value);
  
  FloatType getGaussianWeight(const FloatType x, const FloatType mu, const FloatType sigma);
  

  std::vector<FloatPair> digitizeCurve(const std::vector<FloatPair>& graph, size_t intervalCount,
                                       FloatType sigma = 5.0);
  std::vector<FloatPair> digitizeCurveRobust(const std::vector<FloatPair>& graph, 
                                             size_t intervalCount, FloatType sigma = 5.0);
  
  std::vector<FloatPair> digitizeCurveNew(const std::vector<FloatPair>& graph, FloatType stepSize = 0.03, FloatType averageSize = 0.1);
  
  bool graphDescends(const std::vector<FloatPair>& graph);

  /// Solves an a system of linear equations with the SVD method
  std::valarray<FloatType> solveEquationsSvd(math::Matrix<FloatType>& delta, 
                                             std::valarray<FloatType>& theta,
                                             const FloatType zeroingCutoff = 0.1);

  /// Solves an a system of linear equations with the LU method
  std::valarray<FloatType> solveEquationsLU(math::Matrix<FloatType>& delta, 
                                            std::valarray<FloatType>& theta);  


  larrpack::Polynomial<FloatType> fitStraightLine(const GraphIterator& leftPoint, 
                                                  const GraphIterator& rightPoint);

  FloatType trimToRange(FloatType value, FloatType minValue, FloatType maxValue);

  FloatType getNaN();  
  size_t getNaNSizeT();

  /// A function object to create the probeset composition
  typedef boost::function1<FloatType, std::valarray<FloatType> > MathOptimizationFunction;

  
  std::valarray<FloatType> computeGradientDescent(const std::valarray<FloatType>& initialPoint, 
                                                  const MathOptimizationFunction& f);
  
  // Normal distribution Function
  FloatType normalDistribution(FloatType x, FloatType mu, FloatType sigma);
  
  // Simple Bivariate Distribution
  FloatType simpleBivariateDistribution(FloatType x, FloatType mu, FloatType sigma, FloatType rho);
  
  // Other Integral related function:
  double simpleDistributionIntegral(double x, void *p);
  
  double simpleIntegralFunction(double x, void *p);
  
  double simpleBivariateDistributionIntegral(double x, void *p);
  
  double simpleBivariateIntegralFunction(double x, void *p);
  
  double gcrmaLikeBivariateIntegralFunction(double x, void *p);


} // Namespace
#endif
