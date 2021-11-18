/**
 * @file AveragingProcedure.hpp Declares a family of classes for computing the average 
 *       of some values.
 * @author $Author$ 
 * @author Mario Fasold
 * @date $Date$
 *
 */
#ifndef _AVERAGINGPROCEDURE_
#define _AVERAGINGPROCEDURE_

#include <valarray>
#include <numeric>
#include "math.h"

// // DEBUG
// #include <iostream>
// #include "ValarrayUtil.hpp"

namespace larrpack {

/**
 * @class ValarrayFunctor
 * @brief Abstract class for function objects (=functors) that 
 * works on valarrays. A functor can be regarded as an  elegant 
 * way to implement the strategy pattern in c++. The functor 
 * can be used like a function, i.e. functor(paramters).
 *         
 * @note The main purpose of this class is to make the design 
 *       and the documentation (->class hierarchy) more obvious.
 * @note The behaviour of the templated functions has not been evaluated with
 *       types other than floating-point types. 
 * @note Functions do not perform size(lists) > 0 test.
 *
 * @todo Add a destinct result template variable. (i.e. result type of 
 *       Median(int) != int -> Classes should inherit from ValarrayFunctor<Arg,Res>)
 */
template<class ArgumentType, class ReturnType>
class ValarrayFunctor : public std::unary_function<std::valarray<ArgumentType>, ReturnType>
{
};



/**
 * @class Mean
 * @brief Function object to compute (arithmetic) mean
 * \f$ \frac{1}{n} \sum_{i=1}^n value_i \f$.
 */
template<class T>
class Mean : public ValarrayFunctor<T,T>
{
public:
  /**
   * Computes the mean \f$ \frac{1}{n} \sum_{i=1}^n value_i \f$
   *
   * @param values List of values.
   * @return Mean of values.
   *
   */
  T operator() (const std::valarray<T>& values) const
  {
    return values.sum() / (T)values.size();
  }
};

/**
 * @class Median 
 * @brief Computes the median.
 * 
 */
template<class T>
class Median : public ValarrayFunctor<T,T>
{
public:
  /**
   * Compute the median.
   *
   * @param values List of values.
   * @return Median of values.
   *
   */
  T operator() (const std::valarray<T>& values) const
  {
    // Create a copy of the list and sort it
    std::valarray<T> sortedValues(values);
    size_t valueCount = sortedValues.size();
    std::sort(&sortedValues[0], &sortedValues[valueCount]); // Can use pointer as iterator

    // Return middle element if unequal number of elements
    if (valueCount % 2 == 1) {
      return sortedValues[(valueCount - 1)/ 2];
    }

    // Return average of the two middle elements otherwise
    return (sortedValues[(valueCount/2) - 1] + sortedValues[(valueCount/2)]) / (T)2.0;
  }
};


/**
 * @class OnestepTukeyBiweight
 * @brief Computes One-Step Tukey Biweight.
 * 
 * @see http://www.affymetrix.com/support/technical/whitepapers/sadd_whitepaper.pdf
 */
template<class T>
class OnestepTukeyBiweight : public ValarrayFunctor<T,T>
{
private:
  /// Tuning constant
  double mC;

  /// Small constant that prevents division by zero
  double mEpsilon;
public:
  /// Constructor
  explicit OnestepTukeyBiweight(const double c = 5.0, const double epsilon = 0.0001) : 
    mC(c), mEpsilon(epsilon) {}

  /**
   * Computes One-Step Tukey Biweight.
   *
   * @see http://www.affymetrix.com/support/technical/whitepapers/sadd_whitepaper.pdf
   *
   * @param values List of values.
   * @return Biweight of values.
   *
   * @todo Compare with results from statistic programs.
   */
  T operator() (std::valarray<T>& values) 
  {
    // Calculate median
    Median<T> medianCalculator;
    T median = medianCalculator(values);
    
    // Calculate S = median of distances of every value from median
    T S = medianCalculator(abs(values - median));
    //T S = medianCalculator(values - median);

    // Calculate u and weight(u)
    std::valarray<T> u = (values - median) / (T)(mC*S + mEpsilon);

    std::valarray<T> weight((T)0, values.size()); // Initialize with 0
    for(size_t i = 0; i < values.size(); ++i) {
      if (fabs(u[i]) <= 1) {
        weight[i] = (1 - u[i]*u[i]) * (1 - u[i]*u[i]);
      }
    }
    
    // Calculate weighted average
    T weightedSum = 0;
    for(size_t i = 0; i < values.size(); ++i) {
      weightedSum += weight[i] * values[i];
    }    
    return weightedSum / weight.sum();
  }
};

}
#endif
