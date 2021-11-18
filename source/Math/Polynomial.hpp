/**
 * @file Polynomial.hpp Contains a simple class to compute polynomials.
 * @author $Author: SVN $ 
 * @author Mario Fasold
 * @date $Date: 2007-05-10 17:55:05 +0200 (Do, 10 Mai 2007) $
 *
 */
#ifndef _POLYNIMIAL_
#define _POLYNIMIAL_

#include <sstream>

namespace larrpack {
  /**
   * @class Polynomial
   * @brief Data structure for an univariate polynomial 
   * \f$ y(x) = a_0*x^n + ... + a_{n-2}*x^2 + a_{n-1}*x + a_n \f$.
   * 
   * Note the reverse indexation of the coefficients.
   *
   * @note 
   */
  template<class T> 
  class Polynomial
  {
  public:
    /**
     * Default constructor.
     */
    Polynomial(size_t order=0) : mCoefficients(order + 1)
    {}

    /**
     * Constructor to inititalaze any polynomial.
     *
     * @param coefficients Coefficients \f$ a_i, i=0..n \f$. 
     */
    Polynomial(const std::valarray<T>& coefficients) 
      : mCoefficients(coefficients)
    {}

    /**
     * Copy-constructor.
     *
     * @param p Another polynomial.
     */
    Polynomial(const Polynomial& p)
      : mCoefficients(p.mCoefficients)
    {}
      
//     /**
//      * Constructor to inititalaze a linear form polynomial 
//      * \f$ y(x) = a*x + b* \f$.
//      *
//      * @param a a.
//      * @param b b.
//      */
//     Polynomial(const T a, const T b)
//       : mCoefficients(2)
//     {
//       mCoefficients[0] = a;
//       mCoefficients[1] = b;
//     }

//     /**
//      * Constructor to inititalaze a quadratic polynomial 
//      * \f$ y(x) = a*x^2 + b*x + c \f$.
//      *
//      * @param a a.
//      * @param b b.
//      * @param c c.
//      */
//     Polynomial(const T a, const T b, const T c)
//       : mCoefficients(3)
//     {
//       std::cout << mCoefficients.size() << std::endl;
//       mCoefficients[0] = a;
//       mCoefficients[1] = b;
//       mCoefficients[2] = c;
//     }

    /**
     * Valarray does not support operator = for different sized arrays,
     * we have to resize them first.
     *
     */
    void operator=(const Polynomial& p) 
    {
      mCoefficients.resize(p.mCoefficients.size());
      mCoefficients = p.mCoefficients;
    }

//     T getCoefficient(const size_t index) const
//     {
//       return mCoefficients[index];
//     }

    T operator[] (const size_t index) const
    {
      return mCoefficients[index];
    }

    /**
     * Returns the polynomial as a string.
     *
     */
    std::string toString() const
    {
      size_t n = mCoefficients.size();
      std::ostringstream strStream;
      size_t i = 0;
      while (1) {
        strStream << mCoefficients[i];

        // Print power if power > 0
        if (i < n - 1) {
          strStream << "*x**" << (n - i - 1);
        }
        
        // Exit loop after last coefficient
        if (n <= ++i) {
          break;
        }

        strStream << " + ";
      }
      return strStream.str();
    }

  private:
    /// Contains coefficients (beginning with mCoefficients[0] = a_0).
    std::valarray<T> mCoefficients;
  };
}
#endif

