/**
 * Defines Utility functions for std::valarray.
 * 
 * @todo Do proper documentation
 */
#ifndef _VALARRAYUTIL_
#define _VALARRAYUTIL_

#include <string>
// #include <cstdlib>
#include <iostream>
#include <fstream>
#include <ostream>
#include <valarray>

#include <boost/lexical_cast.hpp>

/**
 * Prints a valarray on one line, enclosed by curly braces,
 * e.g., "{ 1 2 3 }".
 * 
 * @param out Output stream
 * @param a Valarray to print
 */ 
template<typename T>
void printValarray(std::ostream& out,
                    const std::valarray<T>& a)
{
  out << '{';
  for (std::size_t i = 0; i < a.size(); ++i)
    out << ' ' << a[i];
  out << " }";
}

/**
 * Prints a valarray, each value on one line
 * 
 * @param out Output stream
 * @param a Valarray to print
 */ 
template<typename T>
void printValarrayLinebyline(std::ostream& out,
                             const std::valarray<T>& a)
{
  for (std::size_t i = 0; i < a.size(); ++i)
    out << a[i] << std::endl;
}

// Print a slice_array, gslice_array, etc., by converting
// to a valarray. Converting a valarray to a valarray is
// wasteful, but harmless for these simple examples.

// The commented out ode interacts badly with piping other things
// e.g. strings...
// template<template<typename> class U, typename T>
// std::ostream& operator<<(std::ostream& out, const U<T>& x)
// {
//   printValarray(out, static_cast<const std::valarray<T> >(x));
//   return out;
// }
template<typename T>
std::ostream& operator<<(std::ostream& out, const std::valarray<T>& x)
{
  printValarray(out, static_cast<const std::valarray<T> >(x));
  return out;
}


/**
 * Computes exp10 for all elements of a valarray.
 *
 * @param in Input valarray.
 *
 * @return Valarray with all elements exp10-ed.
 */
template<typename T>
std::valarray<T> exp10(const std::valarray<T>& in)
{
  std::valarray<T> out(0.0, in.size());
  for (std::size_t i = 0; i < in.size(); ++i){
   out[i] =  exp10(in[i]);
  }
  return out;
}
    

/**
 * Computes the sum of to valarray vectors
 *
 * @param in1 fisrt valarray Vector, in2 second valarray Vector
 *
 * @return valarray Vector with the pairwise sum.
 */
//template<typename T>
//std::vector< std::valarray<T>  > operator+(const std::vector< std::valarray<T> >& in1, const std::vector< std::valarray<T> >& in2)
//{
//	assert(in1.size() == in2.size());
//	std::vector< std::valarray<T> > out(in1.size());
////  std::valarray<T> out(0.0, in.size());
//	for (std::size_t i = 0; i < in1.size(); ++i){
//		out[i] =  in1[i] + in2[i] ;
//	}
//	return out;
//}




/**
 * Gives back a flatend vector when given a vector of valarrays
 *
 * @param in vector of valarrays
 *
 * @return flatened vector
 */
template<typename T>
std::vector<T> valarrayVectorFlat(const std::vector< std::valarray<T> >& in)
{
  std::vector<T> out;
  for (std::size_t i = 0; i < in.size(); ++i){
	  backInsertValarrayElements(in[i], out);
  }
  return out;
}


  
/**
 * Converts a vector of any numerical type to a valarray.
 *
 * @param vec The input vector.
 * @return Valarray with the same values as the input vector.
 *
 * @todo Test
 */
template<class T> 
inline std::valarray<T> convertVectorToValarray(const std::vector<T>& vec) 
{    
  return std::valarray<T>(&vec[0], vec.size());
}

/**
 * Append all elements of an valarray to a vector.
 *
 * @todo Generalize for Back Insertion Sequence instead vector.
 *
 * @param v Valarray with elements to append.
 * @param target Vector to append elements to
 */
template<class T> 
void backInsertValarrayElements(const std::valarray<T> v, std::vector<T>& target)
{
  for (size_t i = 0; i < v.size(); ++i) {
    target.push_back(v[i]);
  }
}


/**
 * Function object to append all elements of an valarray to a vector
 * 
 * @note Does the same as backInsertValarrayElements but as a functor.
 */
template<class T> 
class BackInsertValarrayElements
{
public:
  /// Constructor
  BackInsertValarrayElements(std::vector<T>& target) 
    : mTarget(target) {}

  /// Performs append operation
  void operator() (std::valarray<T>& v)
  {
    backInsertValarrayElements(v, mTarget);
//     for (size_t i = 0; i < v.size(); ++i) {
//       mTarget.push_back(v[i]);
//     }
  }
private:
  std::vector<T>& mTarget;
};

/**
 * Reads stream into a valarray, each value on one line
 * 
 * @param in Input stream
 * @todo Make generic
 */ 
//template<typename T>
inline std::valarray<double> readValarrayLinebyline(std::ifstream& in)
{
  std::vector<double> v;
  std::string currentLine;
  getline(in, currentLine);
  try {
    v.push_back(boost::lexical_cast<double>(currentLine));
  }
  catch(boost::bad_lexical_cast &) {
    // Ignore lines that can't be interpreted
  }
  
  while (!in.eof()) {
    getline(in, currentLine);
    try {
      v.push_back(boost::lexical_cast<double>(currentLine));
    }
    catch(boost::bad_lexical_cast &) {
      // Ignore lines that can't be interpreted.
    }
  }
  std::valarray<double> back = convertVectorToValarray(v);
  return back;
}


/**
 * Multiplies a vectors of valarrays with a scalar 
 * 
 * @param v vector of valarrays
 * @param s scalar multiplier
 * @return New vector v of valarrays, with v[i] = a[i] * s  for all i
 */
template<typename T>
std::vector< std::valarray<T> > operator*(const std::vector< std::valarray<T> >& a, 
                                          double  s)
{
  std::vector< std::valarray<T> > result;
  for (size_t i = 0; i < a.size(); ++i) {
    result.push_back(a[i] * s);
  }
  return result;
}


/**
 * Adds two vectors of valarrays 
 * 
 * @param a First vector of valarrays
 * @param b Second vector of valarrays
 * @return New vector v of valarrays, with v[i] = a[i] + b[i] for all i
 */
template<typename T>
std::vector< std::valarray<T> > operator+(const std::vector< std::valarray<T> >& a, 
                                          const std::vector< std::valarray<T> >& b)
{
  assert(a.size() == b.size());
  std::vector< std::valarray<T> > result;
  for (size_t i = 0; i < a.size(); ++i) {
    result.push_back(a[i] + b[i]);
  }
  return result;
}


/**
 * Adds two vectors of same size and type
 * 
 * @param a First vector 
 * @param b Second vector
 * @return New vector v, with v[i] = a[i] + b[i] for all i
 */
template<typename T>
std::vector<T> operator+(const std::vector<T>& a, 
                                          const std::vector<T>& b)
{
  assert(a.size() == b.size());
  std::vector<T> result;
  for (size_t i = 0; i < a.size(); ++i) {
    result.push_back(a[i] + b[i]);
  }
  return result;
}

  
  /**
   * Subtracts two vectors of same size and type
   * 
   * @param a First vector 
   * @param b Second vector
   * @return New vector v, with v[i] = a[i] - b[i] for all i
   */
  template<typename T>
  std::vector<T> operator-(const std::vector<T>& a, 
                                            const std::vector<T>& b)
  {
    assert(a.size() == b.size());
    std::vector<T> result;
    for (size_t i = 0; i < a.size(); ++i) {
      result.push_back(a[i] - b[i]);
    }
    return result;
  }
  

/**
 * Subtracts two vectors of valarrays 
 * 
 * @param a First vector of valarrays
 * @param b Second vector of valarrays
 * @return New vector v of valarrays, with v[i] = a[i] + b[i] for all i
 */
template<typename T>
std::vector< std::valarray<T> > operator-(const std::vector< std::valarray<T> >& a, 
                                          const std::vector< std::valarray<T> >& b)
{
  assert(a.size() == b.size());
  std::vector< std::valarray<T> > result;
  for (size_t i = 0; i < a.size(); ++i) {
    result.push_back(a[i] - b[i]);
  }
  return result;
}

/**
 * Divides a vector of valarrays by a scalar
 * 
 * @param a First vector of valarrays
 * @param b A scalar
 * @return New vector v of valarrays, with v[i] = a[i] / b for all i
 */
template<typename T>
std::vector< std::valarray<T> > operator/(const std::vector< std::valarray<T> >& a, 
                                          const T b)
{
  std::vector< std::valarray<T> > result;
  for (size_t i = 0; i < a.size(); ++i) {
    result.push_back(a[i] / b);
  }
  return result;
}

/**
 * Subtracts the vector from a scalar 
 * 
 * @param a scalar
 * @param b Second vector of valarrays
 * @return New vector v of valarrays, with v[i] = a - b[i] for all i
 */
template<typename T>
std::vector< std::valarray<T> > operator-(const T a, const std::vector< std::valarray<T> >& b)
{
  std::vector< std::valarray<T> > result;
  for (size_t i = 0; i < b.size(); ++i) {
    result.push_back(a - b[i]);
  }
  return result;
}

/**
 * Applies exp10 to a vector of valarrays 
 * 
 * @param a First vector of valarrays
 * @param b Second vector of valarrays
 * @return New vector v of valarrays, with v[i] = a[i] + b[i] for all i
 */
template<typename T>
std::vector< std::valarray<T> > exp10(const std::vector< std::valarray<T> >& a)
{
  std::vector< std::valarray<T> > result;
  for (size_t i = 0; i < a.size(); ++i) {
    result.push_back(exp10(a[i]));
  }
  return result;
}



#endif
