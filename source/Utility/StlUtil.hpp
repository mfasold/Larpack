/**
 * @file StlUtil.hpp Defines additional functions for STL types
 * @author $Author$ 
 * @author Mario Fasold
 * @date $Date$
 */
#ifndef _STLUTIL_
#define _STLUTIL_

#include <vector>
#include <cassert>
#include <functional>
#include <iostream>
#include <fstream>
#include <boost/function.hpp>

namespace larrpack {
  /**
   * Returns a slice of a vector. Using this function for vectors with large objects or with  large slices
   * is not recommended, since the vector slice objects are copied at least twice (frist, to create 
   * a temporal slice vector, second, during return).
   *
   *
   * @note This function is located in the header file to avoid linking errors with function
   *       templates.
   * @see  http://www.parashift.com/c++-faq-lite/templates.html#faq-35.13
   *
   * @note I am aware that stl::vector is not designed for this kind of access. Iterators are a much 
   * better way to do this. However, since some routins heavily work with vector indices, this function
   * is useful until a better design is available.
   *
   * @param entireVector Input vector from which to cut the slice out
   * @param firstIndex First element of the slice 
   * @param lastIndex Last element of the slice 
   *
   * @return A new vector with elements firstIndex..lastIndex from the input vector.
   * 
   */
  template<class T> std::vector<T> getVectorSlice(const std::vector<T>& entireVector, size_t firstIndex, size_t lastIndex)
  {
    // Check bounds
    assert(firstIndex >= 0 && lastIndex >= 0);
    assert(lastIndex - firstIndex >= 0);
    assert(entireVector.size() > lastIndex);

    // Copy splice to temporary vector and return it
    std::vector<T> tmpVector(lastIndex - firstIndex + 1);
    for (size_t i = 0; i <= lastIndex - firstIndex; ++i) {
      tmpVector[i] = entireVector[firstIndex + i];
    }
    return tmpVector;
  }


  /**
   * Returns a new vector containing only elements with index i from vector a, 
   * if and only if a function f is true 
   * for element with index i in b. That is f(b[i]) is true.
   *
   * @note We could rewrite this function using a FunctionObject and the STL erase/remove
   *       functions.
   *
   * @todo Make it more general.
   * @todo Rename function.
   *
   * @param a Vector to filter elements from
   * @param b Testing vector
   * @param f Testing function
   */
  typedef boost::function<bool (double)> DoublePredicate;
  template<class T> std::vector<T> filterParallel(const std::vector<T>& a, const std::vector<double>& b, DoublePredicate f)
  {
    assert(a.size() == b.size());
    std::vector<T> filtered;
    for (size_t i = 0; i < b.size(); ++i) {
      if (f(b[i])) {
        filtered.push_back(a[i]);
      }
    }
    return filtered;
  }


  /**
   * Template functor that returns false for Xth element
   * and true otherwise. That can be used to delete 
   * everything but the  Xth element of a vector.
   *
   * @note It seems as the unary negator only works for
   *       const operator(), which is not possible here. 
   */
  template<class T> class NotEveryXthElement : public std::unary_function<T, bool>
  {
  public:
    NotEveryXthElement(const int x) : mX(x) { counter = 0; }
    bool operator() (const T& t) { return !(++counter % mX == 0); }
  private:
    int mX;
    long counter;
  };
  
  /**
   * Overrides << to send a pair to any stream.
   *
   * @param out Target stream.
   * @param x The pair to be sent to the stream.
   */
  template<typename T>
  std::ostream& operator<<(std::ostream& out, const std::pair<T,T>& x)
  {
    out << x.first << "\t" << x.second;
    return out;
  }

  /**
   * Templatized functor to print a std::pair. Can be used to print
   * a list of pairs using for_each. It takes advantage of the
   * type conversion ability of streams.
   * 
   */
  template<class PairType> 
  class PairPrinter
  {
  public:
    /**
     * Constructor
     * @param out An output stream.
     */ 
    PairPrinter(std::ostream& out = std::cout) : mOut(out) {}
    void operator() (const PairType& p) { mOut << p.first << "\t" << p.second << std::endl; }
  protected:
    /// Internal reference to the output stream
    std::ostream& mOut;
  };

  /**
   * Templated function that prints a vector of pairs to a file.
   *
   * @param pairList List of pairs
   * @param filename of the output file.
   * @param header Optional commentary first line 
   */
  template<class vectorType> 
  void printVectorPairsToFile(const vectorType& pairList, const std::string filename, const std::string header = "") {
    std::ofstream dataFile;
    dataFile.open(filename.c_str());

    // Write the data using PairPrinter
    if (header.size() > 0) {
      dataFile << header << std::endl; 
    }
    for_each(pairList.begin(), pairList.end(), PairPrinter<typename vectorType::value_type>(dataFile));

    dataFile.close();
  }

  /**
   * Templated functor to add a pair a to pair b, elementwise. In other 
   * words, it returns a pair c with c.first = a.first + b.first
   * and c.second = a.second + b.second.
   * Usage example: PairType c = AddToPair<PairType>(b)(a);
   *
   */
  template<class PairType> 
  class AddToPair : public std::unary_function<PairType, PairType>
  {
  public:
    /**
     * Constructor
     *
     * @param b Second pair.
     */ 
    AddToPair(const PairType& b) : mB(b) {}
    /**
     * Performs addition
     *
     * @param a First pair.
     * @return a + b.
     */ 
    PairType operator() (const PairType& a) 
    { return PairType(a.first + mB.first, a.second + mB.second); }
  protected:
    /// Internal reference to b
    const PairType mB;
  };
  
  // @Note: I'd like to use a templated function instead of the above object.
  // However, i don"t know how to get a pointer to a templatized function.
  // Running "bind2nd(ptr_fun(addToPair),b)" gives
  // error: no matching function for call to ‘ptr_fun(<unresolved overloaded function type>)’
  //   template<class PairType> 
  //   PairType addToPair(PairType a, PairType b) 
  //   {
  //     return PairType(a.first + b.first, a.second + b.second);
  //   }

}

#endif
