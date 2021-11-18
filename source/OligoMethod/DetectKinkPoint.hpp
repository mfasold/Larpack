/**
 * @file DetectKinkPoint.hpp Methods to detect the kink point in the hookcurve.
 * @author $Author: mario $ 
 * @author Mario Fasold
 * @date $Date: 2007-04-20 16:55:02 +0200 (Fri, 20 Apr 2007) $
 */
#ifndef _DETECTKINKPOINT_
#define _DETECTKINKPOINT_

#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include "HookcurveAnalyzer.hpp"
#include "Polynomial.hpp"

namespace larrpack {
  class DetectKinkPoint; // Forward declaration

  /// Two polynomials are returned by the least square fitting methods
  typedef std::pair<Polynomial<IntensityType>, Polynomial<IntensityType> > TwoPolynomials;

  /// A function object for the following methods
  typedef boost::function<IntensityPair (IntensityMapping&)> DetectKinkPointFunction;
  // Not all compilers like the declaration above. Alternatively, try:
  // typedef boost::function1<IntensityPair, IntensityMapping& > DetectKinkPointFunction;

  /// Pointer to one of the methods
  /// @note Using a pointer instead of the functor allows us to call
  /// setExportFilename for all methods. As a drawback, the function 
  /// call has the uncommon form (*detectKinkPoint)(graph). 
  typedef boost::shared_ptr<DetectKinkPoint> DetectKinkPointPtr;

  /**
   * Abstract interface for methods to detect a hookcurve's kink 
   * point, that is the point that separates non-specific and specific
   * binding.
   *
   */
  class DetectKinkPoint
  {
  public:
    DetectKinkPoint(const std::string filename = "");
    
    /// Checks if the correction point is valid (within the margins of the graph)
    /// If it is valid that same point is returned, if not it is set to the 0.05*graph.size() probesets average sumlogi
    IntensityPair checkAndCorrectKinkPoint(IntensityPair ip, const IntensityMapping& graph);
    
    /// returns the name of the detect kink point method.
    std::string getMethodName();
    
    /**
     * Excecutes the kink point detection.
     *
     * @param graph Points of a graph.
     * @return Kink point between specific and non-specific binding.
     */
    virtual IntensityPair operator() (const IntensityMapping& graph) = 0;

    virtual void setExportFilename(const std::string filename = ""); 
    
    //virtual void printGeneralGnuplotCommands(const IntensityMapping& graph, const ofstream& gnuplotFile);
    

    /// Virtual classes need virtual destructors
    virtual ~DetectKinkPoint() {};

  protected:
    /// Shall method operator() document detection to a file?
    bool mDoExport;    
    
        
     /// Name of the method.
    std::string mMethodName;

    /// Name of the file to export
    std::string mExportFilename;
  };

  /**
   * @class FitStraightLineAndParabola
   * @brief Computes the kink point by looking for the best fit of a straight 
   * line to the points \f$ 1..i \f$ and a parabola to the points \f$ i..n \f$ 
   * for all \f$ i=1..n \f$.
   * 
   * @see FitStraightLineAndParabola::fitLineAndParabola
   */
  class FitStraightLineAndParabola : public DetectKinkPoint
  {
  public:
    FitStraightLineAndParabola(const std::string filename = "");
    IntensityPair operator() (const IntensityMapping& graph);
    std::string getMethodSpecificGnuplotCommands();
  private:
    std::pair<Polynomial<IntensityType>, Polynomial<IntensityType> > mFittedPolynomials;
    //std::pair<Polynomial<IntensityType>, Polynomial<IntensityType> > getPolynomials();
    static TwoPolynomials fitLineAndParabola(const IntensityMapping& probeData);
    static IntensityType calculateSmallestIntersectionPoint(const Polynomial<IntensityType>& line, const Polynomial<IntensityType>& parabola);
  };

  /**
   * @class FirstDerivativeAnalysis
   * @brief Computes the maximum of the derivative of the digitized hookcurve.
   * The part left from that maximum is fitted by two straight lines. 
   * (See FitStraightLineAndParabola)
   * 
   */
  class FirstDerivativeAnalysis : public DetectKinkPoint
  {
    typedef IntensityMapping::const_iterator GraphIterator;
    public:
      FirstDerivativeAnalysis(const std::string filename = "");
      IntensityPair operator() (const IntensityMapping& graph);
  };

  /**
   * @class FirstDerivativeIntersectionWithZeroAnalysis
   * @brief This is an extensio to the FirstDerivativeAnalysis. The algorithm is the same
   * except not the intersection between two straight lines is used, but the intersection
   * of the second line with (delta) y = 0.C
   * 
   */
  class FirstDerivativeIntersectionWithZeroAnalysis : public DetectKinkPoint
  {
    typedef IntensityMapping::const_iterator GraphIterator;
  public:
    FirstDerivativeIntersectionWithZeroAnalysis(const std::string filename = "");
    IntensityPair operator() (const IntensityMapping& graph);
  };

  
  /**
   * @class SecondDerivativeAnalysis
   * @brief 
   * 
   */
  class SecondDerivativeAnalysis : public DetectKinkPoint
  {
    typedef IntensityMapping::const_iterator GraphIterator;
    public:
      SecondDerivativeAnalysis(const std::string filename = "");
      IntensityPair operator() (const IntensityMapping& graph);
  };
  
  /**
   * @class FitTwoStraightLines
   * @brief Computes the kink point by fitting two lines to the points
   * \f$ i..j \f$ and \f$ j..k \f$, respectively, for all \f$ j=i..k \f$. 
   *
   * As a first step, we 
   * divide the hookcurve into c intervals and calculate the ascent of 
   * each (using a fitted straight line). Let I be the interval with the
   * steepest ascent, then i is the first point of interval I-1 and k the 
   * last point in I.
   *
   */
  class FitTwoStraightLines : public DetectKinkPoint
  {    
    typedef IntensityMapping::const_iterator GraphIterator;

  public:
    FitTwoStraightLines(const size_t intervalCount, const IntensityType robustnessFactor = 0.6, 
                        const std::string filename = "");
    IntensityPair operator() (const IntensityMapping& graph);

  private:
    static std::vector< std::pair<GraphIterator, IntensityType> > getIntervalsAndAscents(const GraphIterator& leftPoint, 
                                                                                          const GraphIterator& rightPoint,
                                                                                          const size_t intervalCount);

    static IntensityType getPointBeforeSteepestAscent(const IntensityMapping& graphPoints, 
                                                       size_t intervalCount);

  public:
    static TwoPolynomials fitTwoLines(const GraphIterator& leftPoint, 
                                      const GraphIterator& rightPoint);

  private:
    /// Number of intervals the hookcurve is devided into
    size_t mIntervalCount;

    /// The algorithm uses left interval as maximum as long it's ascent >= rubustnessFactor*currentAscent
    IntensityType mRobustnessFactor;
  };

  /**
   * @class FitDigitizedPlot
   * Runs one of the other detection methods but digitizes the plot 
   * to a lower number of points before. 
   *
   * Implements the decorator design pattern.
   *
   */
  class FitDigitizedPlot : public DetectKinkPoint
  {
  public:
    FitDigitizedPlot(DetectKinkPointPtr detectKinkPointFunctor, 
                     const size_t digitizeIntervals);
    IntensityPair operator() (const IntensityMapping& graph);

    void setExportFilename(const std::string filename = ""); 
  private:
    /// The function that detects the kink point, e.g. FitTwoLines
    DetectKinkPointPtr mDetectKinkPoint;

    /// Number of points the plot is digitized to
    size_t mDigitizeIntervals;
  };

}

#endif
