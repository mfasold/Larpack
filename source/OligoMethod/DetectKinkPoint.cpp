/**
 * @file DetectKinkPoint.cpp Defines methods to detect the kink point in the hookcurve.
 * @author $Author: mario $ 
 * @author Mario Fasold
 * @date $Date: 2007-04-20 16:55:02 +0200 (Fri, 20 Apr 2007) $
 */

#include "DataCollector.hpp"
#include "DetectKinkPoint.hpp"
#include "MathUtil.hpp"
#include "StlUtil.hpp"
#include "ValarrayUtil.hpp"
#include "StringUtil.hpp"


using namespace std;
using namespace larrpack;
using namespace MathUtil;

/**
 * Constructor. Lets you optionally provide a filename to enable
 * debugging.
 *
 * @param filename File to export results of the detection.
 */
DetectKinkPoint::DetectKinkPoint(const std::string filename)
 //: mDoExport(!filename.empty()), mExportFilename(filename)
{
  setExportFilename(filename);
}

/**
 * Returns a secription of the detect kink point method.
 * 
 */
string DetectKinkPoint::getMethodName() 
{
  return mMethodName;
}

/**
 * Sets the name of the export file. If no name is set, nothing
 * will be exported.
 *
 * @param filename File to export results of the detection.
 */
void DetectKinkPoint::setExportFilename(const std::string filename)
{
  if (!filename.empty()) {
    mExportFilename = filename;
    mDoExport = true;
  }
  else {
    mDoExport = false;
  }
}

/**
 * Checks if the calculated intersection point is within the margins and corrects it some plausible value
 * if it is outside.
 * 
 * @param ip Calculated Intersection Point
 * @param graph The probesets sorted by SumLogI
 * 
 */
IntensityPair DetectKinkPoint::checkAndCorrectKinkPoint(IntensityPair ip, const IntensityMapping& graph)
{
  //size_t smallestIndex = graph()
  size_t index = (size_t) (0.05 * graph.size());
  //size_t index2 = (size_t) (0.97 * graph.size());
  IntensityType marginLeft   = graph[index].first;//graph.begin()->first;
  //IntensityType marginRight  = graph[index2].first;//(graph.end()-1)->first;
//   cout << "smallest : " << smallestSumLogI << endl;
//   cout << "biggest : " << biggestSumLogI << endl;
  if (ip.first < marginLeft) {
    //cout << "Intersection: " << ip.first << " margin: " << marginLeft << endl; 
    return graph[index];
  }
  else {
    return ip;
  }
}



/**
 * Constructor. Lets you optionally provide a filename to enable
 * debugging.
 *
 * @param filename File to export results of the detection.
 */
FitStraightLineAndParabola::FitStraightLineAndParabola(const std::string filename)
  : DetectKinkPoint(filename)
{
  mMethodName = "Line and parabola fit";
}

/**
 * Executes the kink point detection.
 *
 * @param graph The (x,y) pairs in the graph sorted ascending in x.
 * @return Kink point between specific and non-specific binding.
 */
IntensityPair FitStraightLineAndParabola::operator() (const IntensityMapping& graph)
{
//   std::pair<Polynomial<IntensityType>, Polynomial<IntensityType> > fittedPolynomials = fitLineAndParabola(graph);
  mFittedPolynomials = fitLineAndParabola(graph);
  
  
  IntensityType intersectionPoint = calculateSmallestIntersectionPoint(mFittedPolynomials.first,
                                                                        mFittedPolynomials.second);
  // The corresponding y value is calculated by the linear function of the fitting.
  IntensityType intersectionPointYvalue = mFittedPolynomials.first[0] * intersectionPoint + 
    mFittedPolynomials.first[1]; // y = a*x + b
  
  if (mDoExport) {
    // Export results of hookcurve fitting      
    // This plots the (usually) digitized hookcurve to file:
    printVectorPairsToFile(graph, mExportFilename); 

    string filenameBase = stringutil::splitString(mExportFilename, ".")[0];    
    string additionalCommands = "";

    DataCollector::instance().insert(filenameBase+"-intersectionPointXvalue" , intersectionPoint);
    DataCollector::instance().insert(filenameBase+"-intersectionPointYvalue" , intersectionPointYvalue);
    
    additionalCommands += ", " + mFittedPolynomials.first.toString() + " title 'Line Fit' lt 5 lw 2";
    additionalCommands += ", " + mFittedPolynomials.second.toString() + " title '' lt 5 lw 2";
    DataCollector::instance().insert(filenameBase+"-plottingCommands" , additionalCommands);
  }
  IntensityPair ip = IntensityPair(intersectionPoint, intersectionPointYvalue);
  ip = checkAndCorrectKinkPoint(ip, graph);
  return ip;
}


/**
 * Calculates the best least squares fitting of a linear and a quadratic function,
 * such that the linear function fits datapoints \f$ 0..x_g \f$ and the quadratic 
 * function fits points \f$ x_g +1..n \f$. The dataset are x/y pairs sorted by 
 * the x value. The routine uses SVD (MathUtil::solveEquationsSvd) for solving 
 * the equations. 
 *
 * @param probeData List of points sorted by first member (x-value).
 *
 * @todo Remove debug code.
 *
 * @return Parameters for both optimal fitting functions. 
 */
TwoPolynomials  FitStraightLineAndParabola::fitLineAndParabola(const IntensityMapping& probeData)
{
  // Allocate matrices and vectors for linear (2x2,2x1) and quadratic (3x3, 3x1)
  // least squares matrix and initialize with zeros (calloc instead alloc)
  math::Matrix<IntensityType> deltaLinear(2,2, 0.0);
  IntensityArray thetaLinear(0.0, 2);
  math::Matrix<IntensityType> deltaQuadratic(3,3, 0.0);
  IntensityArray thetaQuadratic(0.0, 3);
 
  // Initialize quadratic fit for x_g = 0
  // (-> quadratic fit from 0 to n)
  for (size_t i = 0; i < probeData.size(); ++i) {
    IntensityType x = probeData[i].first;
    IntensityType y = probeData[i].second;
    deltaQuadratic[0][0] +=  (x*x*x*x);
    deltaQuadratic[0][1] +=  (x*x*x);
    deltaQuadratic[0][2] +=  (x*x);

    deltaQuadratic[1][0] += (x*x*x);
    deltaQuadratic[1][1] += (x*x);
    deltaQuadratic[1][2] += (x);

    deltaQuadratic[2][0] += (x*x);
    deltaQuadratic[2][1] += (x);
    deltaQuadratic[2][2] += (1);

    thetaQuadratic[0] += (y*x*x);
    thetaQuadratic[1] += (y*x);
    thetaQuadratic[2] += (y);
  }

  // Set minimal error value to INF and index to -1
  IntensityType minimalTotalError = numeric_limits<IntensityType>::infinity(); // GSL_POSINF;  
 
  // Remember best linear and quadratic function
  Polynomial<IntensityType> bestLinearFit(1);
  Polynomial<IntensityType> bestQuadraticFit(2);
 
  //  For all x_g from 0 to n - 2
  for (size_t i = 0; i <= probeData.size() - 2; ++i) {
    // x,y are the values of current x_g
    IntensityType x = probeData[i].first;
    IntensityType y = probeData[i].second;

    // Add values for linear fit (see Preibisch, S.42)
    deltaLinear[ 0][ 0] += (x*x);
    deltaLinear[ 0][ 1] += (x);
    deltaLinear[ 1][ 0] += (x);
    deltaLinear[ 1][ 1] += (1);

    thetaLinear[ 0] += (y*x);
    thetaLinear[ 1] += (y);

    // Substract values for quadratic fit (see Preibisch, S.42)
    deltaQuadratic[0][ 0] -=  (x*x*x*x);
    deltaQuadratic[0][ 1] -=  (x*x*x);
    deltaQuadratic[0][ 2] -=  (x*x);

    deltaQuadratic[ 1][ 0] -= (x*x*x);
    deltaQuadratic[ 1][ 1] -= (x*x);
    deltaQuadratic[ 1][ 2] -= (x);

    deltaQuadratic[ 2][ 0] -= (x*x);
    deltaQuadratic[ 2][ 1] -= (x);
    deltaQuadratic[ 2][ 2] -= (1);

    thetaQuadratic[ 0] -= (y*x*x);
    thetaQuadratic[ 1] -= (y*x);
    thetaQuadratic[ 2] -= (y);

    // Solve linear equations
    IntensityArray straightLineFit = MathUtil::solveEquationsSvd(deltaLinear, thetaLinear, 0);
    IntensityArray parabolaFit = MathUtil::solveEquationsSvd(deltaQuadratic, thetaQuadratic, 0);
//     cout << "deltaLinear: " << deltaLinear << endl;
//     cout << "thetaLinear: " << thetaLinear << endl;
//     cout << "deltaQuadratic: " << deltaQuadratic << endl;
//     cout << "thetaQuadratic: " << thetaQuadratic << endl;
//     IntensityArray straightLineFit = MathUtil::solveEquationsLU(deltaLinear, thetaLinear);
//     IntensityArray parabolaFit = MathUtil::solveEquationsLU(deltaQuadratic, thetaQuadratic);

    // Calculate error for linear fit
    IntensityType errorLinear = 0;
    for (size_t j = 0; j <= i; ++j) {
      IntensityType dataX = probeData[j].first;
      IntensityType dataY = straightLineFit[0]*dataX + straightLineFit[1];  // y = a*x + b
      errorLinear += (dataY - probeData[j].second)*(dataY - probeData[j].second);
    }

    // Calculate error for quadratic fit
    IntensityType errorQuadratic = 0;
    for (size_t j = i + 1; j < probeData.size(); ++j) {
      IntensityType dataX = probeData[j].first;
      IntensityType dataY = parabolaFit[0]*dataX*dataX + parabolaFit[1]*dataX  + parabolaFit[2]; // y = a*x² + bx + c
      errorQuadratic += (dataY - probeData[j].second)*(dataY - probeData[j].second);
    }

    // Calculate weighted total error
    IntensityType totalError = errorLinear/(i + 1) + errorQuadratic/(probeData.size() - i - 1);
    
    // If current total error smaller than current minimal error, reset minimal and index
//     cout << "i: " << i << endl;
//     cout << "errorLinear: " << errorLinear << endl;
//     cout << "errorQuadratic: " << errorQuadratic << endl;
//     cout << "straightLineFit: " << straightLineFit << endl;
//     cout << "parabolaFit: " << parabolaFit << endl;
    if (totalError < minimalTotalError) {
      minimalTotalError = totalError;
      
      // Remember minimal values
      bestLinearFit = straightLineFit;
      bestQuadraticFit = parabolaFit; // = Polynomial<IntensityType>(parabolaFit);
//       cout << "minimalTotalError: " << minimalTotalError << endl;
    }
  }

  return make_pair(bestLinearFit, bestQuadraticFit);
} 

/**
 * Calculates the leftmost intersections point of the linear and
 * the quadratic function, if one exists. In other words, it looks
 * for the smallest x, such that \f$ mx + n = ax^2 + bx + c \f$. That is
 * \f$ x_{1/2} = - \frac{p}{2} \pm \sqrt{\frac{p^2}{4} - q} \f$ 
 * with \f$ p = \frac{b-m}{a}, q = \frac{c-n}{a} \f$.
 * Since p~q~a in the order of magnitude and they tend to be very
 * small, one might use the alternative form 
 * \f$ x_{1/2} = - \frac{p}{2a} \pm \frac{1}{a} \sqrt{\frac{p^2}{4} - q*a} \f$
 * with \f$ p = b-m, q = c-n \f$
 * where we subtract similar values under the squareroot.
 *
 * @note Numeric behaviour untested.
 * 
 * @param line Line polynomial.
 * @param parabola Parabola polynomial.
 *
 * @return Leftmost intersection point
 */
IntensityType FitStraightLineAndParabola::calculateSmallestIntersectionPoint(const Polynomial<IntensityType>& line, const Polynomial<IntensityType>& parabola)
{
  IntensityType p = parabola[1] - line[0];
  IntensityType q = parabola[2] - line[1];
  IntensityType a = parabola[0];

  // double radicand = p*p/4 - a*q;
  IntensityType radicand = (p/a)*(p/a)/4 - (q/a);
  
  // Return NaN if no intersection point can be calculated
  if (radicand < 0 || a == 0) {
    return MathUtil::getNaN();
  }
  
  // Return smallest x (this includes the case radicant == 0)
  // return (-p/(2*a) - sqrt(radicand)) / a;
  return (-p/(2*a) - sqrt(radicand));
}



/**
 * Constructor. Lets you optionally provide a filename to enable
 * debugging.
 * 
 * @param filename File to export results of the detection
 * 
 * 
 */
FirstDerivativeAnalysis::FirstDerivativeAnalysis(const std::string filename)
  : DetectKinkPoint(filename)
{
  mMethodName = "Two Line Fit with help of first derivative";
}


/**
 * Executes the kink point detection.
 *
 * @param graph The (x,y)-pairs in the graph sorted ascending in x.
 * @return Kink point between specific and non-specific binding.
 */
IntensityPair FirstDerivativeAnalysis::operator() (const IntensityMapping& graph)
{
  assert(graph.size() > 0);
//  assert(mIntervalCount > 0);

  IntensityMapping derivationPrime;
//   IntensityMapping derivationSecond;
  
  // Get interval borders and ascents
  GraphIterator leftPoint  = graph.begin();
 // vector< pair<GraphIterator, IntensityType> > intervalsAndAscends = getIntervalsAndAscents(leftPoint, graph.end() - 1, mIntervalCount);
  
  size_t slidingWindowSize = 7;
  size_t shift = slidingWindowSize/2;
//   vector< pair<IntensityType, GraphIterator> > maxSlope;
  GraphIterator maximal = graph.begin();
  IntensityType maximalSlope = -50.0; // Some small initial value
  bool calculatedDerivative = true;
  if (graph.size() < slidingWindowSize + 2*shift) { // If the graph is too small we can't calculate a robust derivative, so we better don't
    maximal = graph.end(); // This causes the algorithm to fit two line onto the whole range of the graph.
    calculatedDerivative = false;
  }
  else {
    for (GraphIterator currentPoint = graph.begin(); currentPoint + slidingWindowSize < graph.end(); ++currentPoint) {
      Polynomial<IntensityType> currentLine = MathUtil::fitStraightLine(currentPoint, currentPoint + slidingWindowSize);  
      derivationPrime.push_back(IntensityPair((currentPoint+shift)->first, currentLine[0]));
      if (currentLine[0] > maximalSlope) {
        maximalSlope = currentLine[0];
        maximal = currentPoint + 2*shift; // End of the sliding window in focus (We want to have the point somewhat right of the maximum
      }
  //     maxSlope.push_back(pair<IntensityType, GraphIterator>(currentLine[0], (currentPoint + 4)));
  //     cout << (currentPoint + 4)->first << "\t" << currentLine[0] << endl;
    }
  }
  TwoPolynomials fittedLines = FitTwoStraightLines::fitTwoLines(leftPoint, maximal);
  
  
  // Compute intersection
  IntensityType intersectionPoint = (fittedLines.first[1] - fittedLines.second[1]) / 
      (fittedLines.second[0] - fittedLines.first[0]);
  IntensityType intersectionPointYvalue = fittedLines.first[0]*intersectionPoint + fittedLines.first[1]; // y = a*x + b
  
  
    // Export results
  if (mDoExport) {
    printVectorPairsToFile(graph, mExportFilename); 
    string filenameBase = stringutil::splitString(mExportFilename, ".")[0];
    string derivativeFileName = filenameBase+".firstDerivative";
    if (calculatedDerivative) printVectorPairsToFile(derivationPrime, derivativeFileName);
    
    string additionalCommands = "";
    

    DataCollector::instance().insert(filenameBase+"-intersectionPointXvalue" , intersectionPoint);
    DataCollector::instance().insert(filenameBase+"-intersectionPointYvalue" , intersectionPointYvalue);
    
    additionalCommands += ", " + fittedLines.first.toString()  + " title 'Line Fit' lt 5 lw 2 "; 
    additionalCommands += ", " + fittedLines.second.toString() + " title '' lt 5 lw 2";
    
    //// @Debug prints a gnuplot file that diretly can produce a graph of the data and the fitted polynomials
    string commandsFileName = filenameBase+".gnuplot";
    ofstream gnuplotFile;
    gnuplotFile.open(commandsFileName.c_str());
    string additionalCommands2 = "set term png large size 800,600\n set output \""+filenameBase+".png\"\nset grid\nshow grid";
    additionalCommands2 += "\nplot \""+ mExportFilename+"\"";
    additionalCommands2 += ", " + fittedLines.first.toString()  + " title 'Line Fit' lt 5 lw 2 "; 
    additionalCommands2 += ", " + fittedLines.second.toString() + " title '' lt 5 lw 2";
    if (calculatedDerivative)  additionalCommands2 += ", \"" + derivativeFileName + "\" using 1:2 w l lw 0.3 lt 7 title '" +  "1st derivative" + "'";
    gnuplotFile << additionalCommands2;
    gnuplotFile.close();
    /////
    
    if (calculatedDerivative)  additionalCommands += ", \"" + derivativeFileName + "\" using 1:2 w l lw 0.3 lt 7 title '" +  "1st derivative" + "'";
    
    DataCollector::instance().insert(filenameBase+"-plottingCommands" , additionalCommands);
    
  }
  IntensityPair ip = IntensityPair(intersectionPoint, intersectionPointYvalue);
  ip = checkAndCorrectKinkPoint(ip, graph);
  return ip;
}


/**
 * Constructor. Lets you optionally provide a filename to enable
 * debugging.
 * 
 * @param filename File to export results of the detection
 * 
 * 
 */
FirstDerivativeIntersectionWithZeroAnalysis::FirstDerivativeIntersectionWithZeroAnalysis(const std::string filename)
  : DetectKinkPoint(filename)
{
  mMethodName = "Two Line Fit with help of first derivative, using intersecion with y=0.";
}


/**
 * Executes the kink point detection.
 *
 * @param graph The (x,y)-pairs in the graph sorted ascending in x.
 * @return Kink point between specific and non-specific binding.
 */
IntensityPair FirstDerivativeIntersectionWithZeroAnalysis::operator() (const IntensityMapping& graph)
{
  assert(graph.size() > 0);
  //  assert(mIntervalCount > 0);

  IntensityMapping derivationPrime;
  //   IntensityMapping derivationSecond;
  
  // Get interval borders and ascents
  GraphIterator leftPoint  = graph.begin();
  // vector< pair<GraphIterator, IntensityType> > intervalsAndAscends = getIntervalsAndAscents(leftPoint, graph.end() - 1, mIntervalCount);
  
  size_t slidingWindowSize = 7;
  size_t shift = slidingWindowSize/2;
  //   vector< pair<IntensityType, GraphIterator> > maxSlope;
  GraphIterator maximal = graph.begin();
  IntensityType maximalSlope = -50.0; // Some small initial value
  bool calculatedDerivative = true;
  if (graph.size() < slidingWindowSize + 2*shift) { // If the graph is too small we can't calculate a robust derivative, so we better don't
    maximal = graph.end(); // This causes the algorithm to fit two line onto the whole range of the graph.
    calculatedDerivative = false;
  }
  else {
    for (GraphIterator currentPoint = graph.begin(); currentPoint + slidingWindowSize < graph.end(); ++currentPoint) {
      Polynomial<IntensityType> currentLine = MathUtil::fitStraightLine(currentPoint, currentPoint + slidingWindowSize);  
      derivationPrime.push_back(IntensityPair((currentPoint+shift)->first, currentLine[0]));
      if (currentLine[0] > maximalSlope) {
        maximalSlope = currentLine[0];
        maximal = currentPoint + 2*shift; // End of the sliding window in focus (We want to have the point somewhat right of the maximum
      }
      //     maxSlope.push_back(pair<IntensityType, GraphIterator>(currentLine[0], (currentPoint + 4)));
      //     cout << (currentPoint + 4)->first << "\t" << currentLine[0] << endl;
    }
  }
  TwoPolynomials fittedLines = FitTwoStraightLines::fitTwoLines(leftPoint, maximal);
  
  
  // Compute intersection of second line with y=0
  IntensityType intersectionPoint = - fittedLines.second[1] / fittedLines.second[0]; // y != 0 -> x = -b/a
  IntensityType intersectionPointYvalue = 0; 
  
  
  // Export results
  if (mDoExport) {
    printVectorPairsToFile(graph, mExportFilename); 
    string filenameBase = stringutil::splitString(mExportFilename, ".")[0];
    string derivativeFileName = filenameBase+".firstDerivative";
    if (calculatedDerivative) printVectorPairsToFile(derivationPrime, derivativeFileName);
    
    string additionalCommands = "";
    

    DataCollector::instance().insert(filenameBase+"-intersectionPointXvalue" , intersectionPoint);
    DataCollector::instance().insert(filenameBase+"-intersectionPointYvalue" , intersectionPointYvalue);
    
    additionalCommands += ", " + fittedLines.first.toString()  + " title 'Line Fit' lt 5 lw 2 "; 
    additionalCommands += ", " + fittedLines.second.toString() + " title '' lt 5 lw 2";
    
    //// @Debug prints a gnuplot file that diretly can produce a graph of the data and the fitted polynomials
    string commandsFileName = filenameBase+".gnuplot";
    ofstream gnuplotFile;
    gnuplotFile.open(commandsFileName.c_str());
    string additionalCommands2 = "set term png large size 800,600\n set output \""+filenameBase+".png\"\nset grid\nshow grid";
    additionalCommands2 += "\nplot \""+ mExportFilename+"\"";
    additionalCommands2 += ", " + fittedLines.first.toString()  + " title 'Line Fit' lt 5 lw 2 "; 
    additionalCommands2 += ", " + fittedLines.second.toString() + " title '' lt 5 lw 2";
    if (calculatedDerivative)  additionalCommands2 += ", \"" + derivativeFileName + "\" using 1:2 w l lw 0.3 lt 7 title '" +  "1st derivative" + "'";
    gnuplotFile << additionalCommands2;
    gnuplotFile.close();
    /////
    
    if (calculatedDerivative)  additionalCommands += ", \"" + derivativeFileName + "\" using 1:2 w l lw 0.3 lt 7 title '" +  "1st derivative" + "'";
    
    DataCollector::instance().insert(filenameBase+"-plottingCommands" , additionalCommands);
    
  }
  IntensityPair ip = IntensityPair(intersectionPoint, intersectionPointYvalue);
  ip = checkAndCorrectKinkPoint(ip, graph);
  return ip;
}



/**
 * Constructor. Lets you optionally provide a filename to enable
 * debugging.
 * 
 * @param filename File to export results of the detection
 * 
 * 
 */
SecondDerivativeAnalysis::SecondDerivativeAnalysis(const std::string filename)
  : DetectKinkPoint(filename)
{
  mMethodName = "Intersection point is at the global maximum of 2nd derivative";
}


/**
 * Executes the kink point detection.
 *
 * @param graph The (x,y)-pairs in the graph sorted ascending in x.
 * @return Kink point between specific and non-specific binding.
 */
IntensityPair SecondDerivativeAnalysis::operator() (const IntensityMapping& graph)
{
  assert(graph.size() > 0);
//  assert(mIntervalCount > 0);

  IntensityMapping derivationPrime;
  IntensityMapping derivationSecond;
  
  // Get interval borders and ascents
  GraphIterator leftPoint  = graph.begin();
 // vector< pair<GraphIterator, IntensityType> > intervalsAndAscends = getIntervalsAndAscents(leftPoint, graph.end() - 1, mIntervalCount);
  
  size_t slidingWindowSize = 5;
  size_t shift = slidingWindowSize/2;
//   vector< pair<IntensityType, GraphIterator> > maxSlope;
  GraphIterator maximal = graph.begin();
  IntensityType maximalSlope = -50.0;
  for (GraphIterator currentPoint = graph.begin(); currentPoint + slidingWindowSize < graph.end(); ++currentPoint) {
    Polynomial<IntensityType> currentLine = MathUtil::fitStraightLine(currentPoint, currentPoint + slidingWindowSize);  
    derivationPrime.push_back(IntensityPair((currentPoint+shift)->first, currentLine[0]));
//     if (currentLine[0] > maximalSlope) {
//       maximalSlope = currentLine[0];
//       maximal = currentPoint;
//     }
//     maxSlope.push_back(pair<IntensityType, GraphIterator>(currentLine[0], (currentPoint + 4)));
//     cout << (currentPoint + 4)->first << "\t" << currentLine[0] << endl;
  }
  for (GraphIterator currentPoint = derivationPrime.begin(); currentPoint + slidingWindowSize < derivationPrime.end(); ++currentPoint) {
    Polynomial<IntensityType> currentLine = MathUtil::fitStraightLine(currentPoint, currentPoint + slidingWindowSize);  
    derivationSecond.push_back(IntensityPair((currentPoint+shift)->first, currentLine[0]));
    if (currentLine[0] > maximalSlope) {
      maximalSlope = currentLine[0];
      maximal = currentPoint;
    }
//     maxSlope.push_back(pair<IntensityType, GraphIterator>(currentLine[0], (currentPoint + 4)));
//     cout << (currentPoint + 4)->first << "\t" << currentLine[0] << endl;
  }
  
  
//   pair<IntensityType, GraphIterator> maxi = max(maxSlope);

  // Get find leftmost robust maximum Index 
  // @todo Clean this up.
//   size_t maxIndex = intervalsAndAscends.size() - 1;
//   for (int i = intervalsAndAscends.size() - 2; i >= 0; --i) {
//     if (atan(intervalsAndAscends[i].second) > mRobustnessFactor*atan(intervalsAndAscends[maxIndex].second)) {
//       maxIndex = i;
//     }
//   } 

  // Set left index  
//   size_t middleIntervalIndex = 0;
//   if (maxIndex >= 2) { // Left point does not change iff the steepest ascent is within the first two intervals
//     middleIntervalIndex = maxIndex-1; 
//     leftPoint = intervalsAndAscends[middleIntervalIndex-1].first;
//   } 
  //   
//   assert(middleIntervalIndex + 1 < intervalsAndAscends.size()); 

  // Fit two straight lines
//   Polynomial<IntensityType> leftIntervalLine = MathUtil::fitStraightLine(leftPoint, 
//                                                                         intervalsAndAscends[middleIntervalIndex].first);
//   Polynomial<IntensityType> rightIntervalLine = MathUtil::fitStraightLine(intervalsAndAscends[middleIntervalIndex].first, 
//                                                                          intervalsAndAscends[middleIntervalIndex + 1].first);

  
  
//  TwoPolynomials fittedLines = fitTwoLines(leftPoint, maximal);

  // Compute intersection
  IntensityType intersectionPoint = maximal->first;/*(fittedLines.first[1] - fittedLines.second[1]) / 
      (fittedLines.second[0] - fittedLines.first[0]);
  IntensityType intersectionPointYvalue = fittedLines.first[0]*intersectionPoint + fittedLines.first[1]; // y = a*x + b*/
  IntensityType intersectionPointYvalue = 0.0;
//   cout << "intersectionPoint: " << intersectionPoint << endl;


  // Export results
  if (mDoExport) {
    // Export results of hookcurve fitting      
    printVectorPairsToFile(graph, mExportFilename); 
    string filenameBase = stringutil::splitString(mExportFilename, ".")[0];
    string firstDerivativeFileName = filenameBase+".firstDerivative";
    string secondDerivativeFileName = filenameBase+".secondDerivative";
    printVectorPairsToFile(derivationPrime, firstDerivativeFileName); 
    printVectorPairsToFile(derivationSecond, secondDerivativeFileName); 
    
    string additionalCommands = "";

    DataCollector::instance().insert(filenameBase+"-intersectionPointXvalue" , intersectionPoint);
    DataCollector::instance().insert(filenameBase+"-intersectionPointYvalue" , intersectionPointYvalue);
    
    additionalCommands += ", \"" + firstDerivativeFileName + "\" using 1:2 w l lt 6 lw 0.5 title '" +  "1st derivative" + "'";
    additionalCommands += ", \"" + secondDerivativeFileName + "\" using 1:2 w l lt 7 lw 0.5 title '" + "2nd derivative" + "'";
    DataCollector::instance().insert(filenameBase+"-plottingCommands" , additionalCommands);
  }
  IntensityPair ip = IntensityPair(intersectionPoint, intersectionPointYvalue);
  ip = checkAndCorrectKinkPoint(ip, graph);
  return ip;
}



/**
 * Constructor. Lets you optionally provide a filename to enable
 * debugging.
 *
 * @param intervalCount Number of intervals the hookcurve shall be 
 *  devided into.
 * @param robustnessFactor Parameter that controls the robust search
 *  for the *leftmost* interval i with steepest ascent \f$ \alpha_i \f$. 
 *  That is
 *  \f$ \min_{i} \alpha_i \ge robustnessFactor * \alpha_{max} \f$.
 * @param filename File to export results of the detection.
 */
FitTwoStraightLines::FitTwoStraightLines(const size_t intervalCount, 
                                         const IntensityType robustnessFactor, 
                                         const std::string filename)
  : DetectKinkPoint(filename), 
    mIntervalCount(intervalCount), 
    mRobustnessFactor(robustnessFactor)
{
  mMethodName = "Two lines are fitted in a plausible interval";
}

/**
 * Executes the kink point detection.
 *
 * @param graph The (x,y)-pairs in the graph sorted ascending in x.
 * @return Kink point between specific and non-specific binding.
 */
IntensityPair FitTwoStraightLines::operator() (const IntensityMapping& graph)
{
  assert(graph.size() > 0);
  assert(mIntervalCount > 0);

  // Get interval borders and ascents
  GraphIterator leftPoint  = graph.begin();
  vector< pair<GraphIterator, IntensityType> > intervalsAndAscends = getIntervalsAndAscents(leftPoint, graph.end() - 1, mIntervalCount);

  // Get find leftmost robust maximum Index 
  // @todo Clean this up.
  size_t maxIndex = intervalsAndAscends.size() - 1;
  for (int i = intervalsAndAscends.size() - 2; i >= 0; --i) {
    if (atan(intervalsAndAscends[i].second) > mRobustnessFactor*atan(intervalsAndAscends[maxIndex].second)) {
      maxIndex = i;
    }
  } 

  // Set left index  
  size_t middleIntervalIndex = 0;
  if (maxIndex >= 2) { // Left point does not change iff the steepest ascent is within the first two intervals
    middleIntervalIndex = maxIndex-1; 
    leftPoint = intervalsAndAscends[middleIntervalIndex-1].first;
  } 
  
  assert(middleIntervalIndex + 1 < intervalsAndAscends.size()); 

  // Fit two straight lines
//   Polynomial<IntensityType> leftIntervalLine = MathUtil::fitStraightLine(leftPoint, 
//                                                                         intervalsAndAscends[middleIntervalIndex].first);
//   Polynomial<IntensityType> rightIntervalLine = MathUtil::fitStraightLine(intervalsAndAscends[middleIntervalIndex].first, 
//                                                                          intervalsAndAscends[middleIntervalIndex + 1].first);

  TwoPolynomials fittedLines = fitTwoLines(leftPoint, intervalsAndAscends[middleIntervalIndex + 1].first);

  // Compute intersection
  IntensityType intersectionPoint = (fittedLines.first[1] - fittedLines.second[1]) / 
      (fittedLines.second[0] - fittedLines.first[0]);
  IntensityType intersectionPointYvalue = fittedLines.first[0]*intersectionPoint + fittedLines.first[1]; // y = a*x + b
  
  // Export results
  if (mDoExport) {
    // Export results of hookcurve fitting      
    printVectorPairsToFile(graph, mExportFilename); 
    string filenameBase = stringutil::splitString(mExportFilename, ".")[0];

    string additionalCommands = "";

    DataCollector::instance().insert(filenameBase+"-intersectionPointXvalue" , intersectionPoint);
    DataCollector::instance().insert(filenameBase+"-intersectionPointYvalue" , intersectionPointYvalue);
    
    additionalCommands += ", " + fittedLines.first.toString() + " title 'Line Fit' lt 5 lw 2";
    additionalCommands += ", " + fittedLines.second.toString() + " title '' lt 5 lw 3";
    DataCollector::instance().insert(filenameBase+"-plottingCommands" , additionalCommands);
  }
  IntensityPair ip = IntensityPair(intersectionPoint, intersectionPointYvalue);
  ip = checkAndCorrectKinkPoint(ip, graph);
  return ip;
}


/**
 * 
 * 
 */
// FitTwoStraightLines::FitTwoStraightLinesByDerivation(const size_t slidingWindowSize, 
//                                          const std::string filename)
//   : DetectKinkPoint(filename), 
//     slidingWindowSize(slidingWindowSize)
// {
// }

// IntensityPair FitTwoStraightLinesByDerivation::operator() (const IntensityMapping& graph)
//   {
//     assert(graph.size() > 0);
//     assert(mIntervalCount > 0);
// 
// // Get interval borders and ascents
//     GraphIterator leftPoint  = graph.begin();
//     for (GraphIterator currentPoint = graph.begin(); currentPoint + slidingWindowSize < graph.end(); ++currentPoint) {
//       Polynomial<IntensityType> currentLine = MathUtil::fitStraightLine(currentPoint, currentPoint + slidingWindowSize);  
//       cout << currentLine[0] << endl;
//     }
//     
//         
//     vector< pair<GraphIterator, IntensityType> > intervalsAndAscends = getIntervalsAndAscents(leftPoint, graph.end() - 1, mIntervalCount);
// 
// // Export results
//     if (mDoExport) {
//   // Export results of hookcurve fitting      
//       printVectorPairsToFile(graph, mExportFilename); 
// 
//   // Debug for "plot "testplot_shuffle.dat" using 1:2, "breakpointDebug.dat" using 1:2"
//       ofstream gnuplotFile;  
//       gnuplotFile.open((mExportFilename + string(".gnuplot")).c_str());
// 
//       gnuplotFile << "set title 'intersectionPoint: " << intersectionPoint << "/" << intersectionPointYvalue << "'\n";
//       gnuplotFile << "set grid\n";
//       gnuplotFile << "show grid\n";
//       gnuplotFile << "set xrange[0:5]\nset yrange[-0.2:1]\n";
//       gnuplotFile << "set xlabel '0.5(log(PM)+Log(MM))'\nset ylabel 'log(PM)-log(MM)'\n";
//       gnuplotFile << "set arrow from " << intersectionPoint << ", graph 0 to " << intersectionPoint << " , graph 1 nohead lt -1\n";
//       gnuplotFile << "plot " << fittedLines.first.toString() << " title 'Straight Line Fit 1' " << ", " << fittedLines.second.toString() << " title 'Straight Line Fit 2'";
//       gnuplotFile << ", \"" << mExportFilename << "\" using 1:2 w lines title '" <<  mExportFilename << "'" << endl;
//       gnuplotFile << "pause -1 \"Hit return to continue\"" << endl;   
//       gnuplotFile.close();
//     }
//     IntensityPair ip = IntensityPair(intersectionPoint, intersectionPointYvalue);
//     ip = checkAndCorrectKinkPoint(ip, graph);
//     return ip;
//   }
    
    
    
 
/**
 * Devides a graph consisting of points leftPoint..rightPoint into intervalCount 
 * intervals and returns a list of pairs containing the last point and the 
 * ascent for each interval.
 * 
 * @param leftPoint First point of the graph
 * @param rightPoint Last point of the graph
 * @param intervalCount Number of intervals.
 */
std::vector< std::pair<GraphIterator, IntensityType> > FitTwoStraightLines::getIntervalsAndAscents(const GraphIterator& leftPoint, 
                                                                                               const GraphIterator& rightPoint,
                                                                                               const size_t intervalCount)
{
  vector< pair<GraphIterator, IntensityType> > intervalsAndAscends;
  // Add first interval border to avoid special cases
  intervalsAndAscends.push_back(make_pair(leftPoint, 0.0)); 

  IntensityType intervalSize = (rightPoint->first - leftPoint->first) / intervalCount;
  GraphIterator currentPoint = leftPoint;    
  for (size_t i = 1; i < intervalCount; ++i) {      
    // Goto first point in next interval
    while ((++currentPoint)->first < leftPoint->first + i*intervalSize) {} 

    // Calculate ascent between first and last point of interval
//     IntensityType currentAscent = (currentPoint->second - intervalsAndAscends.back().first->second) / 
//       (currentPoint->first - intervalsAndAscends.back().first->first);

    IntensityType currentAscent = MathUtil::fitStraightLine(intervalsAndAscends.back().first, currentPoint + 1)[0];

    // cout << currentAscent << " " << currentPoint->first << " " << atan(currentAscent) << " " << endl;

    // Add pair <last interval point (border), ascent> to list
    intervalsAndAscends.push_back(make_pair(currentPoint, currentAscent));  
  }
  // Add ascent for rightmost interval
  IntensityType lastAscent = (rightPoint->second - intervalsAndAscends.back().first->second) / 
    (rightPoint->first - intervalsAndAscends.back().first->first);
  intervalsAndAscends.push_back(make_pair(rightPoint, lastAscent));       
  //cout << lastAscent << " " << currentPoint->first << " " << atan(lastAscent) << " " << endl;

  // Delete first interval (which was a dummy element)
  intervalsAndAscends.erase(intervalsAndAscends.begin(), intervalsAndAscends.begin() + 1);

  return intervalsAndAscends;
}


/**
 * Calculates the best least squares fitting of two linear functions,
 * such that the first fits datapoints \f$ 0..x_g \f$ and the second fits points 
 * \f$ x_g +1..n \f$. The dataset are x/y pairs sorted by the x value. The routine uses 
 * SVD (MathUtil::solveEquationsSvd) for solving the equations. 
 *
 * @param beginPoint First element of a sorted list of points.
 * @param endPoint End element of a sorted list of points.
 *
 * @return Parameters for both optimal fitting functions. 
 */
TwoPolynomials FitTwoStraightLines::fitTwoLines(const GraphIterator& beginPoint, 
                                                const GraphIterator& endPoint)
{
  // Allocate matrices and vectors for linear (2x2,2x1) 
  // least squares matrix and initialize with zeros (calloc instead alloc)
  math::Matrix<IntensityType> deltaFirst(2,2, 0.0), deltaSecond(2,2, 0.0);
  IntensityArray thetaFirst(0.0, 2), thetaSecond(0.0, 2);
 
  // Initialize linear fit for x_g = 0 (second line from 0 to n)
  size_t pointCount = endPoint - beginPoint;
  for (GraphIterator currentPoint = beginPoint; currentPoint != endPoint; ++currentPoint) {
    IntensityType x = currentPoint->first;
    IntensityType y = currentPoint->second;
    deltaSecond[ 0][ 0] += (x*x);
    deltaSecond[ 0][ 1] += (x);
    deltaSecond[ 1][ 0] += (x);
    deltaSecond[ 1][ 1] += (1);

    thetaSecond[ 0] += (y*x);
    thetaSecond[ 1] += (y);
  }

  // Set minimal error value to INF
  IntensityType minimalTotalError = numeric_limits<IntensityType>::infinity(); // GSL_POSINF;  

  // Remember best linear and quadratic function
  Polynomial<IntensityType> bestFirstFit(1);
  Polynomial<IntensityType> bestSecondFit(2);

  //  For all x_g from 0 to n - 2
  for (GraphIterator currentPoint = beginPoint; currentPoint != endPoint - 1; ++currentPoint) {
     // x,y are the values of current x_g
    IntensityType x = currentPoint->first;
    IntensityType y = currentPoint->second;

    // Add values for first linear fit (see Preibisch, S.42)
    deltaFirst[ 0][ 0] += (x*x);
    deltaFirst[ 0][ 1] += (x);
    deltaFirst[ 1][ 0] += (x);
    deltaFirst[ 1][ 1] += (1);

    thetaFirst[ 0] += (y*x);
    thetaFirst[ 1] += (y);

    // Substract values from second linear fit (see Preibisch, S.42)
    deltaSecond[ 0][ 0] -= (x*x);
    deltaSecond[ 0][ 1] -= (x);
    deltaSecond[ 1][ 0] -= (x);
    deltaSecond[ 1][ 1] -= (1);

    thetaSecond[ 0] -= (y*x);
    thetaSecond[ 1] -= (y);

    // Solve linear equations
    IntensityArray firstFit = MathUtil::solveEquationsSvd(deltaFirst, thetaFirst, 0);
    IntensityArray secondFit = MathUtil::solveEquationsSvd(deltaSecond, thetaSecond, 0);
//     IntensityArray firstFit = MathUtil::solveEquationsLU(deltaFirst, thetaFirst);
//     IntensityArray secondFit = MathUtil::solveEquationsLU(deltaSecond, thetaSecond);

    // Calculate error for first linear fit
    IntensityType errorFirst = 0;
    for (GraphIterator errorPoint = beginPoint; errorPoint <= currentPoint; ++errorPoint) {
      IntensityType dataX = errorPoint->first;
      IntensityType dataY = firstFit[0]*dataX + firstFit[1];  // y = a*x + b
      errorFirst += (dataY - errorPoint->second)*(dataY - errorPoint->second);
    }

    // Calculate error for first fit
    IntensityType errorSecond = 0;
    for (GraphIterator errorPoint = currentPoint + 1; errorPoint != endPoint; ++errorPoint) {
      IntensityType dataX = errorPoint->first;
      IntensityType dataY = secondFit[0]*dataX + secondFit[1]; // y = a*x + b
      errorSecond += (dataY - errorPoint->second)*(dataY - errorPoint->second);
    }

    // Calculate weighted total error
    size_t currentIndex = currentPoint - beginPoint + 1;
    IntensityType totalError = errorFirst/(currentIndex) + errorSecond/(pointCount - currentIndex);
    
    // If current total error smaller than current minimal error, reset minimal and index
    if (totalError < minimalTotalError) {
      minimalTotalError = totalError;
      
      // Remember minimal values
      bestFirstFit = firstFit;
      bestSecondFit = secondFit; // = Polynomial<IntensityType>(secondFit);
    }
  }

  return make_pair(bestFirstFit, bestSecondFit);
} 

/**
 * Initializer.
 *
 * @param detectKinkPointFunctor function that detects the kink 
 *        point, e.g. FitTwoLines
 * @param digitizeIntervals Number of points the plot is digitized to
 *
 */
FitDigitizedPlot::FitDigitizedPlot(DetectKinkPointPtr detectKinkPointFunctor, 
                 const size_t digitizeIntervals) 
  : mDetectKinkPoint(detectKinkPointFunctor), mDigitizeIntervals(digitizeIntervals)
{
  mMethodName = detectKinkPointFunctor->getMethodName();
}

/**
 * Runs one of the other detection methods but digitizes the plot 
 * to a lower number of points before. 
 *
 * @param graph The (x,y)-pairs in the graph sorted ascending in x.
 * @return Kink point between specific and non-specific binding.
 *
 */
IntensityPair FitDigitizedPlot::operator() (const IntensityMapping& graph)
{
  return (*mDetectKinkPoint)(MathUtil::digitizeCurveRobust(graph,mDigitizeIntervals));
}

/**
 * Sets the name of the export file for the detectKinkPoint functor. If no name 
 * is set, nothing will be exported.
 *
 * @param filename File to export results of the detection.
 */
void FitDigitizedPlot::setExportFilename(const std::string filename)
{
  mDetectKinkPoint->setExportFilename(filename);
}
