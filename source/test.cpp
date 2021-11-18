/**
 * @file test.cpp
 * @author $Author$ 
 * @date $Date$
 *  
 * @note The testing suite currently is outdated.
 */
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
using boost::unit_test::test_suite;

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>

// #include "OpticalBackground.hpp"
#include "MathUtil.hpp"
#include "Matrix.hpp"
#include "ValarrayUtil.hpp"
#include "AveragingProcedure.hpp"
#include "PmMmProbe.hpp"
#include "SensitivityProfile.hpp"
#include "DetectKinkPoint.hpp"


#include <boost/function.hpp>

using namespace std;
using namespace larrpack;

/**
 * Tests functions for optical background calculation
 *
 * @todo Find out how to test private functions...
 */
void test_optical_background()
{
//   // Test getZone routine
//   BOOST_CHECK_EQUAL( 0, OpticalBackground::getZone(0,0,7,7,3)); 
//   BOOST_CHECK_EQUAL( 0, OpticalBackground::getZone(0,1,7,7,3));
//   BOOST_CHECK_EQUAL( 0, OpticalBackground::getZone(1,0,7,7,3));
//   BOOST_CHECK_EQUAL( 0, OpticalBackground::getZone(1,1,7,7,3));
//   BOOST_CHECK_EQUAL( 2, OpticalBackground::getZone(4,1,7,7,3));
//   BOOST_CHECK_EQUAL( 2, OpticalBackground::getZone(7,1,7,7,3));
//   BOOST_CHECK_EQUAL( 2, OpticalBackground::getZone(7,1,7,7,3));
//   BOOST_CHECK_EQUAL( 8, OpticalBackground::getZone(4,4,7,7,3));
//   BOOST_CHECK_EQUAL( 8, OpticalBackground::getZone(4,7,7,7,3));
//   BOOST_CHECK_EQUAL( 8, OpticalBackground::getZone(7,7,7,7,3));

//   // Test getZoneCenter routine
//   pair<double,double> center1(0,0);
//   BOOST_CHECK( (pair<double,double>(2.5,9.0)) == OpticalBackground::getZoneCenter(0,1,12,12,2) ); 
//   BOOST_CHECK( (pair<double,double>(24.5, 24.5)) == OpticalBackground::getZoneCenter(0,0,200,200,4) ); 

  
//   cout << OpticalBackground::getZoneCenter(0,0,200,99,4).first << ", " << OpticalBackground::getZoneCenter(1,0,11,100,2).second << endl;
//   cout << 13.5 << endl;



//   // reports 'error in "free_test_function": test 2 == 1 failed'
//   BOOST_CHECK(2 == 1); // non-critical test => continue after failure

}


/**
 * Computes z = x^2 + y^2 to test gradient descent
 */
double quadFunction(const valarray<double>& z) 
{
  return z[0]*z[0] + z[1]*z[1]; // return z^2
}



/**
 * A normal distribuition to test the integration
 */
double normalDistrib(double x, void *p){
  double mu    = 0.0;
  double std   = 2.0;
  return (1.0/sqrt(2 * M_PI * pow(std , 2))) * exp( - (pow((x - mu) , 2) / (2 * pow(std , 2)) ));
}
  
/** 
 * A simple function to test the integration
 */
double intfunc(double x, void *p){
  return 2*x+1;
}

/**
 * Tests mathematical functions
 *
 */
void test_math()
{
  // Test the svd solver
  // ---------------------------------------------------------------------------------------
  math::Matrix<double> delta(3,3);

  delta[0][0] = 1.0; delta[0][1] = 5.0; delta[0][2] = 2.0;
  delta[1][0] = 3.0; delta[1][1] = 7.0; delta[1][2] = 1.0;
  delta[2][0] = 4.0; delta[2][1] = 8.0; delta[2][2] = 15.0;

  valarray<double> theta(3);
  theta[0] = 13.0; theta[1] = 16.0; theta[2] = 61.0;

  valarray<double> expectedSolution(3);
  expectedSolution[0] = 2.0; expectedSolution[1] = 1.0; expectedSolution[2] = 3.0;

  // Test if expected and calculated sulution are equal
  // cout << MathUtil::solveEquationsSvd(delta, theta) - expectedSolution << endl;

  // @todo: Find out why numerical is soo bad now? (It once was fine with epsilon)
  valarray<bool> isSolutionEqual = (MathUtil::solveEquationsSvd(delta, theta) - expectedSolution < 10*numeric_limits<double>::epsilon());
  // cout << valarray<double>(MathUtil::solveEquationsSvd(delta, theta) - expectedSolution) << endl;
   
  BOOST_CHECK( isSolutionEqual.min() == true ); // Since false < true, this means that all el. are true



  // Test the averaging routines
  // --------------------------------------------------------------------------------------
  cout << "Computing Averages" << endl;
  double testvalues[] = {-10.0, 80.0, 30.0, 40.0, 50.0, 20.0};
  size_t maxx = 6;  
  valarray<double> testvalarray(&testvalues[0], (size_t)maxx);
  cout << testvalarray << endl;

  boost::function<double (valarray<double>)> computeAverage = Mean<double>();
  cout << computeAverage(testvalarray) << endl;
  BOOST_CHECK_EQUAL( 35, computeAverage(testvalarray) );

  computeAverage = Median<double>();
  cout << computeAverage(testvalarray) << endl;
  BOOST_CHECK_EQUAL( 35, computeAverage(testvalarray) );

  computeAverage = OnestepTukeyBiweight<double>(5.0, 0.0001);
  cout << computeAverage(testvalarray) << endl;
  BOOST_CHECK_EQUAL( 35, computeAverage(testvalarray) );

  sort(&testvalarray[0], &testvalarray[maxx]);
  cout << testvalarray << endl; 


  // Test gradient descent
  valarray<double> initP(100.0,2);
  valarray<double> expectedSolutionGradient(0.0, 2);
  initP[1] = -5.0;
  std::cout << "Testing Gradient Descend" << std::endl;
  valarray<double> gradientSolution = MathUtil::computeGradientDescent(initP, quadFunction);
  cout << "Gradient solution: " << gradientSolution << " expected: " << expectedSolutionGradient << endl;
  BOOST_CHECK( abs((gradientSolution - expectedSolutionGradient).max()) < 10e-5);
  
  
  // Test moving average
  double a[] = {0.0,1.0,2.0,4.0,8.0,16.0,32.0,64.0,128.0,256.0};
  valarray<double> v (a,10);
  //cout << MathUtil::calculateGeneralizedMovingAverage(v,3) << endl;
  double b[] = { 1, 1, 2.33333, 4.66667, 9.33333, 18.6667, 37.3333, 74.6667, 149.333, 149.333 };
  //double b[] = { 1, 1, 2, 3, 4, 5, 6, 7, 8, 8 };
  valarray<double> result (b,10);
  BOOST_CHECK_EQUAL( result.size(), MathUtil::calculateGeneralizedMovingAverage(v,3).size());
  BOOST_CHECK_CLOSE( result[0], MathUtil::calculateGeneralizedMovingAverage(v,3)[0], 0.001);
  BOOST_CHECK_CLOSE( result[1], MathUtil::calculateGeneralizedMovingAverage(v,3)[1], 0.001);
  BOOST_CHECK_CLOSE( result[2], MathUtil::calculateGeneralizedMovingAverage(v,3)[2], 0.001);
  BOOST_CHECK_CLOSE( result[3], MathUtil::calculateGeneralizedMovingAverage(v,3)[3], 0.001);
  
  // Test glog and gExp
  BOOST_CHECK_CLOSE(2.0, MathUtil::gLog10(100,2), 0.1); // Since we don't expect the exakt log10 values the relative error is allowed to be 0.1%
  BOOST_CHECK_CLOSE(-2.0, MathUtil::gLog10(-100,2), 0.1);  
  
  BOOST_CHECK_CLOSE(100.0, MathUtil::gExp10(2,2), 0.1);
  BOOST_CHECK_CLOSE(-100.0, MathUtil::gExp10(-2,2), 0.1);
  
  
  
  // Test std deviation
  vector<double> vec;
  vec.push_back(1.0);
  vec.push_back(3.0);
  vec.push_back(5.0);
  vec.push_back(7.0);
  BOOST_CHECK_CLOSE(4.0 , MathUtil::calculateMean(vec), 0.001);
  BOOST_CHECK_CLOSE(2.58199, MathUtil::calculateStandardDeviation(vec), 0.001); // Sample standard deviation
  //cout << MathUtil::calculateStandardDeviation(vec) << endl;
  
  // Test simpson integration
  //cout << MathUtil::simpsonSimple(-10, 10, 20, normalDistrib, &a) << endl;
  BOOST_CHECK_CLOSE(1.0, MathUtil::simpsonSimple(-10, 10, 20, normalDistrib, &a), 0.001);
  BOOST_CHECK_CLOSE(2.0, MathUtil::simpsonSimple(0.0, 1.0, 20, intfunc, &a), 0.001);
  
  
  // Test digitzeCurve:
  std::vector< std::pair<double, double> > pairVec;
  pairVec.push_back(pair<double, double>(0.0,0.0));
  pairVec.push_back(pair<double, double>(0.1,0.01));
  pairVec.push_back(pair<double, double>(0.4,0.16));
  pairVec.push_back(pair<double, double>(0.5,0.25));
  pairVec.push_back(pair<double, double>(1.0,1.0));
  pairVec.push_back(pair<double, double>(1.2,1.44));
  pairVec.push_back(pair<double, double>(1.4,1.95));
  pairVec.push_back(pair<double, double>(1.5,2.25));
  pairVec.push_back(pair<double, double>(1.7,2.89));
  pairVec.push_back(pair<double, double>(2.0,4.0));
  pairVec.push_back(pair<double, double>(2.1,4.41));
  pairVec.push_back(pair<double, double>(2.2,4.84));
  
  std::vector< std::pair<double, double> > resVec = MathUtil::digitizeCurve(pairVec, 5, 1.0);
  
  std::vector< std::pair<double, double> > desiredValues;
  desiredValues.push_back(pair<double, double>(0.22,0.056691));
  desiredValues.push_back(pair<double, double>(0.66,0.25));
  desiredValues.push_back(pair<double, double>(1.1,1.22));
  desiredValues.push_back(pair<double, double>(1.54,2.36324));
  desiredValues.push_back(pair<double, double>(1.98,4.41645));
  BOOST_CHECK_EQUAL(resVec.size(), desiredValues.size());
  for (size_t i = 0; i < resVec.size() ; ++i){
    BOOST_CHECK_CLOSE(resVec[i].first, desiredValues[i].first, 0.001);
    BOOST_CHECK_CLOSE(resVec[i].second, desiredValues[i].second, 0.001);
  }
  
  
  
  

  
}


void test_PmMmProbe()
{
  PmMmProbe probe("ChromosomeName", 'A', "CCCCCCCCCCCCCCCCCCCCCCCCC", 5, 400.0, 300.0, pair<int,int>(5,5), pair<int,int>(7,7), "id");
  BOOST_CHECK_CLOSE(0.124939, probe.getDeltaLogI(), 0.001);
  BOOST_CHECK_CLOSE(5.07918, probe.getSumLogI(), 0.001);
  BOOST_CHECK("CCCCCCCCCCCCGCCCCCCCCCCCC" == probe.getSequenceMm());
  BOOST_CHECK("CCCCCCCCCCCCCCCCCCCCCCCCC" == probe.getSequence());
  BOOST_CHECK("id" == probe.getProbesetId());
  double ar[] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7};
  std::valarray<double> a(ar,100);
  BOOST_CHECK_EQUAL(75, probe.getSensitivityPm(a));
  BOOST_CHECK_EQUAL(77, probe.getSensitivityMm(a)); 
}

void test_Fitting(){
    // Test intersection point calculation
  std::vector< std::pair<double, double> > graph;
  // Line vs parabola:
  //Points of Parabola follow: f(x) = -8*(x-0.4)**2  + 0.7;
  //Points of line:      g(x) = x
  // Intersection Point : 0.134/0.134
    //line:
  graph.push_back(pair<double, double>(0.0,0.0));
  graph.push_back(pair<double, double>(0.02,0.02));
  graph.push_back(pair<double, double>(0.05,0.05));
  graph.push_back(pair<double, double>(0.07,0.07));
  graph.push_back(pair<double, double>(0.1,0.1));
  graph.push_back(pair<double, double>(0.12,0.12));
  graph.push_back(pair<double, double>(0.13,0.13));
  
  //parabola:
  graph.push_back(pair<double, double>(0.18,0.313));
  graph.push_back(pair<double, double>(0.2,0.38));
  graph.push_back(pair<double, double>(0.25,0.52));
  graph.push_back(pair<double, double>(0.32,0.6488));
  graph.push_back(pair<double, double>(0.44,0.6872));
  graph.push_back(pair<double, double>(0.5,0.62));
  graph.push_back(pair<double, double>(0.6,0.38));
  
  
  
  //graph.push_back(pair<double, double>(5,0.5));
  FitStraightLineAndParabola  lineAndParabola("testKinkParabola.txt");
  IntensityPair intersection = lineAndParabola(graph);
  cout << "Intersection line/parabola: " << intersection.first << " / " << intersection.second << endl;
  BOOST_CHECK_CLOSE(0.134, intersection.first, 0.1);
  BOOST_CHECK_CLOSE(0.134, intersection.second, 0.1);
  graph.clear();
  
  // Line vs line.
  graph.push_back(pair<double, double>(0.0,0.0));
  graph.push_back(pair<double, double>(1,0.1));
  graph.push_back(pair<double, double>(2,0.2));
  graph.push_back(pair<double, double>(3,0.4));
  graph.push_back(pair<double, double>(3.5,0.5));
  graph.push_back(pair<double, double>(4,0.6));
  graph.push_back(pair<double, double>(5,0.8));
  FitTwoStraightLines twoLines(2, 0.7, "testKinkLines.txt");
  intersection = twoLines(graph);
  cout << "Intersection two lines: " << intersection.first << " / " << intersection.second << endl;
  BOOST_CHECK_CLOSE(2.0, intersection.first, 0.1);
  BOOST_CHECK_CLOSE(0.2,intersection.second, 0.1);
}
 
 void test_Profiles(){
   // SensitivityProfiles:
   cout << "\n\nTest Sensitvity Profile: " << endl;
   double ar[] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7};
   std::valarray<double> a(ar,100);
   SensitivityProfile profile(25, Mean<IntensityType>() , 1);
   profile.setProfile(a);
   cout << "Model Rank: " << profile.getProfile() << endl;
   BOOST_CHECK_EQUAL(75, profile.getSequenceIncrement("CCCCCCCCCCCCCCCCCCCCCCCCC"));
   BOOST_CHECK_EQUAL(1, profile.getModelRank());
   BOOST_CHECK_EQUAL(100, profile.getProfile().size());
   
   /// @todo Testroutine for computeProfile
   //BOOST_CHECK_EQUAL(100, profile.getSigma().size());
   
   //BOOST_CHECK_EQUAL((size_t) 384, profile.getProfileSize(2, 25, 16));
 }
 
 
test_suite* init_unit_test_suite( int, char* [] ) {
    test_suite* test= BOOST_TEST_SUITE( "Unit test example 1" );

    // this example will pass cause we know ahead of time number of expected failures
    test->add( BOOST_TEST_CASE( &test_optical_background ), 1 /* expected one error */ );

    test->add( BOOST_TEST_CASE( &test_math ), 0 );
    
    test->add( BOOST_TEST_CASE( &test_PmMmProbe ), 0 );
    
    test->add( BOOST_TEST_CASE( &test_Fitting), 0 );
    
    test->add( BOOST_TEST_CASE( &test_Profiles), 0 );

    return test;
}

