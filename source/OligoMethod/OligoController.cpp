/**
 * @file OligoController.cpp Defines controller methods to perform various steps of 
 *       the OLIGO algorithm
 * @author $Author: mario $ 
 * @author Mario Fasold
 * @date $Date: 2007-04-20 16:55:02 +0200 (Fri, 20 Apr 2007) $
 */
#include <time.h>
#include <math.h>


#include <algorithm>

#include "OligoController.hpp"
#include "MicroarrayExperiment.hpp"
#include "ProbesetComposition.hpp"
#include "ProbeFilter.hpp"
#include "StlUtil.hpp"
#include "ValarrayUtil.hpp"
#include "StringUtil.hpp"
#include "MathUtil.hpp"
#include "AveragingProcedure.hpp"
#include "SubsetSensitivityProfile.hpp"
#include "DataCollector.hpp"
#include "Chip.hpp"
#include "HookModel.hpp"


using namespace std;
using namespace larrpack;
using namespace boost;

const IntensityType OligoController::kSpecificThresholdIncrement = 0.5f; 

/// Helper function that returns max(v, 2.0)
IntensityType setMinimumTo2(IntensityType v)
{
  return max(v, 2.0);
}

/**
 * Return averaged sumLogI values
 *
 *
 */
vector<IntensityType> OligoController::getAverageSumLogIArray(IntensityMode mode1, IntensityMode mode2)
{ 
  ProbesetIntensitiesArray sumLogIs = getSumLogIArray(mode1, mode2);
  vector<IntensityType> avgSumLogIs(sumLogIs.size());
  transform(sumLogIs.begin(), sumLogIs.end(), avgSumLogIs.begin(), computeAverageIntensity);
  return avgSumLogIs;
}


/*
 * Initializes the arrays with background corrected intensities 
 * and expression measures according to probeset composition.
 *
 * @note This routine usually takes a lot of time.*/
 
void OligoController::initializeIntensityArrays()
{
  // Initialize corrected intensities
  // Background corrected:
  mCorrectedIntensities[kPmI] = Probeset::initializeIntensitiesArray(mChip->getProbesets());
  mCorrectedIntensities[kMmI] = Probeset::initializeIntensitiesArray(mChip->getProbesets());
  
  // These will store the temporarilly sensitivity corrected. In the beginning they are the same
  mCorrectedIntensities[kSensitivityCorrectedPm] = Probeset::initializeIntensitiesArray(mChip->getProbesets());
  mCorrectedIntensities[kSensitivityCorrectedMm] = Probeset::initializeIntensitiesArray(mChip->getProbesets());
  
  // These are now initialised only if they are needed.
//   mCorrectedIntensities[kNsSubstractedFastDiff] = Probeset::initializeIntensitiesArray(mProbesets);
//   mCorrectedIntensities[kNsSubstractedGcrmaDiff] = Probeset::initializeIntensitiesArray(mProbesets);

  // Loop over all probesets 
  for (size_t probesetIndex = 0; probesetIndex < mChip->getProbesets().size(); ++probesetIndex) {
    const Probeset& currentProbeset = mChip->getProbesets()[probesetIndex];
    for (size_t probeIndex = 0; probeIndex < currentProbeset.getSize(); ++probeIndex) {            
      PmMmProbe& currentProbe = *currentProbeset.getProbePtr(probeIndex);
      (*mCorrectedIntensities[kPmI])[probesetIndex][probeIndex] = log10(currentProbe.getPm());
      (*mCorrectedIntensities[kMmI])[probesetIndex][probeIndex] = log10(currentProbe.getMm());
      (*mCorrectedIntensities[kSensitivityCorrectedPm])[probesetIndex][probeIndex] = log10(currentProbe.getPm());
      (*mCorrectedIntensities[kSensitivityCorrectedMm])[probesetIndex][probeIndex] = log10(currentProbe.getMm());
//       (*mCorrectedIntensities[kNsSubstractedFastDiff])[probesetIndex][probeIndex] = 0.0; // Need to have some initial values.
//       (*mCorrectedIntensities[kNsSubstractedGcrmaDiff])[probesetIndex][probeIndex] = 0.0;
      
    }
  }
}

ProbesetIntensitiesArrayPtr OligoController::createPseudoMmArray(const ProbesetIntensitiesArray& pmArray, const vector< vector<size_t> >& gcToIndexMap) {
	vector<IntensityType> gcToPseudoMmIntensity; // Holds for each gc value the median of the pseudMm probes.
	for (size_t i = 0; i < gcToIndexMap.size(); ++i) { // Iterates over all bg probesets
		vector<IntensityType> intensitiesOfSameGc;
		for (size_t j = 0; j < gcToIndexMap[i].size(); ++j ) { // Iterates
			size_t index = gcToIndexMap[i][j];
			for (size_t k = 0; k < pmArray[index].size(); ++k) {
				intensitiesOfSameGc.push_back(pmArray[index][k]);
			}
		}
		gcToPseudoMmIntensity.push_back(MathUtil::calculateMedian(convertVectorToValarray(intensitiesOfSameGc)));
	}
	
	ProbesetIntensitiesArrayPtr pseudoMms = Probeset::initializeIntensitiesArray(mChip->getProbesets()); // pmArray; // copy the array structure.
	for (size_t probesetIndex = 0; probesetIndex < mChip->getProbesets().size(); ++probesetIndex) {
	  const Probeset& currentProbeset = mChip->getProbesets()[probesetIndex];
	  for (size_t probeIndex = 0; probeIndex < currentProbeset.getSize(); ++probeIndex) {
		  string probeSeq = mChip->getProbesets()[probesetIndex].getProbe(probeIndex).getSequence(); 
		  size_t gc = stringutil::countSubstr("G", probeSeq) + stringutil::countSubstr("C", probeSeq);
		  (*pseudoMms)[probesetIndex][probeIndex] = gcToPseudoMmIntensity[gc];
		  //cout << gc << "    "  << gcToPseudoMmIntensity[gc] << endl;
	  }
//	  cout << "test" << endl;
//	  printValarray(cout, (*pseudoMms)[probesetIndex]);
	}
	return pseudoMms;
}


void OligoController::initializeIntensityArraysPmOnly(const vector< vector<size_t> >& gcToIndexMap)
{
  // Initialize corrected intensities
  // Background corrected:
  mCorrectedIntensities[kPmI] = Probeset::initializeIntensitiesArray(mChip->getProbesets());
  //mCorrectedIntensities[kPseudoMmI] = Probeset::initializeIntensitiesArray(mChip->getProbesets());
  
  // These will store the temporarilly sensitivity corrected. In the beginning they are the same
  mCorrectedIntensities[kSensitivityCorrectedPm] = Probeset::initializeIntensitiesArray(mChip->getProbesets());
  
  
  // These are now initialised only if they are needed.
//   mCorrectedIntensities[kNsSubstractedFastDiff] = Probeset::initializeIntensitiesArray(mProbesets);
//   mCorrectedIntensities[kNsSubstractedGcrmaDiff] = Probeset::initializeIntensitiesArray(mProbesets);

  // Loop over all probesets 
  for (size_t probesetIndex = 0; probesetIndex < mChip->getProbesets().size(); ++probesetIndex) {
    const Probeset& currentProbeset = mChip->getProbesets()[probesetIndex];
    for (size_t probeIndex = 0; probeIndex < currentProbeset.getSize(); ++probeIndex) {            
      PmMmProbe& currentProbe = *currentProbeset.getProbePtr(probeIndex);
      (*mCorrectedIntensities[kPmI])[probesetIndex][probeIndex] = log10(currentProbe.getPm());
      //(*mCorrectedIntensities[kMmI])[probesetIndex][probeIndex] = log10(currentProbe.getMm());
      //(*mCorrectedIntensities[kSensitivityCorrectedPm])[probesetIndex][probeIndex] = log10(currentProbe.getPm());
      //(*mCorrectedIntensities[kSensitivityCorrectedMm])[probesetIndex][probeIndex] = log10(currentProbe.getMm());
//       (*mCorrectedIntensities[kNsSubstractedFastDiff])[probesetIndex][probeIndex] = 0.0; // Need to have some initial values.
//       (*mCorrectedIntensities[kNsSubstractedGcrmaDiff])[probesetIndex][probeIndex] = 0.0;
      
    }
  }
  
  mCorrectedIntensities[kPseudoMmI] = createPseudoMmArray(*mCorrectedIntensities[kPmI], gcToIndexMap);
  
  // This takes once more time.... try to just copy it.
  mCorrectedIntensities[kSensitivityCorrectedPseudoMm] = createPseudoMmArray(*mCorrectedIntensities[kPmI], gcToIndexMap);
}



/**
 * Returns a Refernce !! to the vector holding the corrected log intensity values.
 * @todo: Get rid of returning a reference by time !
 * 
 * 
 */
const std::vector<ProbesetIntensitiesArrayPtr>& OligoController::getAllProbesetIntensitiesArrays()
{
	return mCorrectedIntensities;
}


/**
 * Updates those arrays that store the current sensitivity corrected probe intensities.
 * 
 * @param calculationType The type of sensitivity correction (nonspecific, all)
 * 
 * 
 */
void OligoController::updateCorrectedProbes(const HookcurveCalculationType& calculationType, bool print)
{
  if (calculationType == kRawIntensitiesOnly) {
    return;
  }
  
//  else if (calculationType == kAllProfiles) {
//	  vector<IntensityPair> tmpVec;
//	    for (size_t probesetIndex = 0; probesetIndex < mChip->getProbesetCount(); ++probesetIndex) {
//	    	(*mCorrectedIntensities[kSensitivityCorrectedPm])[probesetIndex]      
//	        =  (*mCorrectedIntensities[kPmI])[probesetIndex]    
//	        - mHookModel->calculateIncrements(mChip->getProbesets()[probesetIndex], 
//	        		(*mCorrectedIntensities[intensityModeSelector[kPm].raw])[probesetIndex], 
//	                                         kPm, true);
//	        (*mCorrectedIntensities[kSensitivityCorrectedMm])[probesetIndex]  
//	        = (*mCorrectedIntensities[kMmI])[probesetIndex]  
//	        - mHookModel->calculateIncrements(mChip->getProbesets()[probesetIndex], 
//	        		(*mCorrectedIntensities[intensityModeSelector[kMm].raw])[probesetIndex],  
//	                                         kMm, true);
//	        for (size_t probeIndex = 0; probeIndex < (*mCorrectedIntensities[kPmI])[probesetIndex].size(); ++probeIndex) {
//	        	tmpVec.push_back(IntensityPair((*mCorrectedIntensities[kSensitivityCorrectedPm])[probesetIndex][probesetIndex], (mHookModel->calculateIncrements(mChip->getProbesets()[probesetIndex], 
//		        		(*mCorrectedIntensities[intensityModeSelector[kPm].raw])[probesetIndex],  
//		                                         kPm))[probeIndex] ));
//	        }
//	    }
//	    printVectorPairsToFile(tmpVec, "testplot.txt");
//	 
//  }
  	

  else{ // This udates the provisional sensitivity corrected pm and mm values.

    std::vector<IntensityType> sumLogIs = getAverageSumLogIArray(kSensitivityCorrectedPm, kSensitivityCorrectedMm);
	      
    for (size_t probesetIndex = 0; probesetIndex < mChip->getProbesetCount(); ++probesetIndex) {
    	
    	(*mCorrectedIntensities[kSensitivityCorrectedPm])[probesetIndex]           
        =  (*mCorrectedIntensities[kPmI])[probesetIndex]    
        - mHookModel->calculateIncrements(mChip->getProbesets()[probesetIndex], 
                                          (*mCorrectedIntensities[intensityModeSelector[kPm].raw])[probesetIndex], 
                                          kPm, sumLogIs[probesetIndex]);

    	(*mCorrectedIntensities[kSensitivityCorrectedMm])[probesetIndex]  
        = (*mCorrectedIntensities[kMmI])[probesetIndex]  
        - mHookModel->calculateIncrements(mChip->getProbesets()[probesetIndex], 
        		(*mCorrectedIntensities[intensityModeSelector[kMm].raw])[probesetIndex],  
                                          kMm, sumLogIs[probesetIndex]);
    	
    }
    
  }
}







void OligoController::updateCorrectedProbesPmOnly(const HookcurveCalculationType& calculationType, bool print)
{
  if (calculationType == kRawIntensitiesOnly) {
    return;
  }

  else{ // This udates the provisional sensitivity corrected pm and mm values.
    std::vector<IntensityType> sumLogIs = getAverageSumLogIArray(kSensitivityCorrectedPm, kSensitivityCorrectedPseudoMm);

    for (size_t probesetIndex = 0; probesetIndex < mChip->getProbesetCount(); ++probesetIndex) {
    	
    	(*mCorrectedIntensities[kSensitivityCorrectedPm])[probesetIndex]           
        =  (*mCorrectedIntensities[kPmI])[probesetIndex]    
        - mHookModel->calculateIncrements(mChip->getProbesets()[probesetIndex], 
        		(*mCorrectedIntensities[kPmI])[probesetIndex], 
                                          kPm, sumLogIs[probesetIndex]);

    }
    
  } // Else
}







/**
 * Reads microarray probe data depending on chip type 
 *
 * @return Probes
 */
PmMmProbePtrVector OligoController::importProbeData(IntensityMatrixPtr intensityArray,BoolArrayPtr intensityMaskArray)
{
  // Read each format
  MicroarrayExperiment microExp;
  if (options.chipType != kChiptypeGenomicfile) { 
    // Read celfile - if not intensityArray provided by parameter
    if (!intensityArray) {
      MicroarrayExperiment::readArrayFromCelfile(options.intensityFilename, intensityArray, intensityMaskArray);
    }

    if (intensityArray->rowNum() <= 1) {
      cerr << "Error reading CEL-file." << endl;
      exit(0);      
    }

    DataCollector::instance().insert("chipRowNumber", intensityArray->rowNum());
    DataCollector::instance().insert("chipColumnNumber", intensityArray->columnNum());

    // MicroarrayExperiment::readProbesetIdsFromCdffile(string(argv[1]));
    
    // Perform background correction
    options.backgroundCorrection->computeCorrection(*intensityArray, *intensityMaskArray);

    // @debug Set all intensities <2 to 2
    if (options.universalParam == "setMinimumIntensity") {
      valarray<IntensityType> tmpArray = intensityArray->flat().apply(setMinimumTo2);
      intensityArray->flat() = tmpArray;
    }

    // Read chip specific data
    switch (options.chipType) {
    case kChiptypeGenechip:
      // For genechips: read info from tabular files
      microExp.readFromProbesequenceFile(options.chipDiscriptionFilename, *intensityArray);
      break;
    case kChiptypeTiling: // For tiling arrays: read either PMAP or BPMAP format, depending on filename suffix
      if (stringutil::endsWith(stringutil::getUppercase(options.chipDiscriptionFilename), "BPMAP")) { 
        microExp.readFromBpmapfile(options.chipDiscriptionFilename, *intensityArray);
      } 
      else { // does not end with BPMAP
        microExp.readFromPmapfile(options.chipDiscriptionFilename, *intensityArray);
      }
      break;
    case kChiptypeSnp:
      microExp.readFromSnpSequenceFile(options.chipDiscriptionFilename, *intensityArray);
      break;
    case kChiptypeExon:
      //microExp.readFromExonSequenceFile(options.chipDiscriptionFilename, *intensityArray);
    	microExp.readPmOnlySequenceFile(options.chipDiscriptionFilename, *intensityArray);
      break;
    default:
      cerr << "No valid chip type" << endl;
      exit(0);      
    }
  }
  else {   // read .genomic file (developmental)
    cerr << "Warning: No background correction is performed" << endl;
    microExp.readFromGenomicFile(options.intensityFilename);
  }

  cout << "Reading success" << endl;

  // Extract PM/MM probes only
  PmMmProbePtrVector pmMmProbes = extractPmMmProbesFromProbelist(microExp.getProbes()); 

  // Debug: Remove parts of the probes 
  if (options.probesetSizeFilter > 1000) {
    cout << "Waning: Debug filter applied. Using only every" 
         << (options.probesetSizeFilter - 1000) << "th probe." << endl;
    pmMmProbes.erase(remove_if(pmMmProbes.begin(), pmMmProbes.end(), 
                               NotEveryXthElement<ProbePtr>(options.probesetSizeFilter - 1000)), 
                     pmMmProbes.end());
  }
  
  // Remove invalid probes (with intensity < 2)
  /// @todo: Do not delete them! Ignore them while calculating the SequenceProfiles!
  pmMmProbes.erase(remove_if(pmMmProbes.begin(), pmMmProbes.end(), ProbeintensityIsLessThan(2)),
                   pmMmProbes.end());

  // DEBUG: Remove probes beginning with GGG (index 22 for end...)
//   pmMmProbes.erase(remove_if(pmMmProbes.begin(), pmMmProbes.end(), ProbesequenceHasSubstring("GG",0)),
//                    pmMmProbes.end());
  
  if (options.maskfile != "") {
    // Remove probes that are in the list of probes given by the maskfile.
    ifstream maskfileStream;
    maskfileStream.open((options.maskfile).c_str());
    std::vector<std::string> maskedProbesets;
    string currentLine;
    // We discard the first two lines.
    //    getline(maskfileStream, currentLine);
    // getline(maskfileStream, currentLine);
    while (!maskfileStream.eof()) {
      getline(maskfileStream, currentLine);
      try {
        vector<string> split = stringutil::splitString(boost::lexical_cast<std::string>(currentLine));
        if (split.size() > 0) {
          std::string probesetId 
            = stringutil::splitString(boost::lexical_cast<std::string>(currentLine), " \t")[0];
          maskedProbesets.push_back(probesetId);
        }
      }
      catch(boost::bad_lexical_cast &) {
        // Ignore lines that can't be interpreted.
      }
    }
    // We reomove all probes matching the filter (before probesets are defined)
    pmMmProbes.erase(remove_if(pmMmProbes.begin(), pmMmProbes.end(), 
                               not1(IsInProbesetFilter(maskedProbesets))), 
                     pmMmProbes.end());
  }

//   pmMmProbes.erase(remove_if(pmMmProbes.begin(), pmMmProbes.end(), not1(ProbeChromosomeEquals(options.universalParam))), // "CHR11_st"))),
//                    pmMmProbes.end());

  return pmMmProbes;
}

/**
 * Set up probesets and chip object
 *
 *
 */
ChipPtr OligoController::createChip(PmMmProbePtrVector& probes)
{
  // Select probeset composition for each chip type
  ProbesetCompositionFunction composeProbeset; // Function used for probeset composition
  switch (options.chipType) {
  case kChiptypeGenechip:
    composeProbeset = ComposeProbesetsById();  // expression arrays: compose probesets by their probeset id
    break;
  case kChiptypeGenomicfile: 
  case kChiptypeTiling:     // Compose probesets with shuffling method
  {
	IntensityArrayTable* intensityArrayPtr = &mCorrectedIntensities;
    composeProbeset = ComposeProbesetsByCompleteShufflingMethod(options.maxProbeDistance,
                                                                options.optimalProbesPerSet, 
                                                                options.maxProbesPerSet,
                                                                options.specificityShufflingTreshold, 
                                                                options.hookcurveMovingAverageWindowSize,
                                                                options.shufflingMovingAverageWindowSize, 
  }
    break;
  case kChiptypeSnp:  
    composeProbeset = 
      ComposeProbesetsByGenotypeCall(ComposeProbesetsByGenotypeCall::readGenotypeCalls
                                     (options.genotypeCallFilename, 
                                      DataCollector::instance().getAs<string>("chipId")), kAllele);
    break;
  case kChiptypeExon:  
    composeProbeset = ComposeProbesetsById();
    break;
  }

  
  // Create the Chip object
  mChip = ChipPtr(new Chip(probes, composeProbeset(probes)));
  return mChip;
}


HookModelPtr OligoController::createHookModel()
{
	mHookModel = HookModelPtr(new HookModel());
	return mHookModel;
}


/*
 * Temoprarily method
 * 
 */
HookcurveAnalyzerPtr OligoController::createHookcurveModel()
{
	mHookcurveModel = HookcurveAnalyzerPtr(new HookcurveAnalyzer(*mChip, *mHookModel/*, mCorrectedIntensities*/));
	return mHookcurveModel;
}

/*
 * This function calculates the I_{max} for PM-only arrays
 * 
 * @return I_{max} (as log10 value)
 * 
 * @todo A bit messy ... but using the fit of a straight line onto the mad vs median plot (like here) is provisional anyway... 
 * 
 * 
 */
IntensityType OligoController::calculateImaxPmOnlyArray()
{
	vector<IntensityPair> tmpPlot; // std vs mean plot
	vector<IntensityPair> tmpPlot1; // mad vs median plot
	IntensitiesFunction mean    = Mean<IntensityType>();
	IntensitiesFunction median = Median<IntensityType>();

	  for (size_t probesetIndex = 0; probesetIndex < mCorrectedIntensities[kSensitivityCorrectedPm]->size(); ++probesetIndex) {
	        
		  if ((*mCorrectedIntensities[kSensitivityCorrectedPm])[probesetIndex].size() > 1) { // Only calculate if we have enough probes in the probeset std is inrobust otherwise.
			                                                          
			  tmpPlot.push_back( IntensityPair(mean((*mCorrectedIntensities[kSensitivityCorrectedPm])[probesetIndex] ),
					  MathUtil::calculateStandardDeviation((*mCorrectedIntensities[kSensitivityCorrectedPm])[probesetIndex] )));
			  tmpPlot1.push_back( IntensityPair(median((*mCorrectedIntensities[kSensitivityCorrectedPm])[probesetIndex] ),
			  					  MathUtil::calculateMAD((*mCorrectedIntensities[kSensitivityCorrectedPm])[probesetIndex] )));
		  }
	  } // Probeset Iterator

	  // for the std vs mean plot:
	  sort(tmpPlot.begin(), tmpPlot.end());
	  IntensityMappingPtr intensityPlot =  MathUtil::movingAv2(tmpPlot, (size_t)91);
	  printVectorPairsToFile(*intensityPlot, "stdvsmean-plot.dat");
	  IntensityMapping intensityPlotDigitized = MathUtil::digitizeCurveNew(*intensityPlot);
	  IntensityMapping::const_iterator pt = intensityPlotDigitized.begin();
	  for (IntensityMapping::const_iterator pair = intensityPlotDigitized.begin(); pair < intensityPlotDigitized.end(); ++pair) {
		  if (pair->second > pt->second) pt = pair;
	  }
	  Polynomial<IntensityType> line = MathUtil::fitStraightLine(pt, intensityPlotDigitized.end());
//	  cout << "MAX " << pt->first << "/" << pt->second << endl; 
	  printVectorPairsToFile(intensityPlotDigitized, "stdvsmeanDigitized-plot.dat");
	  ofstream gnuplotFile;
	  gnuplotFile.open("stdvsmeanDigitized-plot.gnuplot");
	  gnuplotFile << "set term png" << endl;
	  gnuplotFile << "set output " << "\"stdvsmeanDigitized-plot.png\"" << endl;
	  gnuplotFile << "plot " << "\"stdvsmeanDigitized-plot.dat\"" << ", " << line.toString() << endl;
	  gnuplotFile.close();

	  // for the mad vs median plot:
	  sort(tmpPlot1.begin(), tmpPlot1.end());
	  IntensityMappingPtr intensityPlot1 =  MathUtil::movingAv2(tmpPlot1, (size_t)91);
	  printVectorPairsToFile(*intensityPlot1, "madvsmedian-plot.dat");
	  IntensityMapping intensityPlot1Digitized = MathUtil::digitizeCurveNew(*intensityPlot1);
	  IntensityMapping::const_iterator pt1 = intensityPlot1Digitized.begin();
	  	  for (IntensityMapping::const_iterator pair = intensityPlot1Digitized.begin(); pair < intensityPlot1Digitized.end(); ++pair) {
	  		  if (pair->second > pt1->second) pt1 = pair;
	  	  }
	  Polynomial<IntensityType> line1 = MathUtil::fitStraightLine(pt1, intensityPlot1Digitized.end());
	  //	  cout << "MAX " << pt->first << "/" << pt->second << endl; 
	  printVectorPairsToFile(intensityPlot1Digitized, "madvsmedianDigitized-plot.dat");
	  ofstream gnuplotFile1;
	  gnuplotFile1.open("madvsmedianDigitized-plot.gnuplot");
	  gnuplotFile1 << "set term png" << endl;
	  gnuplotFile1 << "set output " << "\"madvsmedianDigitized-plot.png\"" << endl;
	  gnuplotFile1 << "plot " << "\"madvsmedianDigitized-plot.dat\"" << ", " << line1.toString() << endl;
	  gnuplotFile1.close();
	
	  // MAD vs median seems to be more robust.
	return - line1[1] / line1[0]; // Returns the intersection point of the calculated line with y(x) = 0 line. (The x axis) 
}


//IntensityType getGC(string seq)
//{
//	return ((double) (countSubstr("G", seq) + countSubstr("C", seq))) / ((double) seq.size()); 
//}

/*
 * Function that gives back a complete ProbesetIntensitiesArray, holding for each probe it's GC contant
 * @Debug  interesting for analytical purposes.
 * 
 */
ProbesetIntensitiesArrayPtr getGCContantArray(const vector<Probeset>& probesets)
{
//	int count = 0;
	ProbesetIntensitiesArrayPtr gcArray = Probeset::initializeIntensitiesArray(probesets);
	for (size_t probesetIndex = 0; probesetIndex < probesets.size(); ++probesetIndex) {
		valarray<IntensityType> v(0.0,probesets[probesetIndex].getSize());
		for (size_t probeIndex = 0; probeIndex < probesets[probesetIndex].getSize(); ++probeIndex) {
			string seq = probesets[probesetIndex].getProbe(probeIndex).getSequence();
			string g = "G";
			string c = "C"; 
			v[probeIndex] = ((double) (stringutil::countSubstr(g, seq) + stringutil::countSubstr(c, seq))) / ((double) seq.size());
		}
		(*gcArray)[probesetIndex] = v;
//		if (count <10){
//			printValarray(cout, (*gcArray)[probesetIndex]);
//			++count;
//		}
	}
	return gcArray;
}


/**
 * Computes the delta-sum plot.
 *
 * @param calculationType 
 */
IntensityMappingPtr OligoController::computeHookcurve(const HookcurveCalculationType& calculationType,
                                                      HookcurveAnalyzerPtr hookcurveAnalyzer)
{
  updateCorrectedProbes(calculationType); // Correct the pm and mm values by the current sensitivityprofiles @todo -> This should better be done in the main routine.
  calculateNsChip(); // Update the intersection point characteristics.

  // @Todo: Prevent from copying the array if we don't need to (In exon arrays we only need the pointer)
//  ProbesetIntensitiesArray sumLogs;
//  if (options.chipType == kChiptypeExon) {
//	  sumLogs = *mCorrectedIntensities[kSensitivityCorrectedPm];
//  }
//  else {
  ProbesetIntensitiesArray sumLogs   =  getSumLogIArray(kSensitivityCorrectedPm, 
	                                                          kSensitivityCorrectedMm);  
//  }
  
  ProbesetIntensitiesArray deltaLogs = getDeltaLogIArray(kSensitivityCorrectedPm, 
                                                         kSensitivityCorrectedMm);
  // @Debug: shows the gc contant as a function of sumLogI
  //ProbesetIntensitiesArrayPtr gcContant = getGCContantArray(mChip->getProbesets());
  //  IntensityMappingPtr tmpdata = hookcurveAnalyzer->getAveragedHookcurvePlot(sumLogs, *gcContant, options.hookcurveMovingAverageWindowSize);
  
  // Save sumLogIs and deltaLogIs to mapping.
  IntensityMappingPtr tmpdata = hookcurveAnalyzer->getAveragedHookcurvePlot(sumLogs, deltaLogs, options.hookcurveMovingAverageWindowSize);

  return tmpdata; // Smart Pointers must not be temporary objects
}









IntensityMappingPtr OligoController::computeHookcurvePmOnly(const HookcurveCalculationType& calculationType,
                                                      HookcurveAnalyzerPtr hookcurveAnalyzer, const vector< vector<size_t> >& gcToIndexMap)
{
	ProbesetIntensitiesArray sumLogs, deltaLogs;
	
  //updateCorrectedProbes(calculationType); // Correct the pm and mm values by the current sensitivityprofiles @todo -> This should better be done in the main routine.
	  if (calculationType == kRawIntensitiesOnly) {// This udates the provisional sensitivity corrected pm values
		  sumLogs   =  getSumLogIArray(kPmI,  kPseudoMmI);
		  deltaLogs  =  getDeltaLogIArray(kPmI,   kPseudoMmI);
	  }
	  else {
      //std::vector<IntensityType> sumLogIs  = getAverageSumLogIArray(kSensitivityCorrectedPm,  kSensitivityCorrectedPseudoMm);

	    // for (size_t probesetIndex = 0; probesetIndex < mChip->getProbesetCount(); ++probesetIndex) {
	    // 	(*mCorrectedIntensities[kSensitivityCorrectedPm])[probesetIndex]           
	    //     =  (*mCorrectedIntensities[kPmI])[probesetIndex]    
	    //     - mHookModel->calculateIncrements(mChip->getProbesets()[probesetIndex], 
	    //     		(*mCorrectedIntensities[kPmI])[probesetIndex], 
      //                                       kPm, sumLogIs[probesetIndex]);
	    // }
	  mCorrectedIntensities[kSensitivityCorrectedPseudoMm] = createPseudoMmArray(*mCorrectedIntensities[kSensitivityCorrectedPm], gcToIndexMap);
	  sumLogs   =  getSumLogIArray(kSensitivityCorrectedPm,  kSensitivityCorrectedPseudoMm);
	  deltaLogs = getDeltaLogIArray(kSensitivityCorrectedPm,   kSensitivityCorrectedPseudoMm);
	  }	    	
	    	

 // calculateNsChip(); // Update the intersection point characteristics.

  // @Debug: shows the gc contant as a function of sumLogI
  //ProbesetIntensitiesArrayPtr gcContant = getGCContantArray(mChip->getProbesets());
  //  IntensityMappingPtr tmpdata = hookcurveAnalyzer->getAveragedHookcurvePlot(sumLogs, *gcContant, options.hookcurveMovingAverageWindowSize);
  
  // Save sumLogIs and deltaLogIs to mapping.
  IntensityMappingPtr tmpdata = hookcurveAnalyzer->getAveragedHookcurvePlot(sumLogs, deltaLogs, options.hookcurveMovingAverageWindowSize);

  return tmpdata; // Smart Pointers must not be temporary objects
}

/*
 * Creates a graph of gc contant vs sumLogI
 * 
 * @debug
 * 
 * 
 */
IntensityMappingPtr OligoController::calculateGCContant(IntensityMode pm, IntensityMode mm, HookcurveAnalyzerPtr hookcurveAnalyzer)
{
  
  ProbesetIntensitiesArray sumLogs   =  getSumLogIArray(pm, 
                                                        mm);
//  ProbesetIntensitiesArray deltaLogs = getDeltaLogIArray(kSensitivityCorrectedPm, 
//                                                         kSensitivityCorrectedMm);
  ProbesetIntensitiesArrayPtr gcContant = getGCContantArray(mChip->getProbesets());
  // Save sumLogIs and gcContant to mapping.
  IntensityMappingPtr tmpdata = hookcurveAnalyzer->getAveragedHookcurvePlot(sumLogs, *gcContant, options.hookcurveMovingAverageWindowSize);
  
  return tmpdata; // Smart Pointers must not be temporary objects
}



/**
 * Computes coordinates of the kink-point
 *computeImaxAndA
 *
 * @todo Do not access hookcurveAnalyzer
 */
IntensityPair OligoController::computeKinkCoordinates(const std::string hookplotFilename, 
                                                      const IntensityMapping& hookcurvePlot,
                                                      HookModelPtr hookModel)
{
  // Calculate kink point
  options.detectKinkPoint->setExportFilename(hookplotFilename);
  IntensityPair nsThreshold = (*options.detectKinkPoint)(hookcurvePlot);
  // We export some additional information to yield a nice representation if we plot the kink point later on with gnuplot.
//   DataCollector::instance().insert(stringutil::splitString(hookplotFilename, ".")[0], options.detectKinkPoint->getMethodSpecificGnuplotCommands());

  // Set model nsThreshold here
  hookModel->setNsThreshold(nsThreshold);

  return nsThreshold;
}

/**
 * Returns the maximum sumLogI value ST such that there are at least 
 * minimumRequiredProbesets probesets with sumLogI >= ST
 *
 */
IntensityType getMaximumSpecificThreshold(const vector<IntensityType>& sumLogI, 
                                          const size_t minimumRequiredProbesets)
{
  vector<IntensityType> sortedSumLogI(sumLogI);
  sort(sortedSumLogI.begin(), sortedSumLogI.end());
  return sortedSumLogI[sortedSumLogI.size() - minimumRequiredProbesets];
}

std::pair<size_t, IntensityType> 
getOptimumRankAndSpecificThreshold(const vector<IntensityType>& sumLogI, 
                                   const IntensityType nsThreshold)
{
  vector<IntensityType> sortedSumLogI(sumLogI);
  sort(sortedSumLogI.begin(), sortedSumLogI.end());
  IntensityType maximumThresholdNn = sortedSumLogI[sortedSumLogI.size() 
                                                   - ProfileFactory::kMinimumSpecificProbesetsForNn];

  if (maximumThresholdNn > nsThreshold + 0.2) {
    return make_pair(2, maximumThresholdNn);
  } 
  else {
    return make_pair(1, sortedSumLogI[sortedSumLogI.size()
                                      - ProfileFactory::kMinimumSpecificProbesetsForN]);
  }
}


void OligoController::computeIntermediateSpecificProfiles(const IntensityType nsThreshold,
                                                          const std::string filenameSuffix, 
                                                          const size_t probesetLimit)
{
  vector<IntensityType> sumLogIs = getAverageSumLogIArray(kSensitivityCorrectedPm, 
                                                          kSensitivityCorrectedMm);

  // Get maximum nsThreshold (and rank) so that there are enough probesets to compute a profile
  size_t modelRank;
  IntensityType specificThreshold;
  tie(modelRank, specificThreshold) = getOptimumRankAndSpecificThreshold(sumLogIs, nsThreshold);

  // The predicate filters probesets on the basis of its *SumLogI* value
  // @note options.profileModelRank is omitted for specific profiles!
  SensitivityProfilePtr sensitivityProfileSPm = 
    ProfileFactory::createSpecificProfile(kProbeSequenceLength, modelRank, options.gstackTuples); 
  SensitivityProfilePtr sensitivityProfileSMm = 
    ProfileFactory::createSpecificProfile(kProbeSequenceLength, modelRank, options.gstackTuples);

  // This part calculates the specific sensitivity profiles as described in Preibisch.
  // Just clone instance with same rank
//   sensitivityProfileSPm = mHookModel->getProfile(kNsPm)->cloneWithNewRank(options.profileModelRank);
//   sensitivityProfileSMm = mHookModel->getProfile(kNsMm)->cloneWithNewRank(options.profileModelRank);
    
//   // Set center position of S Mm to zero.
//   sensitivityProfileSMm->zeroMiddlebase();


  // ------------------------------------------------------------------
  // Alternate method: compute specific profiles from actual specific probes!
  vector<Probeset> probesets = mChip->getProbesets();
  
  // Consider probesets with sumLogI > nsThreshold + 0.5 are considered "specific"
  // but make sure that there are enough specific probesets to compute a profile
  specificThreshold = 
    min(nsThreshold + kSpecificThresholdIncrement, specificThreshold);
  if (specificThreshold != nsThreshold + kSpecificThresholdIncrement) {
    cout << "Specific threshold has been decreased to " << specificThreshold << endl;
    cout << " (used modelRank = " << modelRank << ")" << endl;
  }
  IntensityPredicate hasSpecificSumLogI = bind2nd(greater<IntensityType>(), specificThreshold); 
	vector<Probeset> filteredProbesets = filterParallel(probesets, sumLogIs, hasSpecificSumLogI);

  cout << "\nCalculating S Pm Sensitivity Profile with " << flush;    
  sensitivityProfileSPm->computeProfile(filteredProbesets,  // WARNING: Lots of data (probesets) copied!
                                         filterParallel(*mCorrectedIntensities[kPmI], 
                                                        sumLogIs, hasSpecificSumLogI),
                                         mem_fun_ref(&PmMmProbe::getSequence),
                                         probesetLimit); 
    
  cout << "\nCalculating S Mm Sensitivity Profile " << flush;
  sensitivityProfileSMm->computeProfile(filteredProbesets, 
                                        filterParallel(*mCorrectedIntensities[kMmI], sumLogIs, hasSpecificSumLogI),
                                         mem_fun_ref(&PmMmProbe::getSequenceMm),
                                         probesetLimit);
  sensitivityProfileSPm->exportToDatafile(string("SensitivityProfileSPm-" + filenameSuffix + ".dat"));
  sensitivityProfileSMm->exportToDatafile(string("SensitivityProfileSMm-" + filenameSuffix + ".dat"));
  // -----------------------------------------------------------------

  mHookModel->setProfile(kSPm, sensitivityProfileSPm);
  mHookModel->setProfile(kSMm, sensitivityProfileSMm);  
}





/**
 * 
 * Small changes to the PM/MM Version.
 * @todo: Try to megre with PM/MM Version.
 * 
 */
void OligoController::computeIntermediateSpecificProfilesPmOnly(const IntensityType nsThreshold,
                                                          const std::string filenameSuffix, 
                                                          const size_t probesetLimit)
{
  vector<IntensityType> sumLogIs = getAverageSumLogIArray(kSensitivityCorrectedPm, 
                                                          kSensitivityCorrectedPseudoMm);

  // Get maximum nsThreshold (and rank) so that there are enough probesets to compute a profile
  size_t modelRank;
  IntensityType specificThreshold;
  tie(modelRank, specificThreshold) = getOptimumRankAndSpecificThreshold(sumLogIs, nsThreshold);

  // The predicate filters probesets on the basis of its *SumLogI* value
  // @note options.profileModelRank is omitted for specific profiles!
  SensitivityProfilePtr sensitivityProfileSPm = 
    ProfileFactory::createSpecificProfile(kProbeSequenceLength, modelRank, options.gstackTuples); 
//  SensitivityProfilePtr sensitivityProfileSMm = 
//    ProfileFactory::createSpecificProfile(kProbeSequenceLength, modelRank, options.gstackTuples);

  // This part calculates the specific sensitivity profiles as described in Preibisch.
  // Just clone instance with same rank
//   sensitivityProfileSPm = mHookModel->getProfile(kNsPm)->cloneWithNewRank(options.profileModelRank);
//   sensitivityProfileSMm = mHookModel->getProfile(kNsMm)->cloneWithNewRank(options.profileModelRank);
    
//   // Set center position of S Mm to zero.
//   sensitivityProfileSMm->zeroMiddlebase();


  // ------------------------------------------------------------------
  // Alternate method: compute specific profiles from actual specific probes!
  vector<Probeset> probesets = mChip->getProbesets();
  
  // Consider probesets with sumLogI > nsThreshold + 0.5 are considered "specific"
  // but make sure that there are enough specific probesets to compute a profile
  specificThreshold = 
    min(nsThreshold + kSpecificThresholdIncrement, specificThreshold);
  if (specificThreshold != nsThreshold + kSpecificThresholdIncrement) {
    cout << "Specific threshold has been decreased to " << specificThreshold << endl;
    cout << " (used modelRank = " << modelRank << ")" << endl;
  }
  IntensityPredicate hasSpecificSumLogI = bind2nd(greater<IntensityType>(), specificThreshold); 
	vector<Probeset> filteredProbesets = filterParallel(probesets, sumLogIs, hasSpecificSumLogI);

  cout << "\nCalculating S Pm Sensitivity Profile with " << flush;    
  sensitivityProfileSPm->computeProfile(filteredProbesets,  // WARNING: Lots of data (probesets) copied!
                                         filterParallel(*mCorrectedIntensities[kPmI], 
                                                        sumLogIs, hasSpecificSumLogI),
                                         mem_fun_ref(&PmMmProbe::getSequence),
                                         probesetLimit); 
    
//  cout << "\nCalculating S Mm Sensitivity Profile " << flush;
//  sensitivityProfileSMm->computeProfile(filteredProbesets, 
//                                        filterParallel(*mCorrectedIntensities[kMmI], sumLogIs, hasSpecificSumLogI),
//                                         mem_fun_ref(&PmMmProbe::getSequenceMm),
//                                         probesetLimit);
  sensitivityProfileSPm->exportToDatafile(string("SensitivityProfileSPm-" + filenameSuffix + ".dat"));
//  sensitivityProfileSMm->exportToDatafile(string("SensitivityProfileSMm-" + filenameSuffix + ".dat"));
  // -----------------------------------------------------------------

  mHookModel->setProfile(kSPm, sensitivityProfileSPm);
//  mHookModel->setProfile(kSMm, sensitivityProfileSMm);  
}











/**
 * Calculates those specific profiles that are needed for expression measure calculation.
 * @note This has been moved here from calculateSpecificPortions
 * 
 * @param sequenceLength Lenght of the probe sequences
 * @param modelRank rank of the model (1 for single, 2 for nearest neighbour)
 * @param profileTypes vector of ProfileTypes that describes which measures are to be calculated
 * @param probesetLimit The number of probesets that shall be taken into calculation for sequence profiles. 
 */
void OligoController::calculateSpecificProfiles(IntensityMode pmMode, IntensityMode mmMode, const size_t& sequenceLength,
                                                const std::vector<ExpressionMeasureProfileType> profileTypes,
                                                const size_t& probesetLimit,
                                                const IntensityType& nsThreshold) 
{
  vector<IntensityType> sumLogIs = getAverageSumLogIArray(pmMode, mmMode);

  // Get maximum nsThreshold (and rank) such that there are enough probesets to compute a profile
  size_t modelRank;
  IntensityType specificThreshold;
  tie(modelRank, specificThreshold) = getOptimumRankAndSpecificThreshold(sumLogIs, nsThreshold);

  // @Debug
  modelRank = 1;

  // Consider probesets with sumLogI > nsThreshold + 0.5 are considered "specific"
  // but make sure that there are enough specific probesets to compute a profile
  specificThreshold = 
    min(nsThreshold + kSpecificThresholdIncrement, specificThreshold);
  if (specificThreshold != nsThreshold + kSpecificThresholdIncrement) {
    cout << "Specific threshold has been decreased to " << specificThreshold;
    cout << " (used modelRank = " << modelRank << ")" << endl;
  }
  
  // @Debug
  //specificThreshold = nsThreshold+0.4;
  
  IntensityPredicate hasSpecificSumLogI = bind2nd(greater<IntensityType>(), specificThreshold); 
	vector<Probeset> filteredProbesets = filterParallel(mChip->getProbesets(), sumLogIs, hasSpecificSumLogI);
	
  // For each of the needed profile types
  for (vector<ExpressionMeasureProfileType>::const_iterator currentProfile = profileTypes.begin(); 
       currentProfile < profileTypes.end(); ++currentProfile) {

    // Create profile
    // @note options.profileModelRank is omitted for specific profiles!
    SensitivityProfilePtr profile 
      = ProfileFactory::createSpecificProfile(kProbeSequenceLength, modelRank, options.gstackTuples);

    // Compute profile
    ExpressionMeasureProfileParams params = expressionMeasureProfileParams[*currentProfile];
    profile->computeProfile(filteredProbesets, filterParallel(*mCorrectedIntensities[params.intensityMode], 
                                                              sumLogIs, hasSpecificSumLogI),
                            params.getSequenceFunction, probesetLimit); 

    // Set model profile
    mHookModel->setProfile(params.sProfile, profile);

    // Export full and single nucleotide profiles
	profile->exportToDatafile("SensitivityProfile" + params.filenameSuffix + ".dat");
    profile->cloneWithNewRank(1)->exportToDatafile("SensitivityProfile" + params.filenameSuffix 
                                                   + "-single.dat");
  }
  
}


/**
 * Calculates or updates the characteristics of the NS area.
 * @note: The mean of the ns distribution is caluclated as an average of all ns probes... 
 * while nsChip is the average of all tukey weighted ns probeset.
 * 
 */
void OligoController::calculateNsChip()
{
	   // The hookcurve parameters
	 IntensityType nsChipTemp = 0.0;
	 IntensityType nsPmChipTemp = 0.0;
	 IntensityType nsMmChipTemp = 0.0;
	 IntensityType logBTemp = 0.0;
	 size_t probesetCount = 0;

   ProbesetIntensitiesArray sumLogIs   =  getSumLogIArray(kSensitivityCorrectedPm, kSensitivityCorrectedMm);
   ProbesetIntensitiesArray deltaLogIs = getDeltaLogIArray(kSensitivityCorrectedPm, kSensitivityCorrectedMm);


	 for (size_t probesetIndex = 0; probesetIndex < mCorrectedIntensities[kSensitivityCorrectedPm]->size(); ++probesetIndex) {
		 IntensityArray correctedLogPms = (*(mCorrectedIntensities[kSensitivityCorrectedPm]))[probesetIndex];
		 IntensityArray correctedLogMms = (*(mCorrectedIntensities[kSensitivityCorrectedMm]))[probesetIndex];

		 if (computeAverageIntensity(sumLogIs[probesetIndex]) < mHookModel->getNsThreshold().first) { 
			 nsChipTemp   += computeAverageIntensity(sumLogIs[probesetIndex]);
			 nsPmChipTemp += computeAverageIntensity(correctedLogPms);
			 nsMmChipTemp += computeAverageIntensity(correctedLogMms);
			 logBTemp     += computeAverageIntensity(deltaLogIs[probesetIndex]);
			 probesetCount++;
		 }
	 }// For

  mHookModel->setParam("NsChip", nsChipTemp   / probesetCount);
  mHookModel->setParam("NsPmChip", nsPmChipTemp / probesetCount);
  mHookModel->setParam("NsMmChip", nsMmChipTemp / probesetCount);
  mHookModel->setParam("LogB", logBTemp     /probesetCount);
}




/*
 * Calculates the average intensity of the background probes
 * 
 * Calculates only over the PMs left of the threshold. (In the PM/MM version
 * it is the avergae over PM and MM (or Sum(PM/MM))
 * 
 * 
 */
void OligoController::calculateNsChipPmOnly()
{
	   // The hookcurve parameters
	 IntensityType nsChipTemp = 0.0;
	 IntensityType logBTemp = 0.0;
	 size_t probesetCount = 0;
	 
	 
	 // n

   ProbesetIntensitiesArray sumLogIs   =  getSumLogIArray(kSensitivityCorrectedPm, kSensitivityCorrectedPseudoMm);
   ProbesetIntensitiesArray deltaLogIs = getDeltaLogIArray(kSensitivityCorrectedPm, kSensitivityCorrectedPseudoMm);


	 for (size_t probesetIndex = 0; probesetIndex < mCorrectedIntensities[kSensitivityCorrectedPm]->size(); ++probesetIndex) {
		 IntensityArray correctedLogPms = (*(mCorrectedIntensities[kSensitivityCorrectedPm]))[probesetIndex];

		 if (computeAverageIntensity(sumLogIs[probesetIndex]) < mHookModel->getNsThreshold().first) { 
			 nsChipTemp   += computeAverageIntensity(correctedLogPms);
			 logBTemp     += computeAverageIntensity(deltaLogIs[probesetIndex]);
			 probesetCount++;
		 }
	 }// For

  mHookModel->setParam("NsChip", nsChipTemp   / probesetCount);
  mHookModel->setParam("LogB", logBTemp     /probesetCount);
}





IntensityType getCheckedThreshold(const vector<IntensityType>& sumLogI, const IntensityType nsThreshold,
                                  const size_t modelRank)
{
  vector<IntensityType> sortedSumLogI(sumLogI);
  sort(sortedSumLogI.begin(), sortedSumLogI.end());

  IntensityType minimumThreshold = sortedSumLogI[ProfileFactory::getMinimumProbesetCountToEstimateModel(modelRank)];
  return max(minimumThreshold, nsThreshold);
}

/**
 * Computes or updates the nonspecific sequence profiles.
 * 
 * 
 */
void OligoController::computeSequenceProfiles(const IntensityType uncheckedNsThreshold,
                                              const std::string filenameSuffix, const size_t probesetLimit)
{
  vector<Probeset> probesets = mChip->getProbesets();
  vector<IntensityType> sumLogIs = getAverageSumLogIArray(kSensitivityCorrectedPm, kSensitivityCorrectedMm);

  // Check if there are enough probesets at SumLogI < nsThreshold to estimate a profile
  // @todo: If we actually obtain a new threshold we should list this in a log file.
  IntensityType nsThreshold = getCheckedThreshold(sumLogIs, uncheckedNsThreshold, options.profileModelRank);

  // The predicate filters probesets on the basis of its *SumLogI* value
  IntensityPredicate isLessThanIntersectionSumLogI = bind2nd(less<IntensityType>(), nsThreshold);
    //bind2nd(greater<IntensityType>(), nsThreshold + 0.4);
  
  // Define which sequence profile to use
  SensitivityProfilePtr sensitivityProfileNSPm 
    = ProfileFactory::createNsProfile(kProbeSequenceLength, options.profileModelRank, options.gstackTuples);

  SensitivityProfilePtr sensitivityProfileNSMm 
    = ProfileFactory::createNsProfile(kProbeSequenceLength, options.profileModelRank, options.gstackTuples);
  
  // Read / Compute the profiles 
  if (options.readProfilesFromFile) { // read 
    cout << "Reading Sensitivity Profiles from files in current directory... " << flush;      
    sensitivityProfileNSPm->importFromDatafile(string("SensitivityProfileNsPm-" + filenameSuffix + ".dat"));
    sensitivityProfileNSMm->importFromDatafile(string("SensitivityProfileNsMm-" + filenameSuffix + ".dat"));
    cout << "done" << endl;
  }
  else { // compute
    vector<Probeset> filteredProbesets = filterParallel(probesets, sumLogIs, isLessThanIntersectionSumLogI);

    
    
    cout << "\nCalculating NsPm Sensitivity Profile with " << flush;    
    sensitivityProfileNSPm->computeProfile(filteredProbesets,  // WARNING: Lots of data (probesets) copied!
                                           filterParallel(*mCorrectedIntensities[kPmI], 
                                                          sumLogIs, isLessThanIntersectionSumLogI), //mem_fun_ref(&Probeset::getPmLogIs), 
                                           mem_fun_ref(&PmMmProbe::getSequence),
                                           probesetLimit); 
    
    cout << "\nCalculating NsMm Sensitivity Profile " << flush;
  
    sensitivityProfileNSMm->computeProfile(filteredProbesets, 
                                           filterParallel(*mCorrectedIntensities[kMmI], sumLogIs, isLessThanIntersectionSumLogI),
                                           mem_fun_ref(&PmMmProbe::getSequenceMm),
                                           probesetLimit);
  } // Else
  cout << endl;
    
  // Export solution to file and print gnuplot command
  if (!options.readProfilesFromFile){
	cout << "Exporting single profiles" << endl;
    sensitivityProfileNSPm->exportToDatafile(string("SensitivityProfileNsPm-" + filenameSuffix + ".dat"));
    //sensitivityProfileNSMm->exportToDatafile(string("SensitivityProfileNsMm-" + filenameSuffix + ".dat"));
    // @bug ? Do not use temporary shared pointers
    sensitivityProfileNSPm->cloneWithNewRank(1)->exportToDatafile(string("SensitivityProfileNsPm-" + filenameSuffix + "-single.dat"));
    //sensitivityProfileNSMm->cloneWithNewRank(1)->exportToDatafile(string("SensitivityProfileNsMm-" + filenameSuffix + "-single.dat"));
   }

  // @debug Print sum of squares error
  if (options.isDebugMode) {
    vector<Probeset> probesets = mChip->getProbesets();
    vector<IntensityType> sumLogIs = getAverageSumLogIArray(kSensitivityCorrectedPm, kSensitivityCorrectedMm);
    vector<Probeset> filteredProbesets = filterParallel(probesets, sumLogIs, isLessThanIntersectionSumLogI);

    // Assume these profiles are simple profiles
//    SimpleSensitivityProfilePtr simpleProfileNSPm 
//      = boost::shared_static_cast<SimpleSensitivityProfile>(sensitivityProfileNSPm);
//    SimpleSensitivityProfilePtr simpleProfileNSMm 
//      = boost::shared_static_cast<SimpleSensitivityProfile>(sensitivityProfileNSMm);
    SensitivityProfilePtr simpleProfileNSPm = sensitivityProfileNSPm;
    if (!simpleProfileNSPm) {
      exit(0);
    }

//    cout << "Sum of Squares Error (NsPm - Pm): " 
//         << simpleProfileNSPm->computeSumOfSquares(filteredProbesets, 
//                                                   filterParallel(*mCorrectedIntensities[kPmI], 
//                                                                  sumLogIs, isLessThanIntersectionSumLogI), 
//                                                   mem_fun_ref(&PmMmProbe::getSequence))
//         << endl;
//    cout << "Sum of Squares Error (NsMm - Mm): " 
//         << simpleProfileNSMm->computeSumOfSquares(filteredProbesets, 
//                                                   filterParallel(*mCorrectedIntensities[kMmI], 
//                                                                  sumLogIs, isLessThanIntersectionSumLogI), 
//                                                   mem_fun_ref(&PmMmProbe::getSequenceMm))
//         << endl;
  
    cout << "Export F-Test Sum of Squares Error (NsPm - Pm): " << endl; 
    SimpleSensitivityProfile::exportConditionalSumOfSquares(filteredProbesets, 
                                                            filterParallel(*mCorrectedIntensities[kPmI], 
                                                                           sumLogIs, isLessThanIntersectionSumLogI), 
                                                            mem_fun_ref(&PmMmProbe::getSequence),
                                                            1 // Compute SSE for all triples
                                                            , "SseProfile1.dat", sensitivityProfileNSPm);
    SimpleSensitivityProfile::exportConditionalSumOfSquares(filteredProbesets, 
                                                            filterParallel(*mCorrectedIntensities[kPmI], 
                                                                           sumLogIs, isLessThanIntersectionSumLogI), 
                                                            mem_fun_ref(&PmMmProbe::getSequence),
                                                            2 // Compute SSE for all triples
                                                            , "SseProfile2.dat", sensitivityProfileNSPm);
    SimpleSensitivityProfile::exportConditionalSumOfSquares(filteredProbesets, 
                                                     filterParallel(*mCorrectedIntensities[kPmI], 
                                                                    sumLogIs, isLessThanIntersectionSumLogI), 
                                                     mem_fun_ref(&PmMmProbe::getSequence),
                                                     3 // Compute SSE for all triples
                                                     , "SseProfile.dat", sensitivityProfileNSPm);
//     SimpleSensitivityProfile::exportConditionalSumOfSquares(filteredProbesets, 
//                                                             filterParallel(*mCorrectedIntensities[kMmI], 
//                                                                            sumLogIs, isLessThanIntersectionSumLogI), 
//                                                             mem_fun_ref(&PmMmProbe::getSequenceMm),
//                                                             3 // Compute SSE for all triples
//                                                             , "SseProfileMm.dat", sensitivityProfileNSMm);
    SimpleSensitivityProfile::exportConditionalSumOfSquares(filteredProbesets, 
                                                            filterParallel(*mCorrectedIntensities[kPmI], 
                                                                           sumLogIs, isLessThanIntersectionSumLogI), 
                                                            mem_fun_ref(&PmMmProbe::getSequence),
                                                            4 // Compute SSE for all triples
                                                            // options.profileModelRank
                                                            , "SseProfile4.dat", sensitivityProfileNSPm);

//     SimpleSensitivityProfile::exportConditionalSumOfSquares(filteredProbesets, 
//                                                             filterParallel(*mCorrectedIntensities[kMmI], 
//                                                                            sumLogIs, isLessThanIntersectionSumLogI), 
//                                                             mem_fun_ref(&PmMmProbe::getSequenceMm),
//                                                             4 // Compute SSE for all triples
//                                                             // options.profileModelRank
//                                                             , "SseProfile4.dat", sensitivityProfileNSMm);
  }

  // Apply profile
  mHookModel->setProfile(kNsPm, sensitivityProfileNSPm);
  mHookModel->setProfile(kNsMm, sensitivityProfileNSMm);
}



//vector<Probeset> getBackgroundProbesProbeset(probesets) 
//{
//	PmMmProbePtrVector probeset; // All background probes are pushed into one probeset
//	vector<Probeset> backProbeset;
//	for(vector<Probeset>::const_iterator psi = probesets.begin(); psi != probesets.end(); ++psi) {
//		if (psi->getProbesetId()[0, 15] == "BackgroundProbe"){
//			Probeset ps = *psi;
//			ps.setId
//			backProbeset.push_back(*psi);
//		}
//	}
//	return backProbeset;
//}




/**
 * Computes or updates the nonspecific sequence profiles for PM only chips
 * 
 * @todo Can this somehow be merged with computeSequenceProfiles ??
 * 
 * @param uncheckedNsThreshold The nsThreshold
 * @param filenameSuffix To export the profiles into a profile data and a gnuplot printing commands file
 * @param probesetLimit For using a limited number of probesets to calculate the profiles.
 * 
 */
void OligoController::computeSequenceProfilesPmOnlyChip(IntensityMode pmIntensityType, IntensityMode mmIntensityType, const IntensityType uncheckedNsThreshold,
                                              const std::string filenameSuffix/*, bool renameBackgroundProbeIds*/ /*(debug)e*/, const size_t probesetLimit)
{
	
	// @debug: Choose which probes are considered background for Sensitivity Profile estimation.
	bool controlAntiGenomic = false;
	bool controlIntron           = false;
	bool hookN                     = true;
	
	
  vector<Probeset> probesets = mChip->getProbesets();
//  vector<IntensityType> sumLogIs = getAverageSumLogIArray(kSensitivityCorrectedPm, kSensitivityCorrectedMm);
  // @todo Maybe we should use a averaged pm array instead of a pseudSumLogI array
  vector<IntensityType> sumLogIs = getAverageSumLogIArray(pmIntensityType, mmIntensityType); 

  // Check if there are enough probesets at SumLogI < nsThreshold to estimate a profile
  // @todo: If we actually obtain a new threshold we should list this in a log file.
  IntensityType nsThreshold = getCheckedThreshold(sumLogIs, uncheckedNsThreshold, options.profileModelRank);

  // The predicate filters probesets on the basis of its *SumLogI* value
  IntensityPredicate isLessThanIntersectionSumLogI = bind2nd(less<IntensityType>(), nsThreshold);
  // bind2nd(greater<IntensityType>(), nsThreshold + 0.5);
  
  // Define which sequence profile to use
  SensitivityProfilePtr sensitivityProfileNSPm 
    = ProfileFactory::createNsProfile(kProbeSequenceLength, options.profileModelRank, options.gstackTuples);
  
  // Read / Compute the profiles 
  if (options.readProfilesFromFile) { // read 
    cout << "Reading Sensitivity Profiles from files in current directory... " << flush;      
    sensitivityProfileNSPm->importFromDatafile(string("SensitivityProfileNsPm-" + filenameSuffix + ".dat"));
  }
  else { 
	vector<Probeset> tmpProbesets;
	vector<PmMmProbePtr> probes;
	ProbesetIntensitiesArray intArray;//, tmpArray;
	vector<IntensityType> logInts;
	vector<IntensityType> newSumLogs;
	for (size_t i = 0; i < probesets.size(); ++i) {
		// @debug
		// All background probes are placed in one pseudo probeset by first filling a vector which is then transformed to a probeset 
		// (otherwise we get biased results, since the original probesets hold probes of similar sequence contant)
		if ((probesets[i].getProbesetId().substr(0, 17) == "BackgroundProbeBG") && controlAntiGenomic){
			for (size_t j = 0; j < probesets[i].getSize(); ++j) {
				probes.push_back(probesets[i].getProbePtr(j)); // Combine all probes of type background probe in one vector 
				logInts.push_back((*mCorrectedIntensities[kPmI])[i][j]); // Save the corresponding intensities
			}
		}
		else if ((probesets[i].getProbesetId().substr(0, 17) == "BackgroundProbeIn") && controlIntron){
			for (size_t j = 0; j < probesets[i].getSize(); ++j) {
				probes.push_back(probesets[i].getProbePtr(j)); // Combine all probes of type background probe in one vector 
				logInts.push_back((*mCorrectedIntensities[kPmI])[i][j]); // Save the corresponding intensities
			}
		}
		else {
			tmpProbesets.push_back(probesets[i]); // Adds the probesets that are no background probesets
			intArray.push_back((*mCorrectedIntensities[kPmI])[i]); // Pushes back the corresponding intensities.
			newSumLogs.push_back(sumLogIs[i]);
		}
	}
	
	vector<Probeset> filteredProbesets;
	ProbesetIntensitiesArray filteredArray;
	
	if (hookN) {
	// Filters the probes that are not background concerning hook.
	/*vector<Probeset> */filteredProbesets = filterParallel(tmpProbesets, newSumLogs, isLessThanIntersectionSumLogI); 
	/*ProbesetIntensitiesArray */filteredArray = filterParallel(intArray, newSumLogs, isLessThanIntersectionSumLogI);
	// Debug: Mark the Probesets that are defined nonspecific by hook (only in the run of the corrected sensitivity measures.)
		if (pmIntensityType == kSensitivityCorrectedPm) {
			vector<Probeset>::iterator pset;
			for (pset = filteredProbesets.begin(); pset != filteredProbesets.end(); ++pset) {
				string id = pset->getProbesetId();
				if (id.substr(0,3) !=  "Bac"){ // Exclude probes that are background in another category (those have been marked before)
					// @debug: This changes the probeset ids to mark the Background probes
//					if(renameBackgroundProbeIds){
//						id = "BackgroundProbeHook" + id;
//					}
					for (size_t p = 0; p < pset->getSize(); ++p){
						pset->getProbePtr(p)->setProbesetId(id);
					}
				}
			}
		}
	}
	
	
// @debug: Add the valarray of background probe intensities (only if controlAntigenomic or controlIntron was set true)
	if (probes.size() > 0) {
		filteredArray.push_back(convertVectorToValarray(logInts));
		Probeset set(probes); // Create a probeset of the background probes
		filteredProbesets.push_back(set); // Add the probeset of background probes.
	}
	
    cout << "\nCalculating NsPm Sensitivity Profile with " << flush;
    
    sensitivityProfileNSPm->computeProfile(filteredProbesets,  // WARNING: Lots of data (probesets) copied!
    										filteredArray,  
    										mem_fun_ref(&PmMmProbe::getSequence),
    										probesetLimit);
    
    // @debug:
    // Write the probesequnces of all background probes to a file (used to analyse base/sequence biases)
/*    ofstream nsProbesequences;
    nsProbesequences.open("nsProbesequences.dat");
    for (size_t i = 0; i < filteredProbesets.size(); ++i) {
    	for (size_t j = 0; j < filteredProbesets[i].getSize(); ++j) {
    		nsProbesequences << filteredProbesets[i].getProbePtr(j)->getSequence() << endl;
    	}
    }*/
    
    
    
//    SimpleSensitivityProfile::exportConditionalSumOfSquares(filteredProbesets, 
//    												 tmpArray, 
//                                                     mem_fun_ref(&PmMmProbe::getSequence),
//                                                     3 // Compute SSE for all triples
//                                                     // options.profileModelRank
//                                                     , "SseProfileBGControlProbes.dat", sensitivityProfileNSPm);
    
  //  cout << "done " << endl;    
  } // Else
  cout << endl;
    
  // Export solution to file and print gnuplot command
  if (!options.readProfilesFromFile){
	cout << "Exporting single profiles" << endl;
    sensitivityProfileNSPm->exportToDatafile(string("SensitivityProfileNsPm-" + filenameSuffix + ".dat"));
    // @bug ? Do not use temporary shared pointers
    sensitivityProfileNSPm->cloneWithNewRank(1)->exportToDatafile(string("SensitivityProfileNsPm-" + filenameSuffix + "-single.dat"));
   }

  // @debug Print sum of squares error
  if (options.isDebugMode) {
    vector<Probeset> probesets = mChip->getProbesets();
    vector<IntensityType> sumLogIs = getAverageSumLogIArray(pmIntensityType, mmIntensityType);
    vector<Probeset> filteredProbesets = filterParallel(probesets, sumLogIs, isLessThanIntersectionSumLogI);

    // Assume these profiles are simple profiles
//    SimpleSensitivityProfilePtr simpleProfileNSPm 
//      = boost::shared_static_cast<SimpleSensitivityProfile>(sensitivityProfileNSPm);
//    SimpleSensitivityProfilePtr simpleProfileNSMm 
//      = boost::shared_static_cast<SimpleSensitivityProfile>(sensitivityProfileNSMm);
    SensitivityProfilePtr simpleProfileNSPm = sensitivityProfileNSPm;
    if (!simpleProfileNSPm) {
      exit(0);
    }

//    cout << "Sum of Squares Error (NsPm - Pm): " 
//         << simpleProfileNSPm->computeSumOfSquares(filteredProbesets, 
//                                                   filterParallel(*mCorrectedIntensities[kPmI], 
//                                                                  sumLogIs, isLessThanIntersectionSumLogI), 
//                                                   mem_fun_ref(&PmMmProbe::getSequence))
//         << endl;
//    cout << "Sum of Squares Error (NsMm - Mm): " 
//         << simpleProfileNSMm->computeSumOfSquares(filteredProbesets, 
//                                                   filterParallel(*mCorrectedIntensities[kMmI], 
//                                                                  sumLogIs, isLessThanIntersectionSumLogI), 
//                                                   mem_fun_ref(&PmMmProbe::getSequenceMm))
//         << endl;
  
    cout << "Export F-Test Sum of Squares Error (NsPm - Pm): " << endl; 
    SimpleSensitivityProfile::exportConditionalSumOfSquares(filteredProbesets, 
                                                     filterParallel(*mCorrectedIntensities[kPmI], 
                                                                    sumLogIs, isLessThanIntersectionSumLogI), 
                                                     mem_fun_ref(&PmMmProbe::getSequence),
                                                     3 // Compute SSE for all triples
                                                     // options.profileModelRank
                                                     , "SseProfile.dat", sensitivityProfileNSPm);
    SimpleSensitivityProfile::exportConditionalSumOfSquares(filteredProbesets, 
                                                     filterParallel(*mCorrectedIntensities[kPmI], 
                                                                    sumLogIs, isLessThanIntersectionSumLogI), 
                                                     mem_fun_ref(&PmMmProbe::getSequence),
                                                     4 // Compute SSE for all triples
                                                     // options.profileModelRank
                                                     , "SseProfile4.dat", sensitivityProfileNSPm);

  }

  // Apply profile
  mHookModel->setProfile(kNsPm, sensitivityProfileNSPm);
  
}











/**
 * Writes the graph of the theoretic function into a file
 * 
 * 
 */
void plotTheoreticCurve(IntensityType a, IntensityType f, IntensityPair kinkCoordinates, IntensityType offset)
{
  ofstream theoreticcurvePlot;
  theoreticcurvePlot.open("theoreticCurve.dat");
  theoreticcurvePlot << "#sum\tdelta\tsumTransf\tdeltaTransf" << endl;
  for (IntensityType rPmLog = -3; rPmLog < 4; rPmLog += 0.03) {
    double rPm = pow(rPmLog, 10);
    IntensityType sumLogI = (log10(rPm + 1.0) - log10(1.0 + f*(rPm + 1.0)))/2.0
      + (log10(rPm*a + 1.0) - log10(1.0 + f*(rPm*a + 1.0)))/2.0;
    IntensityType deltaLogI = log10((rPm + 1.0)/(rPm*a + 1.0)) 
      - log10((1.0 + f*(rPm + 1.0))/(1.0 + f*(rPm*a + 1.0)));
    IntensityType sumLogICorr   = sumLogI + kinkCoordinates.first - offset;
    IntensityType deltaLogICorr = deltaLogI /* + kinkCoordinates.second*/;
    theoreticcurvePlot << sumLogI << "\t" << deltaLogI << "\t" << sumLogICorr << "\t" << deltaLogICorr << endl;
  }
}

/**
 * The function computes updated intensities based on calculationType and writes them
 * to a celfile. It uses the information from a given source celfile for any information 
 * additional to the intensities.
 *
 */
void OligoController::writeIntensitiesToCelfile(const std::string sourceFilename, 
                                                const std::string targetFilename,
                                                const ProbesetIntensitiesArray& intensitiesPm,
                                                const ProbesetIntensitiesArray& intensitiesMm)
{
  // Read original intensity file to overwrite only altered probe intensities
  IntensityMatrixPtr intensityArrayPtr;
  BoolArrayPtr intensityMaskArray;
  MicroarrayExperiment::readArrayFromCelfile(sourceFilename, intensityArrayPtr, intensityMaskArray);

  // Initialize probesets
  const ProbesetVector& probesets = mChip->getProbesets();
  
  // Iterate over every probe in each probeset
  for (size_t probesetIndex = 0; probesetIndex < probesets.size(); ++probesetIndex) {
    for (size_t probeIndex = 0; probeIndex < probesets[probesetIndex].getSize(); ++probeIndex) {
      // Write pm and mm intensity
      std::pair<int,int> positionPm = probesets[probesetIndex].getProbe(probeIndex).getPositionPm();
      (*intensityArrayPtr)[positionPm.first][positionPm.second] = intensitiesPm[probesetIndex][probeIndex];
      std::pair<int,int> positionMm = probesets[probesetIndex].getProbe(probeIndex).getPositionMm();
      (*intensityArrayPtr)[positionMm.first][positionMm.second] = intensitiesMm[probesetIndex][probeIndex];
    }
  }

  // Write to celfile
  MicroarrayExperiment::updateCelfileIntensities(sourceFilename, targetFilename, intensityArrayPtr);
}

// @debug
IntensityType probePmFunction(const PmMmProbe& p, vector<string> mRelevantBasetuples, size_t pos) {
  string sequence = p.getSequence();
  for (std::vector<std::string>::const_iterator tuple = mRelevantBasetuples.begin();
       tuple != mRelevantBasetuples.end(); ++tuple) {
    if (sequence.find(*tuple) == pos) {
//     if (sequence.find(*tuple) != std::string::npos) {
//       return 25000 - 1000*sequence.find(*tuple);
      return 100.0;
    }
  }
  return  1.0;
}



/**
 * Calculate the compression factor parameters
 * @todo: Get rid of the whole compression stuff.
 * 
 * 
 */
void OligoController::calculateCompression(){
	  cout << "Correcting the hookcurve: calculate compression Factor parameters,  " << endl;
	  mHookcurveModel->calculateCompressionFactorParameters(*mCorrectedIntensities[kPmI], mem_fun_ref(&PmMmProbe::getSequence));
}



/**
 * The function computes an abritrary function for each probe and writes the resulting integer  
 * values to a celfile. It uses the information from a given source celfile for any information 
 * additional to the intensities.
 *
 * @note For debug purposes
 */
void OligoController::writeDebuginfoToCelfile(const std::string sourceFilename, 
                                              const std::string targetFilename)
{
  vector<string> tuples; // = options.gstackTuples
  tuples.push_back("GGG");

  // Initialize intensity and mask array
  size_t arrayRowCount = DataCollector::instance().getAs<size_t>("chipRowNumber");
  size_t arrayColumnCount = DataCollector::instance().getAs<size_t>("chipColumnNumber");
 
  // for (size_t i = 0; i <= 22; ++i) {
    IntensityMatrixPtr intensityArrayPtr = IntensityMatrixPtr(new IntensityMatrix(arrayRowCount,arrayColumnCount, 0.0));
  
    // Iterate over every probe in each probeset
    const ProbesetVector& probesets = mChip->getProbesets();
    for (size_t probesetIndex = 0; probesetIndex < probesets.size(); ++probesetIndex) {
      for (size_t probeIndex = 0; probeIndex < probesets[probesetIndex].getSize(); ++probeIndex) {
        // Write pm and mm intensity
        std::pair<int,int> positionPm = probesets[probesetIndex].getProbe(probeIndex).getPositionPm();
        (*intensityArrayPtr)[positionPm.first][positionPm.second] = 
          probePmFunction(probesets[probesetIndex].getProbe(probeIndex), tuples, 0); // i);

        // MM be gray!
        std::pair<int,int> positionMm = probesets[probesetIndex].getProbe(probeIndex).getPositionMm();
        (*intensityArrayPtr)[positionMm.first][positionMm.second] = 
          probePmFunction(probesets[probesetIndex].getProbe(probeIndex), tuples, 0); // i);
        //10; // exp10(intensitiesMm[probesetIndex][probeIndex]);
      }
    }

    // Write to celfile
    string tarFilename = string("UncorrectedGGG0") + boost::lexical_cast<std::string>(0) + ".cel";
    MicroarrayExperiment::updateCelfileIntensities(sourceFilename, tarFilename, intensityArrayPtr);
    //  }
}

/**
 * Sums average intensities over selected tuples and writes them into file
 *
 * @note For debug purposes
 */
void OligoController::writeMotifIntensitiesDebug(const IntensityPredicate considerProbeset, 
                                                 const std::string targetFilename,
                                                 const IntensityMode intensityMode, 
                                                 const PmMmProbeSequenceFunction& getProbeSequence) // mem_fun_ref(&PmMmProbe::getSequence)
{

  vector<IntensityType> sumLogIs = getAverageSumLogIArray(kSensitivityCorrectedPm, kSensitivityCorrectedMm);

  vector<Probeset> filteredProbesets = filterParallel(mChip->getProbesets(), sumLogIs, considerProbeset);
  ProbesetIntensitiesArray filteredLogIs = filterParallel(*mCorrectedIntensities[intensityMode], sumLogIs, considerProbeset);

  // vector<Probeset> filteredProbesets = mChip->getProbesets();
  // ProbesetIntensitiesArray filteredLogIs = *mCorrectedIntensities[intensityMode];


  // string motifs[] = {"GG", "GGG", "GGGG", "CC", "CCC", "CCCC", "AA", "AAA", "AAAA", "TT", "TTT", "TTTT", "GGG1", "CCC12", "N"};
  string motifs[] = {"GG", "GGG", "GGGG", "CC", "CCC", "CCCC", "CCTCCC",  "CCGCCT",  "CCCTCC",  "TCCGCC", "AAA", "TTT", "GGG1", "CCC12", "N"};
  
  const size_t motifCount = 3*4 + 3;
  
  vector<size_t> motifCounter(motifCount, 0);
  vector<IntensityType> motifLogISum(motifCount, 0.0);
  

  for (size_t probesetIndex = 0; probesetIndex < filteredProbesets.size(); ++probesetIndex) {
    for (size_t probeIndex = 0; probeIndex < filteredProbesets[probesetIndex].getSize(); ++probeIndex) {
      string sequence = getProbeSequence(filteredProbesets[probesetIndex].getProbe(probeIndex));
      IntensityType logI = filteredLogIs[probesetIndex][probeIndex];
      
      for (size_t i = 0; i < 3*4; ++i) {
        if (sequence.find(motifs[i]) != std::string::npos) {
          motifCounter[i] = motifCounter[i] + 1;
          motifLogISum[i] = motifLogISum[i] + logI;
        }
      }

      size_t i = 3*4 +0;
      if (sequence.find("GGG") == 1) {
        motifCounter[i] = motifCounter[i] + 1;
        motifLogISum[i] = motifLogISum[i] + logI;        
      }
      i = 13;
      if (sequence.find("CCC") == 12) {
        motifCounter[i] = motifCounter[i] + 1;
        motifLogISum[i] = motifLogISum[i] + logI;        
      }
      i = 14;
      motifCounter[i] = motifCounter[i] + 1;
      motifLogISum[i] = motifLogISum[i] + logI;                          
    }
  }

  ofstream oFile;
  oFile.open(targetFilename.c_str());
  for (size_t i = 0; i < 3*4 +3; ++i) {
    oFile << motifs[i] << "\t" << motifCounter[i] << "\t" << motifLogISum[i] << "\t" << (motifLogISum[i]/(IntensityType)motifCounter[i]) << endl;
  }
  oFile.close(); 
}




/**
 * Estimates the parameters a,f of the theoretic function and sets new Imax = 1/f
 *
 *
 */
boost::tuple<IntensityType, IntensityType>  OligoController::computeImaxAndA(const IntensityMappingPtr hookcurvePlot)
{
	IntensityPair intersectionCoordinates = computeKinkCoordinates("", MathUtil::digitizeCurveNew(*hookcurvePlot), mHookModel);
	
  // Subtract max(NS) and log(b)
  transform(hookcurvePlot->begin(), hookcurvePlot->end(), hookcurvePlot->begin(),
            AddToPair<IntensityPair>(IntensityPair(-intersectionCoordinates.first, 0 /*- intersectionCoordinates.second*/)));
  
  cout << "Calculating theoretic curve." << endl;
  IntensityType a, f, offset;
  tie(a,f, offset) = mHookcurveModel->estimateTheoreticFunctionParameters(*hookcurvePlot);

  // Dump plots to a file
  plotTheoreticCurve(a, f, intersectionCoordinates, offset);
//   printVectorPairsToFile(mHookcurveModel->calculateProbesetR(intersectionCoordinates.first, a, f, offset),
//                          "theoreticCurveProbesetR.dat", "#");

  // Now correct Imax and the intensities of the probes:
  if ((log10(1/f) + intersectionCoordinates.first) > mHookModel->getImax()) {
    cout << "Setting Imax from " << mHookModel->getImax() << " to " << flush;
    mHookModel->setImax(log10(1/f) + intersectionCoordinates.first);
    cout << mHookModel->getImax() << endl;
  }
  else {
    cout << "Calculated Imax from theoretic curve is " << (log10(1/f) + intersectionCoordinates.first) 
         << " since it is smaller than the pre Imax (" << mHookModel->getImax() << "). We keep the old one." << endl;
  }
  //cout << " f: " << f << " new Imax: " << log10(1/f) + intersectionCoordinates.first << endl;// << " a " << a << endl;
  return make_tuple(f,a);
}




/**
 * This function coordinates the subtraction of the nonspecific background from each intensity
 * 
 * @todo Clean up this messy function
 */
void OligoController::subtractNBackgroundPmOnly(const IntensityType nsThreshold)
{
	
  IntensityPredicate isLessThanIntersectionSumLogI = bind2nd(less<IntensityType>(), nsThreshold);
  vector<IntensityType> sumLogIs = getAverageSumLogIArray(kSensitivityCorrectedPm, kSensitivityCorrectedPseudoMm);

  DistributionParameters pmDistrib;
    
  // Get the distribution of the nonspecific probe intensities.
  vector<IntensityType>  intensitiesVectorPm = valarrayVectorFlat(filterParallel(*(mCorrectedIntensities[kSensitivityCorrectedPm]), sumLogIs, isLessThanIntersectionSumLogI));
  pmDistrib.mean    = MathUtil::calculateMean(intensitiesVectorPm);
  pmDistrib.sigma   = MathUtil::calculateStandardDeviation(intensitiesVectorPm);
    
  cout << "\nMean Distribution_NsPM: " << pmDistrib.mean << " sigma: " << pmDistrib.sigma << endl;
    
  // Subtract the nonspecific background and write a new ProbesetIntensityArray entry.
  mCorrectedIntensities[kNsSubstractedPm] = Probeset::initializeIntensitiesArray(mChip->getProbesets());
  mCorrectedIntensities[kNsSubstractedPm] 
    =  mHookcurveModel->calculateSpecificPortionsMonovariate(*mCorrectedIntensities[kPmI], /**mCorrectedIntensities[kSensitivityCorrectedPm],*/ 
                                                             *mHookModel->getProfile(kNsPm),
                                                             pmDistrib);  
}




/**
 * This function coordinates the subtraction of the nonspecific background from each intensity
 * 
 * @todo Clean up this messy function
 */
void OligoController::subtractNBackground(const IntensityType nsThreshold)
{
	
	IntensityPredicate isLessThanIntersectionSumLogI = bind2nd(less<IntensityType>(), nsThreshold);
  vector<IntensityType> sumLogIs = getAverageSumLogIArray(kSensitivityCorrectedPm, kMmI); //getAverageSumLogIArray(kSensitivityCorrectedPm, kSensitivityCorrectedMm);

  DistributionParameters pmDistrib, mmDistrib;
  BivariateDistributionParameters bivariateDistributionParameters;
    
  vector<IntensityType>  intensitiesVectorPm = valarrayVectorFlat(filterParallel(*(mCorrectedIntensities[kSensitivityCorrectedPm]), sumLogIs, isLessThanIntersectionSumLogI));
  pmDistrib.mean    = MathUtil::calculateMean(intensitiesVectorPm);
  pmDistrib.sigma   = MathUtil::calculateStandardDeviation(intensitiesVectorPm);
  bivariateDistributionParameters.meanPm = pmDistrib.mean;
  bivariateDistributionParameters.sigmaPm = pmDistrib.sigma;
    
  // Same for MMs
  vector<IntensityType>  intensitiesVectorMm = valarrayVectorFlat(filterParallel(*(mCorrectedIntensities[kSensitivityCorrectedMm]), sumLogIs, isLessThanIntersectionSumLogI));
  mmDistrib.mean    = MathUtil::calculateMean(intensitiesVectorMm);
  mmDistrib.sigma   = MathUtil::calculateStandardDeviation(intensitiesVectorMm);
    
  // The correlation of PMs and MMs is needed for the Diff method.
  bivariateDistributionParameters.correlation = MathUtil::calculateCorrelation(intensitiesVectorPm, intensitiesVectorMm);
    
    
  bivariateDistributionParameters.sigmaAverage =    (pmDistrib.sigma + mmDistrib.sigma) / 2.0;  /// @todo: This is wrong (but was calculated like this before... Use the following way instead:
  //IntensityType sigmaAverage =    MathUtil::calculateStandardDeviation(intensitiesVectorPm + intensitiesVectorMm);
  bivariateDistributionParameters.logB =   mHookModel->getParam("LogB"); // MathUtil::calculateMean(intensitiesVectorPm - intensitiesVectorMm);
    
    
  cout << "\nMean Distribution_NsPM: " << pmDistrib.mean << " sigma: " << pmDistrib.sigma << endl;
  cout << "Mean Distribution_NsMM: " << mmDistrib.mean << " sigma: " << mmDistrib.sigma << endl;
  cout << "Correlation: " << bivariateDistributionParameters.correlation << endl;
  cout << "logB  " << bivariateDistributionParameters.logB << endl; 
    
    
  // Subtract the nonspecific background and write a new ProbesetIntensityArray entry.
  // @todo: Only calculate those that are needed.
  mCorrectedIntensities[kNsSubstractedPm] = Probeset::initializeIntensitiesArray(mChip->getProbesets());
  mCorrectedIntensities[kNsSubstractedPm] 
    =  mHookcurveModel->calculateSpecificPortionsMonovariate(*mCorrectedIntensities[kPmI], /**mCorrectedIntensities[kSensitivityCorrectedPm],*/ 
                                                             *mHookModel->getProfile(kNsPm),
                                                             pmDistrib);
  mCorrectedIntensities[kNsSubstractedMm] = Probeset::initializeIntensitiesArray(mChip->getProbesets());
  mCorrectedIntensities[kNsSubstractedMm] 
    =  mHookcurveModel->calculateSpecificPortionsMonovariate(*mCorrectedIntensities[kMmI], /**mCorrectedIntensities[kSensitivityCorrectedMm],*/ 
                                                             *mHookModel->getProfile(kNsMm),
                                                             mmDistrib, 
                                                             mem_fun_ref(&PmMmProbe::getSequenceMm));
    
  mCorrectedIntensities[kNsSubstractedFastDiff] = Probeset::initializeIntensitiesArray(mChip->getProbesets());
  mCorrectedIntensities[kNsSubstractedFastDiff] 
    =  mHookcurveModel->calculateSpecificPortionsBivariate(*mCorrectedIntensities[kPmI], *mCorrectedIntensities[kMmI],
                                                           *mHookModel->getProfile(kNsPm), *mHookModel->getProfile(kNsMm),
                                                           bivariateDistributionParameters,
                                                           &MathUtil::simpleBivariateIntegralFunction);
    
  mCorrectedIntensities[kNsSubstractedGcrmaDiff] = Probeset::initializeIntensitiesArray(mChip->getProbesets());
  mCorrectedIntensities[kNsSubstractedGcrmaDiff] 
    =  mHookcurveModel->calculateSpecificPortionsBivariate(*mCorrectedIntensities[kPmI], *mCorrectedIntensities[kMmI], 
                                                           /**mCorrectedIntensities[kSensitivityCorrectedPm], *mCorrectedIntensities[kSensitivityCorrectedMm],  */
                                                           *mHookModel->getProfile(kNsPm), *mHookModel->getProfile(kNsMm),
                                                           bivariateDistributionParameters,
                                                           &MathUtil::gcrmaLikeBivariateIntegralFunction);  
}


// We expect the iTs to be sorted.
double getW(const valarray<double>& iTs, const valarray<double>& ws, double iT)
{
	for (size_t i = 0; i < iTs.size(); ++i) {
		if (iT < iTs[i]) {
			cout << iT << "  !=  " << iTs[i] << endl;
			return ws[i];
		}
	}
	return ws[ws.size()-1];
}


void OligoController::inverseWashingOfIntensities()
{
	//ProbesetIntensitiesArray& intensities = mCorrectedIntensities[kPm];
	mCorrectedIntensities[kWashedPm] = Probeset::initializeIntensitiesArray(mChip->getProbesets());
	// Washing function constants: (as in "Washing scaling of microarray expression" P24 line 5)
	double wMax    	= 0.9;
	double wMin      	= 0.06;
	double gamma   	= 1.6;
	double aPrime    	= 0.1;
	double iMaxLog 	= 4.5; //mHookModel->getImax();
	cout << "Imax: " << iMaxLog << endl;
	double iMax         = pow(10,iMaxLog);
	//double plotMax = iMax+1;
	size_t numberOfTableEntries = 1000;
	double stepSize = iMaxLog / numberOfTableEntries; 
	valarray<double> i_0(0.0, numberOfTableEntries);
	valarray<double> wI_0(0.0, numberOfTableEntries);
	valarray<double> iT(0.0, numberOfTableEntries);
	//cout << "Imax " << iMax << ", step: " << stepSize << endl;
	for (size_t i = 0 ; i < numberOfTableEntries; ++i) {
		double i0Log = i * stepSize; // steps are chosen to span the log scale in equidistant steps.
		double i0 = pow(10, i0Log); // This way we have a higher density for smaller values.
		i_0[i]     = i0;
		wI_0[i]  = exp( -1 * pow((aPrime * (( iMax - i0) / i0)), gamma) ) * (wMax - wMin) + wMin; // As in the paper
		iT[i]     = pow(10, i0Log + log10(wI_0[i]));
		//cout << iT[i] << endl;
		//cout << i_0[i] << "     " << wI_0[i] << "     " << iT[i] << endl;
	}
	//exit(0);
	
   // Iterate over every probe in each probeset
   for (size_t probesetIndex = 0; probesetIndex < (*mCorrectedIntensities[kPm]).size(); ++probesetIndex) {
     for (size_t probeIndex = 0; probeIndex < (*mCorrectedIntensities[kPm])[probesetIndex].size(); ++probeIndex) {
       // Write pm and mm intensity
    	 (*mCorrectedIntensities[kWashedPm])[probesetIndex][probeIndex] = getW(iT, wI_0, pow(10, (*mCorrectedIntensities[kPm])[probesetIndex][probeIndex]));
     }
   }
	
}



/**
 * Organizes the calculation of the expression measures. Only those expression measures are calculated that have been chosen by the command line parameter e
 * 
 * For each chosen expression measure a probeset avergaed measure is printed as well as a independent measures for each probe.
 * The inverse Logarithm is used to avoid negative Expression Measures. for probesets also a inverse Glog expression measure is printedd on probeset level.
 * 
 * 
 */
void OligoController::calculateAndPrintExpressionMeasures(ExpressionMeasurePtr exprCalculater)
{
	  UnaryIntensityFunction invLogFunction = (IntensityTypeFunction)exp10;
	  UnaryIntensityFunction invGLogFunction = &MathUtil::gExp10;
	
	  if (options.expressionMeasureIdentifier[0] == '1') {
	   cout << "\n\nCalculating Expression Measures by PM... " << endl;
	   // Pm Expression Measure with invGlog on probeset scale
	   printVectorPairsToFile(exprCalculater->calculateProbesetExpressionMeasure
                            (*mCorrectedIntensities[kNsSubstractedPm], kSPm,
	                          mem_fun_ref(&PmMmProbe::getSequence), invGLogFunction), 
	                          "ExpressionMeasure-Pm.dat", "#\t" + options.intensityFilename);
	   // Pm Expression Measure with invLog on probeset scale and printing results to the kEmPmSingleIntegrationExp10 ProbesetIntensitiesArray
	   printVectorPairsToFile(exprCalculater->calculateProbesetExpressionMeasure
                            (*mCorrectedIntensities[kNsSubstractedPm], kSPm,
                             mem_fun_ref(&PmMmProbe::getSequence), 
                             invLogFunction, kEmPmSingleIntegrationExp10), // Printed to array
	                          "ExpressionMeasure-PmExp10.dat", "#\t" + options.intensityFilename);
	   // Pm Expression Measure with invLog on probe scale
	   printVectorPairsToFile(exprCalculater->calculateProbeExpressionMeasures
                            (*mCorrectedIntensities[kNsSubstractedPm], kSPm, 
	                          mem_fun_ref(&PmMmProbe::getSequence), 
                             invLogFunction, kEmPmSingleIntegrationExp10_Probe),
	                          "ExpressionMeasurePmExp10_Probe.dat", "#\t" + options.intensityFilename);
	   
	  }
	  if (options.expressionMeasureIdentifier[1] == '1') {
	   cout << "Calculating Expression Measures by MM... " << endl;
	   // Mm Expression Measure with invGLog on probeset scale
	   printVectorPairsToFile(exprCalculater->calculateProbesetExpressionMeasure
                            (*mCorrectedIntensities[kNsSubstractedMm], kSMm,
	                          mem_fun_ref(&PmMmProbe::getSequenceMm), invGLogFunction), 
	                          "ExpressionMeasure-Mm.dat", "#\t" + options.intensityFilename);
	   // Mm Expression Measure with invLog on probeset scale
	   printVectorPairsToFile(exprCalculater->calculateProbesetExpressionMeasure
                            (*mCorrectedIntensities[kNsSubstractedMm], kSMm,
                             mem_fun_ref(&PmMmProbe::getSequenceMm), 
                             invLogFunction, kEmMmSingleIntegrationExp10), 
	                          "ExpressionMeasure-MmExp10.dat", "#\t" + options.intensityFilename);
	   // Mm Expression Measure with invGLog on probe scale
	   printVectorPairsToFile(exprCalculater->calculateProbeExpressionMeasures
                            (*mCorrectedIntensities[kNsSubstractedMm], kSMm, 
                             mem_fun_ref(&PmMmProbe::getSequenceMm), 
                             invLogFunction, kEmMmSingleIntegrationExp10_Probe),
	                          "ExpressionMeasureMmExp10_Probe.dat", "#\t" + options.intensityFilename);
	  }
	  if (options.expressionMeasureIdentifier[2] == '1') {
	   cout << "Calculating Expression Measures by FastDiff... " << endl;
	   // Delta single integration Expression Measure with invGLog on probeset scale
	   printVectorPairsToFile(exprCalculater->calculateProbesetExpressionMeasure
                            (*mCorrectedIntensities[kNsSubstractedFastDiff], kSFastDiff, 
	                          mem_fun_ref(&PmMmProbe::getSequence), invGLogFunction),
	                          "ExpressionMeasureFastDiff_Probeset.dat", "#\t" + options.intensityFilename);
	   // Delta single integration Expression Measure with invLog on probeset scale
	   printVectorPairsToFile(exprCalculater->calculateProbesetExpressionMeasure
                            (*mCorrectedIntensities[kNsSubstractedFastDiff], kSFastDiff, 
	                          mem_fun_ref(&PmMmProbe::getSequence), 
                             invLogFunction, kEmDeltaSingleIntegrationExp10), // Printed to array
	                          "ExpressionMeasure-FastDiffExp10.dat", "#\t" + options.intensityFilename);
	   // Delta single integration Expression Measure with invLog on probe scale
	   printVectorPairsToFile(exprCalculater->calculateProbeExpressionMeasures
                            (*mCorrectedIntensities[kNsSubstractedFastDiff], kSFastDiff, 
	                          mem_fun_ref(&PmMmProbe::getSequence), 
                             invLogFunction, kEmDeltaSingleIntegrationExp10_Probe),
	                          "ExpressionMeasureFastDiffExp10_Probe.dat", "#\t" + options.intensityFilename);
	  }
	  if (options.expressionMeasureIdentifier[3] == '1') {
	   cout << "Calculating Expression Measures by GCRMA like Diff... " << endl;
	   // GCRMA like Delta single integration Expression Measure with invGlog on probeset scale
	   printVectorPairsToFile(exprCalculater->calculateProbesetExpressionMeasure
                            (*mCorrectedIntensities[kNsSubstractedGcrmaDiff], kSGcrmaDiff, 
	                          mem_fun_ref(&PmMmProbe::getSequence), invGLogFunction), 
	                          "ExpressionMeasureGCRMADiff_Probeset.dat", "#\t" + options.intensityFilename);
	   // GCRMA like Delta single integration Expression Measure with invLog on probeset scale
	   printVectorPairsToFile(exprCalculater->calculateProbesetExpressionMeasure
                            (*mCorrectedIntensities[kNsSubstractedGcrmaDiff], kSGcrmaDiff, 
	                          mem_fun_ref(&PmMmProbe::getSequence), 
                             invLogFunction, kEmDeltaGcrmaLikeIntegrationExp10), // Printed to array
	                          "ExpressionMeasure-GCRMADiffExp10.dat", "#\t" + options.intensityFilename);
	   // GCRMA like Delta single integration Expression Measure with invLog on probe scale
	   printVectorPairsToFile(exprCalculater->calculateProbeExpressionMeasures
                            (*mCorrectedIntensities[kNsSubstractedGcrmaDiff], kSGcrmaDiff, 
	                          mem_fun_ref(&PmMmProbe::getSequence), 
                             invLogFunction, kEmDeltaGcrmaLikeIntegrationExp10_Probe), 
	                          "ExpressionMeasureGCRMADiffExp10_Probe.dat", "#\t" + options.intensityFilename);
	  }

}



void OligoController::logInitialOptions() 
{
  // Log chip id (file basename without extension)
  // @note Obtaining chip-id from intensity-filename may not be general enough
  string chipname = options.intensityFilename.substr(0, options.intensityFilename.length() - 4); // remove ".cel"
  size_t basenameIndex = chipname.rfind("/");
  if (basenameIndex == string::npos) { // @bug Due to windows compatibility, filenames must not contain "\"
    basenameIndex = chipname.rfind("\\");
  }
  DataCollector::instance().insert("chipId", chipname.substr(basenameIndex + 1));

  // Log starting time of computation
  time_t timeRaw;
  time(&timeRaw);
  DataCollector::instance().insert("startingTime", ctime(&timeRaw));
}

/**
 * Returns the affinity profile for non-specific binding according to parameters
 *
 * @param sequenceLength Length of the sequence.
 * @param modelRank Size of the base tuples within the current model
 * (mononucleotides = 1, dinucleotides = 2, ...)
 * @param relevantBasetuples The base tuples (e.g. ["GGGG","CCCC"]) 
 *        the contributions of which the profile shall model 
 *
 * @return new sequence profile
 */
SensitivityProfilePtr ProfileFactory::createNsProfile(const size_t sequenceLength, 
                                                    const unsigned short modelRank, 
                                                    const std::vector<std::string> gstackTuples)
{
  SensitivityProfilePtr returnProfile; // avoid temporary shared pointers

  if (gstackTuples.empty()) {  // Return simple profiles
    returnProfile = SensitivityProfilePtr(new SimpleSensitivityProfile(sequenceLength, modelRank));  
  }
  else { // use composite profiles (simple and gstacks)
    returnProfile = SensitivityProfilePtr(new AdditivePofile(gstackTuples, sequenceLength, modelRank));
  }
  return returnProfile;
}

/**
 * Returns the respective affinity profile for specific binding according to parameters.
 * Contrary to the non-specific profile, the gstacks-option is omitted (it is known that such a 
 * effect doesnot occur in specific binding) and the model-rank is at most 2 (nearest-neighbor
 * model).
 *
 * @return new sequence profile
 */
SensitivityProfilePtr 
ProfileFactory::createSpecificProfile(const size_t sequenceLength, const unsigned short modelRank, 
                                      const std::vector<std::string> gstackTuples)
{
  SensitivityProfilePtr returnProfile =  // avoid temporary shared pointers
    ProfileFactory::createNsProfile(sequenceLength, min((int) modelRank, 2));
    // ProfileFactory::createNsProfile(sequenceLength, modelRank, gstackTuples);
  return returnProfile;
}

size_t ProfileFactory::getMinimumProbesetCountToEstimateModel(const size_t modelRank) {
  if (modelRank == 1) // N
    return kMinimumSpecificProbesetsForN;
  if (modelRank == 2) // NN
    return kMinimumSpecificProbesetsForNn;
  return 5000; // otherwise
}
