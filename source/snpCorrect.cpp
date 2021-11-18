/**
 * @file snpCorrect.cpp -
 * @author $Author: mario $ 
 * @date $Date: 2008-08-22 12:54:40 +0200 (Fri, 22 Aug 2008) $
 * 
 * The main program for SNP-chip analysis
 */
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator> 

#include "boost/tuple/tuple.hpp"
#include <boost/lexical_cast.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/function.hpp>

#include "DataCollector.hpp"
#include "MicroarrayExperiment.hpp"
#include "Probeset.hpp"
#include "HookcurveAnalyzer.hpp"
#include "StringUtil.hpp"
#include "MathUtil.hpp"
#include "BackgroundSubtraction.hpp"
#include "StlUtil.hpp"
#include "ValarrayUtil.hpp"
#include "ProbeFilter.hpp"
#include "ProbesetComposition.hpp"
#include "SensitivityProfile.hpp"
#include "DetectKinkPoint.hpp"
#include "OligoController.hpp"
#include "HookcurveStatistics.hpp"
#include "Chip.hpp"
#include "HookModel.hpp"

 
#include <functional>
#include <numeric>
#include <valarray>


using namespace std;
using namespace larrpack;
using namespace stringutil;
using namespace boost;
using boost::lexical_cast;

// All data structures needed for an independent hookcurve analysis 
struct HaplotypeData {
  PmMmProbePtrVector probes;
  ChipPtr chip;
  HookModelPtr model;
  HookcurveAnalyzerPtr analyzer;
  IntensityMappingPtr hookcurvePlot;
};

void plotTheoreticCurve(IntensityType a, IntensityType f, IntensityPair kinkCoordinates, IntensityType offset, std::string filename)
{
  ofstream theoreticcurvePlot;
  theoreticcurvePlot.open(filename.c_str());
  theoreticcurvePlot << "#sum\tdelta\tsumTransf\tdeltaTransf" << endl;
  for (IntensityType rPmLog = -3; rPmLog < 4; rPmLog += 0.03) {
    double rPm = pow(rPmLog, 10);
    IntensityType sumLogI = (log10(rPm + 1.0) - log10(1.0 + f*(rPm + 1.0)))/2.0
      + (log10(rPm*a + 1.0) - log10(1.0 + f*(rPm*a + 1.0)))/2.0;
    IntensityType deltaLogI = log10((rPm + 1.0)/(rPm*a + 1.0)) 
      - log10((1.0 + f*(rPm + 1.0))/(1.0 + f*(rPm*a + 1.0)));
    IntensityType sumLogICorr   = sumLogI + kinkCoordinates.first - offset;
    IntensityType deltaLogICorr = deltaLogI + kinkCoordinates.second;
    theoreticcurvePlot << sumLogI << "\t" << deltaLogI << "\t" << sumLogICorr << "\t" << deltaLogICorr << endl;
  }
}

IntensityArray computeTripleContributions(SnpProbesetVector& probesets,
                                          ProbesetIntensitiesArrayPtr intensities,
                                          SimpleSensitivityProfilePtr simpleProfile) {
  size_t baseTupleCount = (size_t)pow((float)sequenceutil::kNucleotideCount, 3);
  IntensityArray tripleProfile(0.0, baseTupleCount);
  IntensityArray motifCount(0.0, baseTupleCount);

  size_t validSequenceIndices = 23;
  IntensityArray tripleProfile2(0.0, baseTupleCount* validSequenceIndices);
  IntensityArray motifCount2(0.0, baseTupleCount * validSequenceIndices);


  // For each probeset
  for (size_t probesetIndex = 0; probesetIndex < probesets.size(); ++probesetIndex) { // 3; ++probesetIndex) { //
    
    const SnpProbeset& currentProbeset = probesets[probesetIndex];

    // Calculate corrected intensity
    // Should i do this outside and pass values?
    // Actually it has proven well to update intensities once, and store array
    
    // First, calculate increment: analog to HookModel::calculateIncrements, not sure if we can use HookModel
    std::vector<std::string> sequences = currentProbeset.getProbeSequences();
    IntensityArray increment(0.0, currentProbeset.getSize());
    for (size_t probeIndex = 0; probeIndex < currentProbeset.getSize(); ++probeIndex) {            
      increment[probeIndex] = simpleProfile->getSequenceIncrement(sequences[probeIndex]);
    }
    
    IntensityArray correctedIntensities = (*intensities)[probesetIndex]; // - increment;
    // cout << "uncorrectedIntensities: " << correctedIntensities << endl;    
    correctedIntensities = correctedIntensities - increment;
    // cout << "correctedIntensities: " << correctedIntensities << endl;

    // Compute contribution due to base triple flanking offset
    IntensityArray tripleDecrement(0.0, currentProbeset.getSize());
    std::valarray<int> offsets(0, currentProbeset.getSize());
    transform(currentProbeset.beginProbe(), currentProbeset.endProbe(), &offsets[0],
              boost::mem_fn(&SnpProbe::getLocation));    
    // cout << "offsets: " << offsets << endl;
    for (size_t probeIndex = 0; probeIndex < currentProbeset.getSize(); ++probeIndex) {
      int l = 12 + offsets[probeIndex]; // probe sequence index 13 + offset
      tripleDecrement[probeIndex] = simpleProfile->getSubTripleIncrement(sequences[probeIndex].substr(l-1, 3), l);
    }

    // Compute practical increment including triple contribution
    IntensityType average = computeAverageIntensity(correctedIntensities);
    // cout << "average: " << average << endl;

    // @frage Vor oder nach dem Abzug des triples mitteln??  
    // IntensityArray deltaA = correctedIntensities /* + tripleDecrement */ - average;
    IntensityArray deltaA = correctedIntensities + tripleDecrement - average;
    // IntensityArray deltaA = correctedIntensities;// + tripleDecrement - average;

    // cout << "deltaA: " << deltaA << endl;

    // Now compute position-independent triple averages
    for (size_t probeIndex = 0; probeIndex < currentProbeset.getSize(); ++probeIndex) {
      int l = 12 + offsets[probeIndex]; // probe sequence index 13 + offset
      size_t tripletIndex = sequenceutil::sequenceToIndex(sequences[probeIndex].substr(l-1, 3));
      tripleProfile[tripletIndex] += deltaA[probeIndex];
      motifCount[tripletIndex] += 1;

      tripleProfile2[tripletIndex * validSequenceIndices + l] += deltaA[probeIndex];
      motifCount2[tripletIndex * validSequenceIndices + l] += 1;     

      // cout << "sequences[probeIndex]: " << sequences[probeIndex] << endl;
      // cout << "sequences[probeIndex].substr(l-1, 3): " << sequences[probeIndex].substr(l-1, 3) << endl;
    }   
    // cout << "motifCount2: " << motifCount2 << endl;
  }
   
  return tripleProfile/motifCount;
  // return tripleProfile2 / motifCount2;
}

int main (int argc, char *argv[])
{  
  // Read options from File/Commandline
  ProgramOptions options = parseCommandlineOptions(argc, argv); 
  OligoController controller(options);
  controller.logInitialOptions();
    
  cout << "Importing Experiment and applying background correction" << endl;


  // Does not work, as it creates PmMmProbes
  // PmMmProbePtrVector probes = controller.importProbeData(); 

  // ------------------------------------------------------------------------------------------
  // Following part is copied from OligoController.importProbeData
  // @todo Rewrite such that function importProbeData can give back SnpProbes!
  // ------------------------------------------------------------------------------------------
  MicroarrayExperiment microExp;
  IntensityMatrixPtr intensityArray = IntensityMatrixPtr();
  BoolArrayPtr intensityMaskArray = BoolArrayPtr();
  MicroarrayExperiment::readArrayFromCelfile(options.intensityFilename, intensityArray, intensityMaskArray);
  DataCollector::instance().insert("chipRowNumber", intensityArray->rowNum());
  DataCollector::instance().insert("chipColumnNumber", intensityArray->columnNum());
  microExp.readFromSnpSequenceFile(options.chipDiscriptionFilename, *intensityArray);
  SnpProbePtrVector probes = extractSnpProbesFromProbelist(microExp.getProbes()); 

  // // Read from CRLMM
  // map<string,char>  initialGenotypes = ComposeProbesetsByGenotypeCall::readGenotypeCalls
  //   (options.genotypeCallFilename, DataCollector::instance().getAs<string>("chipId") + ".CEL");
  // Read from Affy GDAS
  map<string,char>  initialGenotypes = ComposeProbesetsByGenotypeCall::readGenotypeCallsFromGDAS
    (options.genotypeCallFilename, DataCollector::instance().getAs<string>("chipId"));

  SnpProbesetVector probesetsX = ComposeProbesetsByGenotypeCall::composeSnpProbesets(probes, initialGenotypes);

  SensitivityProfilePtr sensitivityProfileNSPm 
    = ProfileFactory::createNsProfile(kProbeSequenceLength, options.profileModelRank, options.gstackTuples);


  // Write all probe data to file
  if (1) {
    ofstream probeInfoFile;
    probeInfoFile.open("CompleteProbeInfo.dat");

    probeInfoFile << "Probe Set ID\t probe x\t probe y\t probe interrogation position\t probe sequence\t target strandedness\t Chromosome\t Call\t PM" << endl;

    SnpProbesetVector probesets = probesetsX;

    // Loop over all probesets
    for (size_t probesetIndex = 0; probesetIndex < probesets.size(); ++probesetIndex) {
      SnpProbeset& currentProbeset = probesets[probesetIndex];
      for (size_t probeIndex = 0; probeIndex < currentProbeset.getSize(); ++probeIndex) {
        SnpProbe& currentProbe = *currentProbeset.getProbePtr(probeIndex);
        probeInfoFile << currentProbe.getProbesetId() << "\t" 
                      << currentProbe.getPositionPm().first << "\t" 
                      << currentProbe.getPositionPm().second << "\t" 
                      << currentProbe.getLocation() << "\t" 
                      << currentProbe.getSequence() << "\t" 
                      << currentProbe.getStrand() << "\t" 
                      << currentProbe.getChromosome() << "\t" 
                      << currentProbeset.getCall() << "\t" 
                      << currentProbe.getPm()
                      << endl;
      }
      currentProbeset.swapCross();
      for (size_t probeIndex = 0; probeIndex < currentProbeset.getSize(); ++probeIndex) {
        SnpProbe& currentProbe = *currentProbeset.getProbePtr(probeIndex);
        probeInfoFile << currentProbe.getProbesetId() << "\t" 
                      << currentProbe.getPositionPm().first << "\t" 
                      << currentProbe.getPositionPm().second << "\t" 
                      << currentProbe.getLocation() << "\t" 
                      << currentProbe.getSequence() << "\t" 
                      << currentProbe.getStrand() << "\t" 
                      << currentProbe.getChromosome() << "\t" 
                      << currentProbeset.getCall() << "\t" 
                      << currentProbe.getPm()
                      << endl;
      }
      currentProbeset.swapCross();
    }

    probeInfoFile.close();
    return(0);
  }

  if (1) {
    // @test: Use only probesets with first allele present
    // @note could be done with STL copy_if
    SnpProbesetVector probesets;
    for (size_t probesetIndex = 0; probesetIndex < probesetsX.size(); ++probesetIndex) {
      const SnpProbeset& currentProbeset = probesetsX[probesetIndex];
      //if (currentProbeset.getCall() == 'A') {
      if (currentProbeset.getCall() != 'H') {
        probesets.push_back(currentProbeset);
      }
    }
 
    // ------------------------------------------------------------------------------------------
    // Following part resembles OligoController.initializeIntensityArrays
    // ------------------------------------------------------------------------------------------
    // Initialize corrected intensities

  ProbesetIntensitiesArrayPtr intensitiesAllele = SnpProbeset::initializeIntensitiesArray(probesets);

  // Loop over all probesets 
  for (size_t probesetIndex = 0; probesetIndex < probesets.size(); ++probesetIndex) {
    const SnpProbeset& currentProbeset = probesets[probesetIndex];
    //if (currentProbeset.getCall() != 'H') {
      for (size_t probeIndex = 0; probeIndex < currentProbeset.getSize(); ++probeIndex) {            
        SnpProbe& currentProbe = *currentProbeset.getProbePtr(probeIndex);
        (*intensitiesAllele)[probesetIndex][probeIndex] = log10(currentProbe.getPm());
      }
      //}
  }

  // Create list of probe sequences
  ProbesetSequencesArray probesetSequences;
  for (SnpProbesetVector::const_iterator pi = probesets.begin(); pi != probesets.end(); ++pi) {
    probesetSequences.push_back(pi->getProbeSequences());
  } 

  cout << "options.probesetLimitForProfileComputation: " << options.probesetLimitForProfileComputation << endl;
  cout << "\nCalculating NsPm Sensitivity Profile with " << flush;    
  sensitivityProfileNSPm->computeProfile(probesetSequences, *intensitiesAllele, 
                                         options.probesetLimitForProfileComputation); 
  sensitivityProfileNSPm->exportToDatafile(string("SensitivityProfileNsPm.dat"));

  SimpleSensitivityProfilePtr simpleProfileNSPm 
    = boost::shared_static_cast<SimpleSensitivityProfile>(sensitivityProfileNSPm);
  // cout << simpleProfileNSPm->getSubTripleIncrement("GGG", 0) << endl;
  // cout << simpleProfileNSPm->getSubTripleIncrement("GGG", 1) << endl;
  // cout << simpleProfileNSPm->getSubTripleIncrement("GGG", 21) << endl;
  // cout << simpleProfileNSPm->getSubTripleIncrement("GGG", 22) << endl;

  // Print interpolated triplet Profile to file
  ofstream sensitivityProfileFile;
  sensitivityProfileFile.open("TripletProfile.dat");
  for (size_t baseIndex = 0; baseIndex < 64; ++baseIndex) {    
    for(size_t pos = 0; pos < 23; ++pos) {
      string baseString = sequenceutil::indexToSequence(baseIndex, 3);
      sensitivityProfileFile << pos + 1 << "\t" << baseString << "\t" << 
        simpleProfileNSPm->getSubTripleIncrement(baseString, pos) << endl;
    }
  }
  sensitivityProfileFile.close();


  // @todo Teste Qualitaet der tripel interpolation indem ich interpoliertes Profil mit
  //       "echtem" 3er Profil vergleiche

  // Next step:
  // write method to, for all triples XXX, get average contribution of triple with mismatch

  // Test: teste positionsabhaengige triplets ( schaue profil an) bevor es position herausgerechnet wird
  // @todo: verwende nur allele probesets!!!!
  IntensityArray tripletProfile = computeTripleContributions(probesets, intensitiesAllele, simpleProfileNSPm);  
  sensitivityProfileFile.open("SpecialTripletProfile.dat");
  // for (size_t baseIndex = 0; baseIndex < 64; ++baseIndex) {    
  //   for(size_t pos = 0; pos < 23; ++pos) {
  //     string baseString = sequenceutil::indexToSequence(baseIndex, 3);
  //     sensitivityProfileFile << pos + 1 << "\t" << baseString << "\t" << 
  //       tripletProfile[baseIndex * 23 + pos] << endl;
  //   }
  // }
  for (size_t baseIndex = 0; baseIndex < 64; ++baseIndex) {    
    string baseString = sequenceutil::indexToSequence(baseIndex, 3);
    sensitivityProfileFile << "13\t" << baseString << "\t" << 
      tripletProfile[baseIndex] << endl;
  }
  sensitivityProfileFile.close();

  }


  // ------------------------------------------------------------------------------------------
  // Now do the same as above, only for Cross-Allele
  // ------------------------------------------------------------------------------------------

  if (1) {
    // @test: Use only probesets with first allele present
    SnpProbesetVector probesets;
    for (size_t probesetIndex = 0; probesetIndex < probesetsX.size(); ++probesetIndex) {
      const SnpProbeset& currentProbeset = probesetsX[probesetIndex];
      //      if (currentProbeset.getCall() == 'A') {
      if (currentProbeset.getCall() != 'H') {
        probesets.push_back(currentProbeset);
        probesets.back().swapCross(); // Swap cross and Allele!!
      }
    }
 
    ProbesetIntensitiesArrayPtr intensitiesAllele = SnpProbeset::initializeIntensitiesArray(probesets);

    // Loop over all probesets 
    for (size_t probesetIndex = 0; probesetIndex < probesets.size(); ++probesetIndex) {
      const SnpProbeset& currentProbeset = probesets[probesetIndex];
      for (size_t probeIndex = 0; probeIndex < currentProbeset.getSize(); ++probeIndex) {            
        SnpProbe& currentProbe = *currentProbeset.getProbePtr(probeIndex);
        (*intensitiesAllele)[probesetIndex][probeIndex] = log10(currentProbe.getPm());
      }
    }

    // Create list of probe sequences
    ProbesetSequencesArray probesetSequences;
    for (SnpProbesetVector::const_iterator pi = probesets.begin(); pi != probesets.end(); ++pi) {
      probesetSequences.push_back(pi->getProbeSequences());
    } 


    SimpleSensitivityProfilePtr simpleProfileNSPm 
      = boost::shared_static_cast<SimpleSensitivityProfile>(sensitivityProfileNSPm);
    // Print interpolated triplet Profile to file
    ofstream sensitivityProfileFile;


    IntensityArray tripletProfile = computeTripleContributions(probesets, intensitiesAllele, simpleProfileNSPm);  
    sensitivityProfileFile.open("SpecialTripletProfileCross.dat");
    for (size_t baseIndex = 0; baseIndex < 64; ++baseIndex) {    
      string baseString = sequenceutil::indexToSequence(baseIndex, 3);
      sensitivityProfileFile << "13\t" << baseString << "\t" << 
        tripletProfile[baseIndex] << endl;
    }
    sensitivityProfileFile.close();


  }






  
  
  return(0);


  /* ---------------------------------------------------------------------------------------------

     It follows the source code of the initial version

     ---------------------------------------------------------------------------------------------
   */

  // // Get all ProbesetIds
  // vector<string> probesetIds;
  // transform(probes.begin(), probes.end(), back_inserter(probesetIds), boost::mem_fn(&PmMmProbe::getProbesetId));
  // vector<string>::iterator it = unique (probesetIds.begin(), probesetIds.end()); 
  // probesetIds.resize( it - probesetIds.begin() ); 
    
  // // Create map of each probeset id to allele A
  // map<string,char>  initialGenotypes;
  // BOOST_FOREACH(string p, probesetIds) {
  //   initialGenotypes[p] = 'A';
  // }
  
  // // Create both probeset sets
  // ProbesetVector probesetsHtA =  ComposeProbesetsByGenotypeCall(initialGenotypes, kAllele)(probes);
  // ProbesetVector probesetsHtB =  ComposeProbesetsByGenotypeCall(initialGenotypes, kCross)(probes);

  // vector<HaplotypeData> ht(2);
  // HaplotypeData htA, htB;
  // ht[kAllele].chip = ChipPtr(new Chip(probes, probesetsHtA);
  // ht[kCross].chip  = ChipPtr(new Chip(probes, probesetsHtB);
  // BOOST_FOREACH(HaplotypeData& h, ht) {
  //   h.model = HookModelPtr(new HookModel());
  //   h.analyzer = HookcurveAnalyzerPtr(new HookcurveAnalyzer(h.chip, h.model));
  // }
  
  // // Compute intensity of those sets

  // // swap elements in lists so that a.pm >= b.pm
  
  // // show hook for both




  // Create two hookcurves based on genotype call
  vector<HaplotypeData> haplotypes(2);
  vector<Haplotype> relevantHaplotypes;
  relevantHaplotypes.push_back(kAllele);
  relevantHaplotypes.push_back(kCross);
  // relevantHaplotypes.push_back(kHetero);

  vector<Haplotype>::const_iterator hapType = relevantHaplotypes.begin(); // kAllele

  // for (vector<Haplotype>::const_iterator hapType = relevantHaplotypes.begin(); hapType != relevantHaplotypes.end(); 
  //      ++hapType) {
  
  // Read initial calls from CLRMM calling
  // @bug Do not add .CEL but be robust in the Reader instead
  // map<string,char>  initialGenotypes = ComposeProbesetsByGenotypeCall::readGenotypeCalls
  //   (options.genotypeCallFilename, DataCollector::instance().getAs<string>("chipId") + ".CEL");
    //(options.genotypeCallFilename, options.intensityFilename);
    //(options.genotypeCallFilename, DataCollector::instance().getAs<string>("chipId"));

  // for (map<string,char>::const_iterator ps = initialGenotypes.begin(); ps != initialGenotypes.end(); ++ps) {
  //   cout << ps->first << "\t" << ps->second << endl;
  // }
  

  // Ab sofort gibt es ein Problem , und ich muss sehen wie ich dies strukturieren kann:
  // der Controller hat Referenzen auf Chip, Model, Hookcurve. Sprich, bei allen Controller
  // funktion koennen diese verwendet werden, daher muss man tierisch aufpassen

  // temporaere loesung: mache controllervariablen schreibbar (evil!!!!!!)
  // langfristig: nehme mehrere Controller statt der Variablen in HaplotypeData ?


  // HaplotypeData& ht = haplotypes[*hapType];
  // ht.chip 
  //   = ChipPtr(new Chip(probes, ComposeProbesetsByGenotypeCall(initialGenotypes, *hapType)(probes)));
  
  // cout << "There are "  << ht.chip->getProbes().size() 
  //      << " probes and " << ht.chip->getProbesets().size() << " probesets."  << endl;

  // controller.mChip = ht.chip; // !!!!!!
  // controller.initializeIntensityArrays(); // !!!!!!!!

  // //HookModelPtr theModel = controller.createHookModel(); // sinnlos!
  // ht.model = controller.createHookModel(); // sinnlos!
  // //ht.model = HookModelPtr(new HookModel());
  // //ht.analyzer = HookcurveAnalyzerPtr(new HookcurveAnalyzer(*ht.chip, *ht.model));

  // ht.analyzer = controller.createHookcurveModel();

  //   // ht.hookcurvePlot = controller.computeHookcurve(kRawIntensitiesOnly, ht.analyzer);
  //   // cout << "hookcurvePlot.size(): " << ht.hookcurvePlot->size() << endl;

  //   // //string filename2 = string("Hookplot-Primary") + boost::lexical_cast<std::string>(*hapType) + ".dat2";
  //   // string filename2 = string("Hookplot-Primary") + ".dat";
  //   // printVectorPairsToFile(*ht.hookcurvePlot, filename2); 
  // // }

  // IntensityType nsThreshold = 20; // unreachable high threshold
  // controller.computeSequenceProfiles(nsThreshold, string("Primary"),
  //                                    options.probesetLimitForProfileComputation);

  // cout << "5" << endl;


  // IntensityPair hookstartCoords 
  //   = controller.computeKinkCoordinates("", haplotypes[kCross].hookcurvePlot, haplotypes[kCross].analyzer);
  // cout << "hookstartCoords.first: " << hookstartCoords.first << endl;
  // cout << "hookstartCoords.second: " << hookstartCoords.second << endl;

  // hookstartCoords.first = (*haplotypes[kCross].hookcurvePlot)[0].first;
  // hookstartCoords.second = (*haplotypes[kCross].hookcurvePlot)[0].second;
  // cout << "hookstartCoords.first: " << hookstartCoords.first << endl;
  // cout << "hookstartCoords.second: " << hookstartCoords.second << endl;

  // // Subtract max(NS) and log(b) (kAllele) from both plots!
  // transform(haplotypes[kAllele].hookcurvePlot->begin(), haplotypes[kAllele].hookcurvePlot->end(), 
  //           haplotypes[kAllele].hookcurvePlot->begin(),
  //           AddToPair<IntensityPair>(IntensityPair(-hookstartCoords.first, 
  //                                                  - hookstartCoords.second)));
  // printVectorPairsToFile(*haplotypes[kAllele].hookcurvePlot, "HookplotAllele-Transformed.dat");
  
  // transform(haplotypes[kCross].hookcurvePlot->begin(), haplotypes[kCross].hookcurvePlot->end(), 
  //           haplotypes[kCross].hookcurvePlot->begin(),
  //           AddToPair<IntensityPair>(IntensityPair(- hookstartCoords.first, 
  //                                                  - hookstartCoords.second)));
  // printVectorPairsToFile(*haplotypes[kCross].hookcurvePlot, "HookplotCross-Transformed.dat");

  
  // cout << "Calculating theoretic curve." << endl;
  // IntensityType a1, f1, a2,offset;
  // tie(a1,f1, a2, offset) = 
  //   HookcurveAnalyzer::estimateTwoTheoreticFunctionParameters((haplotypes[kAllele].hookcurvePlot),
  //                                                             (haplotypes[kCross].hookcurvePlot));
  // cout << "a1: " << a1 << endl; cout << "f1: " << f1 << endl;
  // cout << "a2: " << a2 << endl; 

  // // Dump plots to a file
  // plotTheoreticCurve(a1, f1, hookstartCoords, offset, "theoCurve1.dat");
  // plotTheoreticCurve(a2, f1, hookstartCoords, offset, "theoCurve2.dat");

  // That's the last thing we do (maybe we can put this in the destructor of DataCollector ?)

 ExitProgram:
  DataCollector::instance().writeToFile("dataLogger.log");
  return 0;
}


