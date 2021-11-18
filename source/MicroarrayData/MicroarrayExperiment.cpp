#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <stdio.h>

#include <boost/lexical_cast.hpp>
using boost::lexical_cast;

#include "StringUtil.hpp"
#include "MicroarrayExperiment.hpp"
#include "ValarrayUtil.hpp"
#include "MathUtil.hpp"

#include "BPMAPFileData.h"
#include "CELFileData.h"
#include "CDFFileData.h"
#include "CELFileWriter.h"

using namespace std;
using namespace larrpack;
using namespace affxbpmap;
using namespace affxcel;
using namespace affxcdf;


/**
 * Returns all probes of the experiment.
 *
 * @return All probes.
 *
 * @note mProbes is a vector of pointers, and all elements are copied.
 */
ProbePtrVector MicroarrayExperiment::getProbes() 
{
  return mProbes;
}

/**
 *  Reads microarray experiment data from a .genomic file. 
 *
 * @param filename Name of the .genomic file
 * @param readProbesetComposition If true, reads an assignment of probes to probesets from the file.
 *
 * @todo Replace exceptions with asserts.
 */
void MicroarrayExperiment::readFromGenomicFile(const string& filename, const bool readProbesetComposition)
{
  // define line eintries of genomic file
  enum LineentriesOfGenomicfile {
    kLeChromosome = 0,
    kLePosition,
    kLeSequence,
    kLeStrand,
    kLePmMean, kLePmStddev, kLePmPixels,
    kLeMmMean, kLeMmStddev, kLeMmPixels,
    kLePmX, kLePmY, kLeMmX, kLeMmY,
    kLeProbeset
  }; 

  // Try to open file
  ifstream genomicFile;
	genomicFile.open(filename.c_str());
	if (!genomicFile) {
    // raise file not found exception
    throw invalid_argument("Could not open genomic file.");
	}

  // omit file header
	string currentLine;
  getline(genomicFile, currentLine); // use std::getLine instead ifstram::getLine
                                     // since it gives back string

	// Read file line-by-line
  vector<string> lineTokens;
	while (!genomicFile.eof()) {
    getline(genomicFile, currentLine);

    lineTokens = stringutil::splitString(currentLine, "\t");

    if (lineTokens.size() >= kLeMmMean ) {
      mProbes.push_back(ProbePtr( new PmMmProbe(lineTokens[kLeChromosome], lineTokens[kLeStrand][0], 
                                                lineTokens[kLeSequence], 
                                                lexical_cast<int>(lineTokens[kLePosition]),
                                                lexical_cast<IntensityType>(lineTokens[kLePmMean]), 
                                                lexical_cast<IntensityType>(lineTokens[kLeMmMean]),
                                                pair<int,int>(lexical_cast<int>(lineTokens[kLePmX]),
                                                              lexical_cast<int>(lineTokens[kLePmY])),
                                                pair<int,int>(lexical_cast<int>(lineTokens[kLeMmX]),
                                                              lexical_cast<int>(lineTokens[kLeMmY]))) ));
    }
  }

  genomicFile.close();
}




/**
 * Reads microarray experiment data from a probesequence.txt file. 
 *
 * @param filename Name of the .txt file
 * @param intensities 2-dimensional array of intensity values
 *
 * @todo Replace exceptions with asserts.
 */
void MicroarrayExperiment::readFromProbesequenceFile(const string& filename, const IntensityMatrix& intensities)
{
  // define line eintries of genomic file
  enum LineentriesOfGenomicfile {
    kLeProbesetName = 0,
    kLePmX, kLePmY, 
    kLeProbeInterrogationPosition,
    kLeSequence,
    kLeStrand
  }; 

  // Try to open file
  ifstream probesequenceFile;
	probesequenceFile.open(filename.c_str());
	if (!probesequenceFile) {
    // raise file not found exception
    throw invalid_argument("Could not open probesequence file.");
	}

  // Omit file header
	string currentLine;
  getline(probesequenceFile, currentLine); // use std::getLine instead ifstram::getLine
                                     // since it gives back string

	// Read file line-by-line
  vector<string> lineTokens;
	while (!probesequenceFile.eof()) {
    getline(probesequenceFile, currentLine);

    lineTokens = stringutil::splitString(currentLine, "\t");

    if (lineTokens.size() >= kLeStrand) {
      size_t pmX = lexical_cast<size_t>(lineTokens[kLePmX]);
      size_t pmY = lexical_cast<size_t>(lineTokens[kLePmY]);      
      int interrgation = lexical_cast<int>(lineTokens[kLeProbeInterrogationPosition]);
      mProbes.push_back(ProbePtr( new PmMmProbe("", // Sequence/chromosome unknown
                                                lineTokens[kLeStrand][0], 
                                                // sequenceutil::reverse(sequenceutil::complement(lineTokens[kLeSequence])), 
                                                lineTokens[kLeSequence], 
                                                interrgation, // Resembles genomic gposition 
                                                intensities[pmX][pmY],
                                                intensities[pmX][pmY + 1], // MM coordinate is y+1 !
                                                pair<int,int>(pmX, pmY),
                                                pair<int,int>(pmX, pmY + 1), // MM coordinate is y+1 !
                                                lineTokens[kLeProbesetName]) ));
    }
  }
  probesequenceFile.close();
}

map<size_t,  double> probeIdToIntensityMap;

//IntensityType averageInLogSpace(vector<IntensityType>)
//{
//	vector<IntensityType> logedVector;
//	for
//}

/**
 * Reads microarray experiment data from a exon probesequence.txt file. 
 *
 * @param filename Name of the .txt file
 * @param intensities 2-dimensional array of intensity values
 *
 * @todo Replace exceptions with asserts.
 */
void MicroarrayExperiment::readFromExonSequenceFile(const string& filename, const IntensityMatrix& intensities)
{
  // define line eintries of genomic file
  enum LineentriesOfGenomicfile {
    kLeProbeId = 0,
    kLeExonId,
    kLeProbeX, kLeProbeY, 
    kLeAssembly,
    kLeChromosome,
    kLeStart,
    kLeStop,
    kLeStrand,
    kLeProbeSequence,
    kLeStrandedness,
    kLeCategory
  }; 

  // Try to open file
  ifstream probesequenceFile;
  probesequenceFile.open(filename.c_str());
  if (!probesequenceFile) {
    // raise file not found exception
    throw invalid_argument("Could not open probesequence file.");
  }

  // Omit file header  
  string cLine;
  getline(probesequenceFile, cLine); // use std::getLine instead ifstram::getLine
                                     // since it gives back string

	// Read file line-by-line
  
  vector<string> lineTokens;
  // We first fill up a hash linking probe Ids to their intensities.
  while (!probesequenceFile.eof()) {
    getline(probesequenceFile, cLine);
    if (cLine.size() < 5) {  // Ignore empty lines
    	continue;
    }
    lineTokens = stringutil::splitString(cLine, "\t");
    size_t pmX, pmY;
    pmX = lexical_cast<size_t>(lineTokens[kLeProbeX]);
    pmY = lexical_cast<size_t>(lineTokens[kLeProbeY]);
    probeIdToIntensityMap[lexical_cast<size_t>(lineTokens[kLeProbeId])] = intensities[pmX][pmY];
  }
  
  
  probesequenceFile.close(); // Is there a more convenient way to jump back to the beginning of that file?
  probesequenceFile.open(filename.c_str());
  
  // Omit file header
  string currentLine;
  getline(probesequenceFile, currentLine); // use std::getLine instead ifstram::getLine
                                     // since it gives back string

	// Read file line-by-line
  //vector<string> lineTokens;
  size_t count = 0;
  while (!probesequenceFile.eof()) {
    getline(probesequenceFile, currentLine);
    if (currentLine.size() < 5) {  // Ignore empty lines
    	continue;
    }

    lineTokens = stringutil::splitString(currentLine, "\t");
    if (lineTokens[kLeCategory] != "main") {
    	++count;
    	continue;
    }

    size_t pmX, pmY, position;
    if (lineTokens.size() >= kLeProbeSequence) {
//       cout << lineTokens[kLeProbeX] << "  " << lineTokens[kLeProbeY] << "   " << lineTokens[kLeStart] << endl;
      if (lineTokens[kLeStart][0] == '-') // check how to handle the affy control probes... what do they mean?
        continue;
      
//      try {
      pmX = lexical_cast<size_t>(lineTokens[kLeProbeX]);
      pmY = lexical_cast<size_t>(lineTokens[kLeProbeY]);
      position = lexical_cast<size_t>(lineTokens[kLeStart]);
//      }
//      catch(boost::bad_lexical_cast &) {
//    	  cout << "Line: " << currentLine << endl;
//    	  exit(0);
//      }
      }

      
      
      vector<IntensityType> pseudoMmIntensitiesV;
      for (size_t tokenIndex = (kLeCategory + 1); tokenIndex < lineTokens.size(); ++tokenIndex) {
    	  try {
    	  pseudoMmIntensitiesV.push_back(probeIdToIntensityMap[lexical_cast<size_t>(lineTokens[tokenIndex])]);
    	  }
    	  catch(boost::bad_lexical_cast &) {
    		  cout << "L      :   "  << currentLine << endl;
    		  cout << lineTokens[tokenIndex] << endl;
    		  exit(0);
    	  }
    	  
      }
      // For exons we use "pseudoMMs" calculated of the average of the log intensity values probes in the same GC-Cluster.
      valarray<IntensityType> pseudoMmIntensitiesValarray = convertVectorToValarray(pseudoMmIntensitiesV);
      pseudoMmIntensitiesValarray = log10(pseudoMmIntensitiesValarray);
      IntensityType average = exp10(pseudoMmIntensitiesValarray.sum() / pseudoMmIntensitiesValarray.size());
      

      ProbePtr tmpProbe = ProbePtr(new PmMmProbe(lineTokens[kLeChromosome], 
      //mProbes.push_back(ProbePtr( new PmMmProbe(lineTokens[kLeChromosome], 
                                         lineTokens[kLeStrand][0], 
                                         lineTokens[kLeProbeSequence], 
                                         position, // Resembles genomic position of the start of the probe
                                         intensities[pmX][pmY],
                                         average, //3.0, // We don't have MM probes on Affymetrix Exon Arrays. ///@todo we should not use PmMmProbe here....
                                         pair<int,int>(pmX, pmY),
                                         pair<int,int>(0,0),
                                         lineTokens[kLeExonId]) );
      mProbes.push_back(tmpProbe);
    }
  cout << "Skipped " << count << " times" << endl;
  probesequenceFile.close();
}


//map<size_t,  double> gcToIntensityMap;
//map<size_t,  vector<double> > gcToIntensityVectorMap;

void printvectorMap(const vector< vector<IntensityType> >& mappp) {
	for (size_t i = 0; i < mappp.size(); ++i) {
		cout << "[ ";
		for (size_t j = 0; j < mappp[i].size(); ++j) {
			cout << mappp[i][j] << "  ";
		}
		cout << " ]     ";
	}
	
}


/**
 * Reads microarray experiment data from a exon probesequence.txt file. 
 *
 * @param filename Name of the .txt file
 * @param intensities 2-dimensional array of intensity values
 *
 * @todo Replace exceptions with asserts.
 */
void MicroarrayExperiment::readPmOnlySequenceFile(const string& filename, const IntensityMatrix& intensities)
{
  // define line eintries of genomic file
  enum LineentriesOfGenomicfile {
    kLeProbeId = 0,
    kLeProbesetId,
    kLeProbeX, kLeProbeY, 
    kLeAssembly,
    kLeChromosome,
    kLeStart,
    kLeStop,
    kLeStrand,
    kLeProbeSequence,
    kLeStrandedness,
    kLeCategory
  }; 

  vector< vector<IntensityType> > gcToIntensityVectorMap(26); // The index corresponds to the number of G and Cs.
//  vector< vector<size_t> >            probesetsOfSameGc(26); // For a given gc contant we store the vector of the indeces of the bg probesets having that gc contatnt
  
  // Try to open file
  ifstream probesequenceFile;
  probesequenceFile.open(filename.c_str());
  if (!probesequenceFile) {
    // raise file not found exception
    throw invalid_argument("Could not open probesequence file.");
  }

  // Omit file header  
  string cLine;
  getline(probesequenceFile, cLine); // use std::getLine instead ifstram::getLine
                                     // since it gives back string

	// Read file line-by-line
  
  vector<string> lineTokens;
//  string tmpProbesetId = "";
  vector<IntensityType> intensitiesOfSameGC;
  size_t gcContant;
  //size_t zaehler = 0;
  // We first fill up a hash linking a gc contant to an itensitiy value by calculating the median of all .control->bgp->antigenomic backgroung probes
  // of the same GC contant.
  while (!probesequenceFile.eof()) {
    getline(probesequenceFile, cLine);
    if (cLine.size() < 5) {  // Ignore empty lines
    	continue;
    }
    lineTokens = stringutil::splitString(cLine, "\t");
    if (lineTokens[kLeCategory] != "control->bgp->antigenomic"){
    	continue;
    }
    
//    if (lineTokens[kLeCategory] != "control->bgp->antigenomic" && lineTokens[kLeCategory] != "normgene->intron") {
//    	continue;
//    }

		gcContant = stringutil::countSubstr("G", lineTokens[kLeProbeSequence]) + stringutil::countSubstr("C", lineTokens[kLeProbeSequence]); // calculateGcContant(lineTokens[kLeProbeSequence]);

		// Because we want to calculate the median in the log space, we log10 the Intensity values.
	  gcToIntensityVectorMap[gcContant].push_back(log10(intensities[lexical_cast<size_t>(lineTokens[kLeProbeX])][lexical_cast<size_t>(lineTokens[kLeProbeY])]));
  }
  
  vector<IntensityType> gcContantToMedianIntensity(26); // The index corresponds to the numer of Gs and Cs in the probe sequence.
  for (size_t i = 0; i < 26; ++i) {
	  if (gcToIntensityVectorMap[i].size() > 0) {
		  // The pseudoMM intensities are not hold in the logspace ...so we have to exp10  them.
		  gcContantToMedianIntensity[i] = exp10(MathUtil::calculateMedian(convertVectorToValarray(gcToIntensityVectorMap[i])));
		  //cout << i << "    "  << MathUtil::calculateMedian(convertVectorToValarray(gcToIntensityVectorMap[i])) << endl;
	  }
	  
  }
  
  probesequenceFile.close(); // Is there a more convenient way to jump back to the beginning of that file?
  probesequenceFile.open(filename.c_str());
  
  // Omit file header
  string currentLine;
  getline(probesequenceFile, currentLine); // use std::getLine instead ifstram::getLine
                                     // since it gives back string

	// Read file line-by-line
  //vector<string> lineTokens;
  //size_t count = 0;
  while (!probesequenceFile.eof()) {
    getline(probesequenceFile, currentLine);
    if (currentLine.size() < 5) {  // Ignore empty lines
    	continue;
    }

    lineTokens = stringutil::splitString(currentLine, "\t");
    
    if (lineTokens[kLeCategory] == "control->affx"){ // Avoid these control probes, some habe sequence length of only 23 !
    	continue;
    }
    	
    
    size_t pmX, pmY;
    if (lineTokens.size() >= kLeProbeSequence) {

      pmX 	= lexical_cast<size_t>(lineTokens[kLeProbeX]);
      pmY 	= lexical_cast<size_t>(lineTokens[kLeProbeY]);
      
      size_t position;
      if (lineTokens[kLeStart][0] == '-') {
    	  position = MathUtil::getNaNSizeT();
      }
      else {
    	  position = lexical_cast<size_t>(lineTokens[kLeStart]);
      }
        
      
      size_t gcCont = stringutil::countSubstr("G", lineTokens[kLeProbeSequence]) + stringutil::countSubstr("C", lineTokens[kLeProbeSequence]);
      IntensityType pseudoMmIntensity = gcContantToMedianIntensity[gcCont];
      string probesetId = "";
//      if (lineTokens[kLeCategory] == "control->bgp->antigenomic"  || lineTokens[kLeCategory] == "normgene->intron") {  // We mark the backgroudn control probes by changing their probesetId
      if (lineTokens[kLeCategory] == "control->bgp->antigenomic") {  // We mark the backgroudn control probes by changing their probesetId
    	  probesetId += "BackgroundProbeBG";
      }
      if (lineTokens[kLeCategory] == "normgene->intron") {  
    	  probesetId += "BackgroundProbeIn";
      }      
      probesetId += lineTokens[kLeProbesetId];
//      if (lineTokens[kLeCategory] == "control->bgp->antigenomic")
//    	  cout << probesetId << endl;
      
      
      ProbePtr tmpProbe = ProbePtr(new PmMmProbe(lineTokens[kLeChromosome], 
      //mProbes.push_back(ProbePtr( new PmMmProbe(lineTokens[kLeChromosome], 
                                         lineTokens[kLeStrand][0], 
                                         lineTokens[kLeProbeSequence], 
                                         position, // Resembles genomic position of the start of the probe
                                         intensities[pmX][pmY],
                                         pseudoMmIntensity, 
                                         pair<int,int>(pmX, pmY),
                                         pair<int,int>(0,0), // ...there are no MM probes on PM only chips
                                         probesetId) );
      mProbes.push_back(tmpProbe);
    }
//  cout << "Skipped " << count << " times" << endl;
//  probesequenceFile.close();
} //while (!probesequenceFile.eof()) 
}



/**
 * Reads microarray experiment data from a SNP probesequence.txt file. 
 *
 * @param filename Name of the .txt file
 * @param intensities 2-dimensional array of intensity values
 *
 * @todo Replace exceptions with asserts.
 */
void MicroarrayExperiment::readFromSnpSequenceFile(const std::string& filename, const IntensityMatrix& intensities)
{
  // define line eintries of genomic file
  enum LineentriesOfGenomicfile {
    kLeProbesetName = 0,
    kLePmX, kLePmY, 
    kLeProbeInterrogationPosition,
    kLeSequence,
    kLeStrand, // f or r
    kLeProbeType,
    kLeAllele
  }; 

  // Try to open file
  ifstream probesequenceFile;
	probesequenceFile.open(filename.c_str());
	if (!probesequenceFile) {     // raise file not found exception
    throw invalid_argument("Could not open probesequence file.");
	}

  
	string currentLine;
  // @bug Check if all SNP probesequence files have no header.
  // AFAIK: There is no file header in SNP files
  // // Omit file header
  // getline(probesequenceFile, currentLine);

	// Read file line-by-line
  vector<string> lineTokens;
	while (!probesequenceFile.eof()) {
    getline(probesequenceFile, currentLine);
    lineTokens = stringutil::splitString(currentLine, "\t");

    if (lineTokens.size() >= kLeAllele) {
      size_t pmX = lexical_cast<size_t>(lineTokens[kLePmX]);
      size_t pmY = lexical_cast<size_t>(lineTokens[kLePmY]);      
      int interrgation = lexical_cast<int>(lineTokens[kLeProbeInterrogationPosition]);
      mProbes.push_back(ProbePtr( new SnpProbe("", // Sequence/chromosome unknown
                                               lineTokens[kLeStrand][0], 
                                               lineTokens[kLeSequence], 
                                               interrgation, 
                                               intensities[pmX][pmY],
                                               pair<int,int>(pmX, pmY),
                                               lineTokens[kLeProbesetName],
                                               lineTokens[kLeAllele][0])));

      // // Read all entries for one probeset into set of tokens
      // vector<vector<string> > probesetLineTokens;
      // getline(probesequenceFile, currentLine);
      // lineTokens = stringutil::splitString(currentLine, "\t");

      // while (true) {
      //   string currentSnp = lineTokens[kLeProbesetName];
      //   while (lineTokens[kLeProbesetName] == currentSnp) {
      //     probesetLineTokens.push_back(lineTokens);
      //     getline(probesequenceFile, currentLine);
      //     lineTokens = stringutil::splitString(currentLine, "\t");        
      //   }

      //   // MIND END OF FILE

      //   // Get alleles for probeset, define one as cross and other as allele
      //   char alleleA = probesetLineTokens[0][kLeAllele][0];
      //   char alleleCross = alleleA;
      //   sizt_t i = 1;
      //   while (alleleCross == alleleA) {
      //     alleleCross = probesetLineTokens[i++][kLeAllele][0];
      //   }

      //   // For each allele, search approprate cross allele, create SnpPropePair and delete both from list
      //   while (!probesetLineTokens.empty()) {
      //     size_t alleleIndex = 0;
      //     while (probesetLineTokens[alleleIndex][kLeAllele] != alleleA) {
      //       alleleIndex++;
      //     }
      //     size_t crossIndex = 0;
      //     string interrogationAllele = probesetLineTokens[alleleIndex][kLeProbeInterrogationPosition];
      //     string strandAllele = probesetLineTokens[alleleIndex][kLeStrand];
      //     while (probesetLineTokens[crossIndex][kLeAllele] != alleleCross && 
      //            probesetLineTokens[crossIndex][kLeProbeInterrogationPosition] != interrogationAllele &&
      //            probesetLineTokens[crossIndex][kLeStrand] != strandAllele)  {
      //       crossIndex++;
      //     }         
      //     // @BUG? Double-Check if indices can be found 

      //     // Add snp Pair         
      //     // mProbes.push_back(ProbePtr( new SnpProbe("", // Sequence/chromosome unknown
      //     //                                          lineTokens[kLeStrand][0], 
      //     //                                          lineTokens[kLeSequence], 
      //     //                                          interrgation, 
      //     //                                          intensities[pmX][pmY],
      //     //                                          intensities[pmX][pmY + 1], // MM coordinate is y+1 !
      //     //                                          pair<int,int>(pmX, pmY),
      //     //                                          pair<int,int>(pmX, pmY + 1), // MM coordinate is y+1 !
      //     //                                          lineTokens[kLeProbesetName],
      //     //                                          lineTokens[kLeAllele][0])));

      //     // ...............
          
        
      //     // Erase both elements (from back to front order to not change indexes)
      //     if (crossIndex > alleleIndex) {
      //       probesetLineTokens.erase(probesetLineTokens.begin() + crossIndex);
      //       probesetLineTokens.erase(probesetLineTokens.begin() + alleleIndex);
      //     }
      //     else {
      //       probesetLineTokens.erase(probesetLineTokens.begin() + alleleIndex);
      //       probesetLineTokens.erase(probesetLineTokens.begin() + crossIndex);
      //     }
      //   }
      //          // probesetLineTokens.clear(); // Double-check it's empty
      // }
    }
  }
  probesequenceFile.close();
}



// bool has50PercentGCcontant(string seq)
// {
//   size_t absNo = 0;
//   for (size_t sequenceIndex = 0; sequenceIndex < seq.size(); ++sequenceIndex) {
// //     cout << seq[sequenceIndex] << endl;
//     if ((seq[sequenceIndex] == 'G') or (seq[sequenceIndex] == 'C')) {
// //       cout << &seq[sequenceIndex] << endl;
//       absNo++;
//     }
//   }
//   if ((absNo >= 12) && (absNo >= 14)) return true;
//   return false;
//   
// }



/**
 * Reads microarray experiment data from a .bpmap file. 
 *
 * @param filename Name of the .bpmap file
 * @param intensities 2-dimensional array of intensity values
 * @param probesetIds 2-dimensional array of strings with 
 *       a probeset identifier for each probe.             
 *
 * @todo Read only those intensities, which are not masked!
 */
void MicroarrayExperiment::readFromBpmapfile(const std::string& filename, const IntensityMatrix& intensities, 
                                             const StringArrayPtr& probesetIds) 
{
  // Try to open file
  CBPMAPFileData bpmap;
  bpmap.SetFileName(filename.c_str());

  if (!bpmap.Exists()) {
    cerr << "Fatal Error: Bpmap file does not exist" << endl;
    exit(0);
  }

  bpmap.Read();      
  // cout << "#PMX\tPMY\tMatchscore\tPosition\tSequence\tTopStrand" << endl;
  // 		cout << "#Seqs = " << bpmap.GetNumberSequences() << endl;
  // 		cout << "Version = " << bpmap.GetVersion() << endl << endl;

  // Iterate all sequences in the file (=chromosomes, in most cases)
  for (int sequenceIndex = 0; sequenceIndex < bpmap.GetNumberSequences(); ++sequenceIndex) {
    CGDACSequenceItem sequence;
    bpmap.GetSequenceItem(sequenceIndex, sequence);

    // Iterate over hits in that sequence (a hit is spot on the array,
    // which could either belong to a PM probe, a MM probe or a spiked-in)
    for (int hitIndex = 0; hitIndex < sequence.GetNumberHits(); ++hitIndex) {
      GDACSequenceHitItemType hit;
      sequence.GetHitItem(hitIndex, hit, true);

      // Try to get probeset Id
      string probesetId = "";
      if (probesetIds != NULL) {
        probesetId = (*probesetIds)[hit.PMX][hit.PMY];
      }

      
//        cout << hit.PMProbe[5] << endl;
//       if (has50PercentGCcontant(hit.PMProbe))
//       {
//       
      mProbes.push_back(ProbePtr( new PmMmProbe(sequence.GetName(), hit.TopStrand, 
                                                hit.PMProbe, 
                                                //sequenceutil::reverse(sequenceutil::complement(hit.PMProbe)), 
                                                hit.Position,
                                                intensities[hit.PMX][hit.PMY],
                                                intensities[hit.MMX][hit.MMY],
                                                pair<int,int>(hit.PMX, hit.PMY),
                                                pair<int,int>(hit.MMX, hit.MMY),
                                                probesetId) ));
//       }

      // cout
      //   << hit.PMX << "\t"
      //   << hit.PMY << "\t"
      //   // << hit.MMX << "\t"
      //   // << hit.MMY << "\t"
      //   << hit.MatchScore << "\t"
      //   << hit.Position << "\t"
      //   << hit.PMProbe << "\t"
      //   // << (int) hit.ProbeLength << "\t"
      //   << (int) hit.TopStrand << endl;
    }
  }	
}

/**
 * Reads microarray experiment data from a .pmap file (which is 
 * the text-file variant of a bpmap file). 
 *
 * @param filename Name of the .bpmap file
 * @param intensities 2-dimensional array of intensity values
 * @param probesetIds 2-dimensional array of strings with 
 *       a probeset identifier for each probe.             
 *
 * @todo Replace exceptions with asserts.
 */
void MicroarrayExperiment::readFromPmapfile(const std::string& filename, const IntensityMatrix& intensities, 
                                             const StringArrayPtr& probesetIds) 
{
  // define line eintries of genomic file
  enum LineentriesOfPmapfile {
    kLeSequence = 0,
    kLeStrand,
    kLeName,
    kLeLocation,
    kLePmX, kLePmY, 
    kLeMmX, kLeMmY
  }; 

  // Try to open file
  ifstream pmapFile;
	pmapFile.open(filename.c_str());
	if (!pmapFile) {
    // raise file not found exception
    throw invalid_argument("Could not open pmap file.");
	}

  // Omit file header
	string currentLine;
  getline(pmapFile, currentLine); // use std::getLine instead ifstram::getLine
                                  // since it gives back string

	// Read file line-by-line
  vector<string> lineTokens;
	while (!pmapFile.eof()) {
    getline(pmapFile, currentLine);

    lineTokens = stringutil::splitString(currentLine, "\t");

    if (lineTokens.size() >= kLeMmY) {
      // Cast PM/MM coordinates
      size_t pmX = lexical_cast<size_t>(lineTokens[kLePmX]);
      size_t pmY = lexical_cast<size_t>(lineTokens[kLePmY]);      
      size_t mmX = lexical_cast<size_t>(lineTokens[kLeMmX]);
      size_t mmY = lexical_cast<size_t>(lineTokens[kLeMmY]);      

      // Try to get probeset Id
      string probesetId = "";
      if (probesetIds != NULL) {
        probesetId = (*probesetIds)[pmX][pmY];
      }

      // Add probe
      mProbes.push_back(ProbePtr( new PmMmProbe(lineTokens[kLeName],
                                                lineTokens[kLeStrand][0], 
                                                lineTokens[kLeSequence], 
                                                lexical_cast<size_t>(lineTokens[kLeLocation]),
                                                intensities[pmX][pmY],
                                                intensities[mmX][mmY], 
                                                pair<int,int>(pmX, pmY),
                                                pair<int,int>(mmX, mmY)
                                                ) ));
    }
  }
  pmapFile.close();

}

/**
 * Reads microarray intensities and masking from a CEL file
 *
 * @param filename File to read.
 * @param intensityArrayPtr Return value: pointer to array of intensities
 * @param isMaskedArrayPtr Return value: pointer to bool array, stating 
 *        wheather or not an array spot is masked
 *
 * @note Array is saved in [x][y] format. Change this to c++ general [y][x].
 */
void MicroarrayExperiment::readArrayFromCelfile(const std::string& filename, IntensityMatrixPtr& intensityArrayPtr, 
                                                BoolArrayPtr& isMaskedArrayPtr)
{
  CCELFileData cel;
	cel.SetFileName(filename.c_str());

  if (!cel.Exists()) {
    cerr << "Fatal Error: CEL file does not exist." << endl;
    exit(0);
  }

  cel.Read();
  //     cout << "Version = " << cel.GetVersion() << endl;
  // 		cout << "#Cols = " << cel.GetCols() << endl;
  // 		cout << "#Rows = " << cel.GetRows() << endl;
  // 		cout << "#Total = " << cel.GetNumCells() << endl;
  // 		cout << "Header = " << cel.GetHeaderString() << endl;
  // 		cout << "Alg = " << cel.GetAlg() << endl;
  // 		cout << "Params = " << cel.GetParams() << endl;
  // 		cout << "Array type = " << cel.GetChipType() << endl;
  // 		cout << "Margin = " << cel.GetCellMargin() << endl;
  // 		cout << "#Outliers = " << cel.GetNumOutliers() << endl;
  // 		cout << "#Masked = " << cel.GetNumMasked() << endl;

  intensityArrayPtr = IntensityMatrixPtr(new IntensityMatrix(cel.GetRows(),cel.GetCols()));
  isMaskedArrayPtr = BoolArrayPtr(new BoolArray(cel.GetRows(), cel.GetCols(), true));

  for (int iy = 0; iy < cel.GetRows(); iy++) {
    for (int ix = 0; ix < cel.GetCols(); ix++) {
      (*intensityArrayPtr)[ix][iy] = cel.GetIntensity(ix,iy);
      (*isMaskedArrayPtr)[ix][iy] = cel.IsMasked(ix, iy);
    }
  }
}

/**
 * Writes microarray intensities to a CEL file using another CEL file as
 * template
 *
 * @param filename File to read.
 * @param intensityArrayPtr Return value: pointer to array of intensities
 * @param isMaskedArrayPtr Return value: pointer to bool array, stating 
 *        wheather or not an array spot is masked
 *
 * @todo Find better bugfix than writing and re-reading file...
 * @note Array is saved in [x][y] format. Change this to c++ general [y][x].
 */
void MicroarrayExperiment::updateCelfileIntensities(const std::string sourceFilename, 
                                                    const std::string targetFilename, 
                                                    IntensityMatrixPtr& intensityArrayPtr)
{
  // Read original celfile
  CCELFileWriter cel;
	cel.SetFileName(sourceFilename.c_str());

  if (!cel.Exists()) {
    cerr << "Fatal Error: CEL file does not exist." << endl;
    exit(0);
  }
  cel.Read();
  
  // Bugfix of Affymetrix library. We must re-read binary format celfiles
  // in text formmat to work correctly
  if (cel.GetFileFormat() != CCELFileData::TEXT_CEL) {
    // Write in text-format
    cel.SetFileName("dummy.cel");
    cel.WriteTextCel();

    // Re-read
    cel.SetFileName("dummy.cel");
    cel.Read();    

    // remove dummy.cel
    remove("dummy.cel");
  }
  
  // Update intensities
  for (int iy = 0; iy < cel.GetRows(); iy++) {
    for (int ix = 0; ix < cel.GetCols(); ix++) {
      if (!cel.IsMasked(ix, iy)) {
        cel.SetIntensity(ix,iy, (*intensityArrayPtr)[ix][iy]);
      }
    }
  }

  // Write to new celfile
  cel.SetFileName(targetFilename.c_str());
  cel.WriteTextCel();
}


/**
 * Reads probeset identifiers for each spot on the array from CDF file
 *
 * @param filename File to read.
 * @note Array is saved in [x][y] format. Change this to c++ general [y][x].
 */
StringArrayPtr MicroarrayExperiment::readProbesetIdsFromCdffile(const std::string& filename)
{
  StringArrayPtr probesetnames;

	CCDFFileData cdf;
	cdf.SetFileName(filename.c_str());

  if (!cdf.Exists()) {
    cerr << "Fatal Error: CDF file does not exist." << endl;
    exit(0);
  }

  cdf.Read();

  probesetnames = StringArrayPtr(new StringArray(cdf.GetHeader().GetRows(), 
                                                 cdf.GetHeader().GetCols(), ""));

  CCDFFileHeader header = cdf.GetHeader();   
  CCDFProbeInformation cel;
  CCDFProbeSetInformation set;
  CCDFProbeGroupInformation group;

  for (int ips=0; ips<header.GetNumProbeSets(); ips++)
		{
			string name = cdf.GetProbeSetName(ips);
      // 			cout << endl << "Probe set #" << ips+1 << endl;
      // 			cout << "Name" << " = " << name << endl;
      // 			cout << "Type = " << cdf.GetProbeSetType(ips) << endl;
			cdf.GetProbeSetInformation(ips, set);

      // 			cout
      // 				<< endl
      // 				<< "#lists = " << set.GetNumLists() << endl
      // 				<< "#groups = " << set.GetNumGroups() << endl
      // 				<< "#cells = " << set.GetNumCells() << endl
      // 				<< "#cellsperlist = " << set.GetNumCellsPerList() << endl
      // 				<< "Number = " << set.GetProbeSetNumber() << endl
      // 				<< "Type = " << set.GetProbeSetType() << endl
      // 				<< "Dir = " << set.GetDirection() << endl;

			int n = set.GetNumGroups();
			for (int ig=0; ig<n; ig++)
        {
          set.GetGroupInformation(ig, group);

          // 				cout
          // 					<< endl
          // 					<< "#lists = " << group.GetNumLists() << endl
          // 					<< "#cells = " << group.GetNumCells() << endl
          // 					<< "#cellsperlist = " << group.GetNumCellsPerList() << endl
          // 					<< "start = " << group.GetStart() << endl
          // 					<< "stop = " << group.GetStop() << endl
          // 					<< "Name = " << group.GetName() << endl;

          for (int ic=0; ic<group.GetNumCells(); ic++)
            {
              group.GetCell(ic, cel);
              (*probesetnames)[cel.GetX()][cel.GetY()] = name;
          
              // 					cout
              // 						<< ic << " =\t"
              // 						<< cel.GetX() << "\t"
              // 						<< cel.GetY() << "\t"
              // 						<< cel.GetListIndex() << "\t"
              // 						<< cel.GetExpos() << "\t"
              // 						<< cel.GetPBase() << "\t"
              // 						<< cel.GetTBase() << endl;
            }
        }
		}
  return probesetnames; 
}
