/**
 * @file ProgramOptions.cpp Defines method to parse commandline-options
 * @author $Author$ 
 * @author Mario Fasold
 * @date $Date$
 */
#include "boost/tuple/tuple.hpp"
#include <boost/lexical_cast.hpp>
#include <time.h>
#include <algorithm>

#include "anyoption.h"

#include "ProgramOptions.hpp"
#include "StringUtil.hpp"
#include "MathUtil.hpp"
#include "StlUtil.hpp"
#include "AveragingProcedure.hpp"

using namespace std;
using namespace larrpack;
using namespace boost;
using boost::lexical_cast;


/**
 * Given a type T (e.g. T = int), the functions returns either key casted to type T (e.g. atoi(key))
 * or, if key == NULL, a given default value
 *
 */
template<class T>
T defaultAs(const char* key, const T default_value) {
  if (key == NULL) {
    return default_value;
  }
  else {
    return lexical_cast<T>(key);
  }  
}


/**
 * Parses the command line to extract program options. Also displays a help
 * message if requested and exits if command line erroneous.
 * Defaults are given for unspecified parameters.
 *
 * @param argc Number of command line parameters.
 * @param argv List of command line parameters.
 *
 * @return Program options.
 */
ProgramOptions larrpack::parseCommandlineOptions(int argc, char *argv[]) 
{
  ProgramOptions options;

  // Parse command line options
  AnyOption *opt = new AnyOption();


//   po::options_description debuggingOptions("Debugging Options");
//   debuggingOptions.add_options()
//     ("debug-mode,g", "Print additional information and statistics. Slower.")
//     ("probeset-size-filter", po::value<size_t>()->default_value(0),
//      "When set to >0, only probesets of that particular size are processed")
//     ("universal-param,y", po::value<string>()->default_value(""),
//      "Parameter that might be used for different debug options")
//     ("maskfile,m", po::value<string>()->default_value(""),
//      "The argument defines a maskfile (ascii) that lists the probesets by \n"
//      "their id that are not considered for calculation")
//     ("read-profiles", "Read sequence profiles from files. All profile .dat \n"
//      "files must reside in current working directory.")
//     ;

  opt->addUsage("General Options:");
  opt->setFlag("help", 'h');   
  opt->addUsage(" -h  --help            Produce help message.");

  opt->setOption("config-file");        
  opt->addUsage(" --config-file        Import options from given config file.");
  //   // Maybe import config file
  //   if (vm.count("config-file")) {
  //     cout << "Input file is: "
  //          << vm["config-file"].as<string>() << "\n";

  //     ifstream ifs(vm["config-file"].as<string>().c_str());
  //     store(parse_config_file(ifs, config_file_options), vm);
  // 	/* read options from a  option/resource file with ':' separated opttions or flags, one per line */
  //   opt->processFile( "/home/user/.options" );  

  //   } 

  // log-parameters ?

  opt->setOption("intensity-file", 'i');
  opt->addUsage(" -i  --intensity-file  Name of the file to contain the intensity information (CEL-file)");
  opt->setOption("chip-type", 't');
  opt->addUsage(" -t  --chip-typ        Type of the microarray chip, either genechip or tiling");
  opt->setOption("chip-file", 'c');
  opt->addUsage(" -c  --chip-file       Name of the file to contain chip specific information. That is");
  opt->addUsage("                       either a tabular, PMAP or BPMAP file, depending on chip type.");
  opt->setOption("genotype-file");
  opt->addUsage(" --genotype-file       SNP-chips only: file to contain precomputed genotype calls.");
  opt->setOption("background-subtraction-zones");
  opt->addUsage(" --background-subtraction-zones");
  opt->addUsage("                       Background subtraction devides the chip in zones^2 areas");
  opt->setOption("background-subtraction-ratio");
  opt->addUsage(" --background-subtraction-ratio         ");
  opt->addUsage("                       Multiplicative factor applied to the background value");
  opt->addUsage("                       to be subtracted from each intensity");
  opt->setOption("averaging-method",'a');
  opt->addUsage(" -a  --averaging-method");
  opt->addUsage("                       Method to robustely determine probeset average intensity, either");
  opt->addUsage("                       mean,  median or tukey (=Onestep Tukey Biweight)");
  opt->setOption("tukey-tuning");
  opt->addUsage(" --tukey-tuning        Tuning constant for the onestep Tukey Biweight Method");
  opt->setOption("gLog-parameter", 'l');
  opt->addUsage(" -l  --gLog-parameter  gLog parameter c.");
  opt->setOption("saturation-computation");
  opt->addUsage(" --saturation-computation                    ");
  opt->addUsage("                       Method to compute the saturation term I_max. Fit either the");
  opt->addUsage("                       theoretical hookcurve (0) or use the average of n strongest");
  opt->addUsage("                       probes (1)");
  opt->setOption("probeset-limit");
  opt->addUsage(" --probeset-limit      Maximum number of probes used to compute the sequence");
  opt->addUsage("                       profiles (0 = use all probes) (default: 5000 for expression");
  opt->addUsage("                       arrays, 15000 for tiling arrays)");
  opt->setFlag("print-probe-statistics");
  opt->addUsage(" --print-probe-statistics");
  opt->addUsage("                       Writes statistics and  expressionmeasures similar to");
  opt->addUsage("                       probesetStatistics.dat on probe-level. Requires for ");
  opt->addUsage("                       plenty of disk-space!");            


  opt->setOption("expression-measures", 'e');
  opt->addUsage(" --expression-measures Types of expression measures that to be calculated:");
  //     ("expression-measures,e", po::value<string>()->default_value("1100"),
  //      "What expression measures shall be calculated:\n"
  //      " 1000  : only PM\n"
  //      " 0100  : only MM\n"
  //      " 0010  : Diff (fast)\n"
  //      " 0001  : Diff (gcrma-like)\n"
  //      " 0101  (guess what this does)")
  //     // Shuffling Options
  //     ("shuffling-optimal-probes-per-set", po::value<size_t>()->default_value(5),
  //      "Number of probes that a probeset optimally contains")
  //     ("shuffling-maximal-probes-per-set", po::value<size_t>()->default_value(8),
  //      "Number of probes that a probeset maximally contains")
  //     ("shuffling-max-probe-distance", po::value<int>()->default_value(30),
  //      "Maximum allowed distance between adjacent probes within a probeset")
  //     ("shuffling-window-size", po::value<size_t>()->default_value(5),
  //      "Size of moving average window during shuffling")
  //     ("shuffling-specific-treshold", po::value<double>()->default_value(0.7),
  //      "Probability above which a probe is considered target-present")
  //     ;

  opt->setFlag("update-celfiles");
  opt->addUsage(" --update-celfiles     Writes the updated signal intensities into a CEL-file");
  opt->addUsage("                       that can be used for downstream analysis.");

  opt->addUsage("Hookcurve Options:");
  opt->setOption("detect-ip",'p');
  opt->addUsage(" -p --detect-ip        Method to calculate NS/S boundary (default=1):");
  opt->addUsage("                         0 = line and parabola");
  opt->addUsage("                         1 = two straight lines in convenient intervals");
  opt->addUsage("                         2 = first deviation");
  opt->addUsage("                         3 = maximum of second deviation");
  opt->addUsage("                         4 = first deviation, intersection with y=0");

  opt->setOption("interval-count");
  opt->addUsage(" --interval-count      Specifies the number of intervals for the FitTwoStraightLines");
  opt->addUsage("                       algorithm to detect the kink point. This should usually");
  opt->addUsage("                       not be changed. In some cases it  might be usueful to change");
  opt->addUsage("                       it to some higher value");

  opt->setOption("hookcurve-window-size",'w');
  opt->addUsage(" -w --hookcurve-window-size");
  opt->addUsage("                       For the hookcurve, the number of probesets to average (default=91)");
  
  opt->setOption("sequence-profile-rank",'r');
  opt->addUsage(" -r --sequence-profile-rank");
  opt->addUsage("                       Rank of the sequence profile model, i.e. mononucleotide (1)");
  opt->addUsage("                       dinucleotide (2) or more");
  opt->addUsage("                       ");

  opt->setFlag("use-gstacks-correction");            
  opt->addUsage(" --use-gstacks-correction");            
  opt->addUsage("                       Enable correction for for so-called g-stacks effects");

  opt->setOption("digitize-intervals",'d');
  opt->addUsage(" -d --digitize-intervals");
  opt->addUsage("                       Number of intervals the hookcurve is digitized to fit the");
  opt->addUsage("                       functions to detect kink point. 0 eqals no digitization");

  opt->setOption("set-uncorrected-kinkpoint");
  opt->addUsage(" --set-uncorrected-kinkpoint");
  opt->addUsage("                       Set kink point of uncorrected hookcurve (0 to compute)");

  opt->setOption("set-corrected-kinkpoint");
  opt->addUsage(" --set-corrected-kinkpoint");
  opt->addUsage("                       Set kink point of corrected hookcurve (0 to compute)");

  opt->setOption("set-saturation-scaling");
  opt->addUsage(" --set-saturation-scaling");
  opt->addUsage("                       Set a factor s that re-scales I_max = s * I_max.");

  opt->setFlag("debug-mode",'g');            
  opt->setOption("probeset-size-filter");
  opt->setOption("universal-param",'y');
  opt->setOption("maskfile",'m');
  opt->setFlag("read-profiles");            

//   opt->setOption("");
//   opt->addUsage(" --         ");

//   opt->setOption("size", 's'); 
//   opt->setOption("name");      
//   opt->setFlag('c');            

	/* go through the command line and get the options  */
  opt->processCommandArgs(argc, argv);
    
  // Print usage if no options or help screem wanted
	if(!opt->hasOptions() || (opt->getFlag("help") || opt->getFlag('h')) ) { 
    opt->printUsage();
    delete opt;
		exit(0);
	}

  // Define Array Type
  if (opt->getValue("chip-type") != NULL || opt->getValue('t') != NULL) {
    map<string, ChipType> chipSwitch; 
    chipSwitch["GENECHIP"] = kChiptypeGenechip;
    chipSwitch["TILING"] = kChiptypeTiling;
    chipSwitch["SNP"] = kChiptypeSnp;
    chipSwitch["EXON"] = kChiptypeExon;
    chipSwitch["GENOMIC"] = kChiptypeGenomicfile;
    string paramValue = stringutil::getUppercase(opt->getValue('t'));
    if (chipSwitch.count(paramValue)) { // Check if parameter valid
      options.chipType = chipSwitch[paramValue];
    }
    else {
      cout << paramValue << " is no valid chip type. Select another one." << endl;
      exit(0);
    }
  } 
  else { // Set default type
    options.chipType = kChiptypeGenechip; 
  }


  // Get filename for CEL-file
  if(opt->getValue("intensity-file") != NULL || opt->getValue('i') != NULL) {
    options.intensityFilename = opt->getValue('i');
  }
  else {
    cout << "Providing an intensity file is mandatory." << endl;
    exit(0);
  }

  // Get filename for tabular file
  if(opt->getValue("chip-file") != NULL || opt->getValue('c') != NULL) {
    options.chipDiscriptionFilename = opt->getValue('c');
  }
  else {
    cout << "Providing a chip discription file is mandatory." << endl;
    exit(0);
  }

  // Get genotype call file for SNP chips
  if (options.chipType == kChiptypeSnp) {
    if (opt->getValue("genotype-file") != NULL) {
      options.genotypeCallFilename = opt->getValue("genotype-file");
    }
    else {
      cout << "Providing a genotype call file is mandatory." << endl;
      exit(0);
    }  
  }
  
  // Apply options for background correction (currently only one method available)
  options.backgroundCorrection.reset
    (new BackgroundSubtraction(defaultAs<size_t>(opt->getValue("background-subtraction-zones"), 4),
                               defaultAs<double>(opt->getValue("background-subtraction-ratio"), 0.5)));

  // Define which averaging method to use
  map<string, IntensitiesFunction> averageSwitch; 
  averageSwitch["MEAN"] = Mean<IntensityType>();
  averageSwitch["MEDIAN"] = Median<IntensityType>();
  averageSwitch["TUKEY"] = OnestepTukeyBiweight<IntensityType>(defaultAs<double>(opt->getValue("tukey-tuning"), 5.0), 0.0001);
  string paramValue = stringutil::getUppercase(defaultAs<string>(opt->getValue("averaging-method"), "tukey"));
  if (averageSwitch.count(paramValue)) { // Check if parameter valid
    computeAverageIntensity = averageSwitch[paramValue];
  }
  else {
    cout << defaultAs<string>(opt->getValue("averaging-method"), "tukey") << " is no valid averaging-method. ";
    cout << "Select another one." << endl;
    exit(0);
  }

  // Define procedure to find kink point of hookcurve
  int detectIp = defaultAs<int>(opt->getValue("detect-ip"), 2);
  if (detectIp == 0) {
    options.detectKinkPoint.reset(new FitStraightLineAndParabola());
  } 
  else if (detectIp == 1) {
    size_t intervalCount = defaultAs<size_t>(opt->getValue("interval-count"), 5); 
    options.detectKinkPoint.reset(new FitTwoStraightLines(intervalCount, 0.5));
  } 
  else if (detectIp >= 2 && ((defaultAs<size_t>(opt->getValue("digitize-intervals"), 100) < 30) || (defaultAs<size_t>(opt->getValue("digitize-intervals"), 100) > 100))) {
    cout << "For detect-ip of 2, 3 and 4 (given: " << detectIp << ") digitize-intervals has to be between 30 and 100 ";
    cout << "(given " << defaultAs<size_t>(opt->getValue("digitize-intervals"), 100) << ")." << endl;
    exit(0);
  }
  else  if (detectIp == 2) {
    options.detectKinkPoint.reset(new FirstDerivativeAnalysis());
  }
  else if (detectIp == 3) {
    options.detectKinkPoint.reset(new SecondDerivativeAnalysis());
  }
  else if (detectIp == 4) {
    options.detectKinkPoint.reset(new FirstDerivativeIntersectionWithZeroAnalysis());
  }
  else {
    cout << "Valid values for detect-ip are 0-4, given:" << detectIp << endl;
    exit(0);
  }

  // If set, define number of digitize intervals by applying decorator
  /// @note remove this: The new digitize method makes this obsolete
  //  if (defaultAs<size_t>(opt->getValue("digitize-intervals"), ) > 0) {
  //    DetectKinkPointPtr tmp = options.detectKinkPoint;
  //    options.detectKinkPoint.reset(new FitDigitizedPlot(tmp, defaultAs<size_t>(opt->getValue("digitize-intervals"), )));
  //  }

  // Get kink point computation mode / value
  options.uncorrectedKinkPoint = defaultAs<double>(opt->getValue("set-uncorrected-kinkpoint"), 0.0);
  options.computeUncorrectedKinkPoint = true;
  if (options.uncorrectedKinkPoint > 0) { // Don't computate kinkpoint if valid value >0 given
    options.computeUncorrectedKinkPoint = false;
  }

  options.correctedKinkPoint = defaultAs<double>(opt->getValue("set-corrected-kinkpoint"), 0.0);
  options.computeCorrectedKinkPoint = true;
  if (options.correctedKinkPoint > 0) { // Don't computate kinkpoint if value set to >0
    options.computeCorrectedKinkPoint = false;
  }

  options.ImaxScalingFactor = defaultAs<double>(opt->getValue("set-saturation-scaling"), 1.0);

  // Get rank of sensitivity profile
  options.profileModelRank = defaultAs<size_t>(opt->getValue("sequence-profile-rank"), 2);

  if (1) { // opt->getValue("expression-measures") != NULL) {
	  options.expressionMeasureIdentifier = defaultAs<string>(opt->getValue("expression-measures"), "1100");
	  for (size_t i = 0; i < options.expressionMeasureIdentifier.size(); ++i) {
		  if (!(options.expressionMeasureIdentifier[i] == '1' || options.expressionMeasureIdentifier[i] == '0')) {
			  cout << "Invalid option for expression-measures." << options.expressionMeasureIdentifier << endl;
			  exit(0);
		  }
	  }
	  if (options.expressionMeasureIdentifier.size() < 4) {
		  for (size_t var = options.expressionMeasureIdentifier.size(); var < 4; ++var) {
			  options.expressionMeasureIdentifier += "0";
		  }
	  }

    // Collect profiles that need to be computed 
    string em = options.expressionMeasureIdentifier;
    if (em[0] == '1') // calculatePm
      options.requiredProfileTypes.push_back(kEmProfilePm);
    if (em[1] == '1') // calculateMm
      options.requiredProfileTypes.push_back(kEmProfileMm);
    if (em[2] == '1') // calculateDiffFast
      options.requiredProfileTypes.push_back(kEmProfileDiffFast);
    if (em[3] == '1') // calculateDiffGCRMA
      options.requiredProfileTypes.push_back(kEmProfileGcrmaFast);
  } 

  // Set tuples for gstacks correction
  if (opt->getFlag("use-gstacks-correction")) {
	  options.gstackTuples.push_back("GGG");
    //    options.gstackTuples.push_back("GGGG");
    //    options.gstackTuples.push_back("GGGC");
    //    options.gstackTuples.push_back("GGGA");
    //    options.gstackTuples.push_back("GGGT");
    //    options.gstackTuples.push_back("CCCC");
    //    options.gstackTuples.push_back("GCCC");
    //    options.gstackTuples.push_back("CCCG");
	  
    //     options.gstackTuples.push_back("GGG");
    //     options.gstackTuples.push_back("GCC");
    //     options.gstackTuples.push_back("CCC");
  }

  // Get gookcurve moving average window size
  options.hookcurveMovingAverageWindowSize = defaultAs<size_t>(opt->getValue("hookcurve-window-size"), 91);

  // Set global gLog parameter
  MathUtil::gGLogParameter = defaultAs<double>(opt->getValue("gLog-parameter"), 4);

  // Set saturation computation method
  options.fitTheoreticHookcurve = true;
  if (defaultAs<size_t>(opt->getValue("saturation-computation"), 0) == 1) {
    options.fitTheoreticHookcurve = false;
  }

  /// Maximum probesets to use for the computation of sequence profiles
  options.probesetLimitForProfileComputation = defaultAs<size_t>(opt->getValue("probeset-limit"), 5000);
  
  // @bug We assume 5000 is default, and then set values for genechip/tiling
  // -> Setting the value to 5000 by hand (tiling) is impossible / gives false results
  if (options.probesetLimitForProfileComputation == 5000) {
    // For tiling and genechips different numbers are adequate.
    if (options.chipType == kChiptypeGenechip) {
      options.probesetLimitForProfileComputation = 5000;
    }
    else { //  if (options.chipType == kChiptypeTiling) {
      options.probesetLimitForProfileComputation = 15000;
    }
  }

  // Get shuffling options
  options.optimalProbesPerSet = defaultAs<size_t>(opt->getValue("shuffling-optimal-probes-per-set"), 5);
  options.maxProbesPerSet = defaultAs<size_t>(opt->getValue("shuffling-maximal-probes-per-set"), 8);
  options.maxProbeDistance = defaultAs<int>(opt->getValue("shuffling-max-probe-distance"), 30);
  options.shufflingMovingAverageWindowSize = defaultAs<size_t>(opt->getValue("shuffling-window-size"), 5);
  options.specificityShufflingTreshold = defaultAs<double>(opt->getValue("shuffling-specific-treshold"), 0.7);

  // Set debugging options
  options.probesetSizeFilter = defaultAs<size_t>(opt->getValue("probeset-size-filter"), 0);
  options.universalParam = defaultAs<string>(opt->getValue("universal-param"), "");
  options.maskfile = defaultAs<string>(opt->getValue("maskfile"), "");

  if (opt->getFlag("print-probe-statistics")) {
    options.flags.insert("print-probe-statistics");
  }

  if (opt->getFlag("update-celfiles")) {
    options.flags.insert("update-celfiles");
  }


  // Log parameters to file
  // if (opt->getFlag("log-parameters")) {  // @debug Log parameters by default
  ofstream logfile;
  logfile.open("parameters.log");
  // First write used program version
  // @note: revision is updated only if this files is commited!!
  logfile << "# Used program revision $Rev$ builded" << endl;
  logfile << "# $Date$ " << endl;
  opt->printOptions(logfile);
  logfile.close();


  // Debug options 
  options.isDebugMode = (bool) opt->getFlag("debug-mode");
  options.readProfilesFromFile = (bool) opt->getFlag("read-profiles");
   
  // Probes with 1/2(PM+MM) distant to the breakpoint by the cutOff
  // are considered to calculate the specific sensitivity profiles.
  options.specificityCutoff = 0.7;

  return options;
}
