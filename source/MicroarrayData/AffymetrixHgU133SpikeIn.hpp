/**
 * @file AffymetrixHgU133SpikeIn.hpp Defines functions to obtain spike-in concentrations 
 *       of the HG-U 133 dataset
 * @author $Author: mario $ 
 * @author Mario Fasold
 * @author Stephan Preibisch
 * @date $Date: 2007-04-20 16:55:02 +0200 (Fri, 20 Apr 2007) $
 */


namespace AffymetrixHgU133SpikeIn {
    /**
     * All available concentrations[pM]
     */
    static double picoMolars[] = {0,0.125,0.25,0.5,1,2,4,8,16,32,64,128,256,512};

    /**
     * Static arrays mapping Geware Psetcode to Affymetrix GroupID
     */
    const std::string  pSetCodes[]   = {"203508_at","204563_at","204513_s_at",
                                        "204205_at","204959_at","207655_s_at",
                                        "204836_at","205291_at","209795_at",
                                        "207777_s_at","204912_at","205569_at",
                                        "207160_at","205692_s_at","212827_at",
                                        "209606_at","205267_at","204417_at",
                                        "205398_s_at","209734_at","209354_at",
                                        "206060_s_at","205790_at","200665_s_at",
                                        "207641_at","207540_s_at","204430_s_at",
                                        "203471_s_at","204951_at","207968_s_at",
                                        "AFFX-r2-TagA_at","AFFX-r2-TagB_at","AFFX-r2-TagC_at",
                                        "AFFX-r2-TagD_at","AFFX-r2-TagE_at","AFFX-r2-TagF_at",
                                        "AFFX-r2-TagG_at","AFFX-r2-TagH_at","AFFX-DapX-3_at",
                                        "AFFX-LysX-3_at","AFFX-PheX-3_at","AFFX-ThrX-3_at"};

    const int groupIDs[]   = {11,12,13,
                              21,22,23,
                              31,31,33,
                              41,42,43,
                              51,52,53,
                              61,62,63,
                              71,72,73,
                              81,82,83,
                              91,92,93,
                              101,102,103,
                              111,112,113,
                              121,122,123,
                              131,132,133,
                              141,142,143};

    /**
     * Static array, maps Geware ExpID to Affymetrix Latin Square ExperimentID called LSExpID
     *
     * LSExpID: concat ExperimentID and ReplicatID
     * e.g.: 11,12,13 are all experiment 1, replicat 1,2,3
     *
     * if experimentID is questionned, the number without replicat is relevant
     * e.g.: 1 is experiment 1
     */
    const std::string  experiments[] = {"12_13_02_U133A_Mer_Latin_Square_Expt1_R1",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt1_R2",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt1_R3",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt2_R1",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt2_R2",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt2_R3",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt3_R1",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt3_R2",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt3_R3",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt4_R1",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt4_R2",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt4_R3",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt5_R1",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt5_R2",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt5_R3",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt6_R1",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt6_R2",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt6_R3",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt7_R1",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt7_R2",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt7_R3",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt8_R1",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt8_R2",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt8_R3",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt9_R1",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt9_R2",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt9_R3",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt10_R1",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt10_R2",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt10_R3",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt11_R1",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt11_R2",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt11_R3",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt12_R1",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt12_R2",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt12_R3",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt13_R1",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt13_R2",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt13_R3",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt14_R1",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt14_R2",
                                        "12_13_02_U133A_Mer_Latin_Square_Expt14_R3"};

    const int LSExpIDs[] = {11,12,13,
                            21,22,23,
                            31,31,33,
                            41,42,43,
                            51,52,53,
                            61,62,63,
                            71,72,73,
                            81,82,83,
                            91,92,93,
                            101,102,103,
                            111,112,113,
                            121,122,123,
                            131,132,133,
                            141,142,143};



//   /**
//    * @class AffymetrixHgU133SpikeIn
//    * @brief Provides functions to obtain spike-in concentrations 
//    *       of the HG-U 133A dataset
//    * @note Adapted from JAVA-code written by Stephan Preibisch.
//    *
//    */
//   class 
//   {
//   public:


    /**
     * Returns Affymetrix GroupID for given ProbeSetCode
     *
     * @param pSetCode String - Geware ProbeSetCode
     * @return int Affymetrix GroupID or -1 if it is none of the Spiked-in Genes
     */
    int getGroupID(std::string pSetCode)
    {
      for (size_t i = 0; i < sizeof(groupIDs)/sizeof(int); i++)
        if (stringutil::getUppercase(pSetCode) == stringutil::getUppercase(pSetCodes[i]))
          return groupIDs[i];

      return -1;
    }

    /**
     * Returns wheather a ProbeSet contains spiked-in probes
     *
     * @param pSetCode String - Geware ProbeSetCode
     * @return boolean true if it is spiked-in gene
     */    
    bool isSpikedIn(std::string pSetCode)
    {
      if (getGroupID(pSetCode) != -1)
        return true;
      else
        return false;
    }

    /**
     * Returns Concentration for given LSExpID(Affymetrix) and given GroupID(Affymetrix)
     *
     * @param LSExpID int - Affymetrix Latin Square ExperimentID
     * @param GroupID int - Affymetrix GroupID
     * @return double - Concentration or -1 for invalid parameters
     */
    double getConcentration(int LSExpID, int GroupID)
    {
      int exp, group, conc;

      exp = LSExpID / 10;
      group = GroupID / 10;

      conc = (exp + group - 2) % 14;

      if (conc < 0 || conc >= (int)(sizeof(picoMolars)/sizeof(double))) {
        return -1;
      }
      else {
        return picoMolars[conc];
      }
    }

    /**
     * Returns Affymetrix LSExpID(Affymetrix) for given Affymetrix Experiment Name (if you know it you can)
     *
     * @param expFileName String - Affymetrix Experiment Name (e.g. 12_13_02_U133A_Mer_Latin_Square_Expt5_R2 )
     * @return int - LSExpID(Affymetrix) or -1 if it is no Latin Square Experiment or wrong file name
     */
    int getLSExpID(std::string expFileName)
    {
      for (int i = 0; i < (int)sizeof(LSExpIDs)/sizeof(int); i++)
        if (stringutil::getUppercase(expFileName) == stringutil::getUppercase(experiments[i]) ||
            stringutil::getUppercase(expFileName) == stringutil::getUppercase(experiments[i]+".CEL"))
          return LSExpIDs[i];

      return -1;
    }

    /**
     * Returns Affymetrix Experiment Name for given LSExpID(Affymetrix)
     *
     * @param LSExpID int - LSExpID(Affymetrix)
     * @return String - Affymetrix Experiment Name or null if it LatinSquareExperimentID out of range
     */
    std::string getExpID(int LSExpID)
    {
      for (int i = 0; i < (int)sizeof(LSExpIDs)/sizeof(int); i++)
        if (LSExpIDs[i] == LSExpID)
          return experiments[i];

      return NULL;
    }

    /**
     * Returns probe set code for given groupID(Affymetrix)
     *
     * @param groupID int - groupID(Affymetrix)
     * @return String - probe set code or null if groupID out of range
     */
    std::string getPSetCode(int groupID)
    {
      for (int i = 0; i < (int)sizeof(groupIDs)/sizeof(int); i++)
        if (groupIDs[i] == groupID)
          return pSetCodes[i];

      return NULL;
    }  

    /**
     * Returns Concentration for given LSExpID(Affymetrix) and given GroupID(Affymetrix)
     *
     * @param experiment Affymetrix latin square experiment filename
     * @param probeset Affymetrix probeset name
     * @return double - Concentration or -1 for invalid parameters
     */
    double getConcentration(std::string experiment,  std::string probeset)
    {
      int LSExpID = getLSExpID(experiment);

      if (isSpikedIn(probeset)) {
        int groupID = getGroupID(probeset);
        return getConcentration(LSExpID, groupID);
      }
      return -1;
    }
}
