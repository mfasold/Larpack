/**
 * @file SnpProbe.hpp Defines data type for a snp probe centric view.
 * @author $Author: kasper $ 
 * @author Mario Fasold
 * @date $Date: 2008-04-17 16:12:26 +0200 (Thu, 17 Apr 2008) $
 */
#ifndef _SNPPROBE_
#define _SNPPROBE_

#include <string>
#include <vector>
#include <valarray>
#include <boost/shared_ptr.hpp>

#include "SequenceUtil.hpp"
#include "PmMmProbe.hpp"

namespace larrpack {  
  class SnpProbe; // Forward definition to declare the next types

  /// Pointer to a SnpProbe
  typedef boost::shared_ptr<SnpProbe> SnpProbePtr;

  /// List of probes
  typedef std::vector<SnpProbePtr> SnpProbePtrVector;

  enum Haplotype {
    kAllele,
    kCross,
    kHetero
  };

  /**
   * @class SnpProbe
   * @brief Represents a single probe which has perfect match
   *        as well as mismatch intensity values.
   * 
   * SnpProbe is used for Affymetrix single nucleotide polymorphism (SNP)-
   * analysis chips. It contains extra information for the targeted allele.
   *
   */
  // class SnpProbe : public PmMmProbe
  // {
  // public:
  //   SnpProbe(const std::string& chromosome, const char strand, const std::string& sequence, 
  //            const int location, const IntensityType pm, const IntensityType mm, 
  //            const std::pair<int,int>& pmPosition, const std::pair<int,int>& mmPosition,
  //            const std::string& probesetId, char allele);

  //   // @name Getter and setter routines
  //   // @{

  //   /// Returns the allele
  //   char getAllele() const
  //   { return mAllele; }

  //   // @}

  // private:
  //   char mAllele; 
  // };

  class SnpProbe : public Probe
  {
  public:
    SnpProbe(const std::string& chromosome, const char strand, const std::string& sequence, 
             const int location, const IntensityType pm, const std::pair<int,int>& pmPosition,
             const std::string& probesetId, char allele);

    // @name Getter and setter routines
    // @{
    /// Returns the allele
    char getAllele() const
    { return mAllele; }
    // @}
    std::string toString();


  private:
    char mAllele; 
  };

  SnpProbePtrVector extractSnpProbesFromProbelist(const ProbePtrVector& probes);
  // SnpProbePtrVector extractSnpProbesFromProbelist(const ProbePtrVector& probes);
}
#endif
