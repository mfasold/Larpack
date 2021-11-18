/**
 * @file SequenceUtil.hpp Defines functions and constants to work with biological sequences.
 * @author $Author$ 
 * @author Mario Fasold
 * @date $Date$
 */
#ifndef _SEQUENCEUTIL_
#define _SEQUENCEUTIL_

#include <cassert>
#include <cmath>
#include <vector>
#include <string>

#include "Matrix.hpp"

namespace sequenceutil {
  // Global constants 
  typedef std::vector<std::string> StringVector;

  /// Number of nucleotides
  const size_t kNucleotideCount = 4;

  /// Assumed length of probes
  const size_t kProbeSequenceLength = 25;

  /**
   * Returns 0 if base = A, 1 if base = C, 2 if base = G, 3 if base = T;
   *
   * @todo This function is very, very performace critical. Use a table mapping
   * ascii code to index to further improve performace.
   *
   * @param base Character in {A,C,G,T}
   * @return Integer in 0..3.
   *
   * @todo Source out to .cpp file
   */
  inline size_t baseToIndex(char base)
  {
    switch (base) {
    case 'A': 
      return 0;
    case 'C':
      return 1;    
    case 'G':
      return 2;
    case 'T':
      return 3;
    }
    return -1;
  }

  /**
   * This function is a generalization of baseToIndex to sequences. E.g. for
   * a sequence of length 3 to it yields values from 0 ("AAA") to 63 ("TTT"). 
   *
   * @note This function is performance critical. Speeding up the code ot usage
   * of this function, e.g. by using c-strings and avoiding string constructing,
   * should give a performance boost.
   *
   * @param sequence Base sequence {A,C,G,T}^n
   * @return Integer in \f$ 0..4^n - 1 \f$
   *
   * @todo Source out to .cpp file
   */
  inline int sequenceToIndex(const std::string& sequence)
  {
    // Beginning from the last char in the sequence, we add up
    // baseToIndex(char(i)) * 4^(length - i - 1)
    int currentMultiplier = 1;
    int index = 0;
    for(size_t i = sequence.size(); i > 0; --i ) {
      index += currentMultiplier * baseToIndex(sequence[i - 1]);
      currentMultiplier <<= 2; // Multiply by 4
    }
    return index;
  }

  
  const char indexToBaseMap[] = "ACGT";

  /**
   * Returns A if index = 0, C if index = 1, G if index = 2, T if index = 3;
   *
   * @param index integer in {0,1,2,3}
   * @return Base in {A,C,G,T}
   *
   * @todo Source out to .cpp file
   */
  inline char indexToBase(int index)
  {
    assert(index < 4 && index >= 0);
    return indexToBaseMap[(size_t)index];
  }

  /**
   * Converts an index to a nuleotide sequence given a sequence length.
   * In other words, it is the extension of indexToBase to a sequence
   * of arbitrary length. 
   *
   * @param index An interger in \f$ 0..4^modelRank-1 \f$
   * @param modelRank Number of bases within the current model
   * (mononucleotides = 1, dinucleotides = 2, ...)
   *
   * @return Nucleotide sequence {A,C,G,T}^n
   *
   * @todo Source out to .cpp file
   */
  inline std::string indexToSequence(const int index, const unsigned short modelRank)
  {
    assert(index < pow(4.0, modelRank) && index >= 0);
    std::string sequence(modelRank,' '); // Init with length modelRank

    // Iterate from last base to first
    int remainingIndex = index;
    for(size_t i = modelRank; i > 0; --i) {
      // The current base is obtained using modulo 4
      sequence[i-1] = indexToBase(remainingIndex % 4);

      // Step to next index by deviding by 4
      remainingIndex = remainingIndex / 4;
    }
    return sequence;
  }


  /// Returns complementary base
  int complement(int base);

  /**
   * Returns the sequence index of the middle base.
   *
   * @param sequenceLength Length of the sequence.
   * @return Position of the middle base.
   */
  inline size_t getMiddlebaseIndex(size_t sequenceLength)
  {
    // Use integer division for both equal and odd sequence lengths
    return (sequenceLength - 1) / 2; // -1 since string index begins with zero 
  }
  
  // Reverts the middle base
  std::string complementMiddlebase(const std::string sequence);

  // Reverts the complete Sequence
  std::string complement(const std::string& sequence);

  // Reverses the base sequence
  std::string reverse(const std::string& sequence);

  /// 
  /**
   * Returns the number of possible nucleiotide-tuples with 
   * the current model rank r, which is \f$ 4^r \f$ for 4
   * nucleotides.
   *
   * @param modelRank Number of bases within the current model
   * (mononucleotide = 1, dinucleotide = 2, ...)
   */
  inline size_t getBasetupleCount(size_t modelRank)
  { 
    return (size_t) pow((float) sequenceutil::kNucleotideCount, (int) modelRank); 
  }

  /**
   * Returns number of base position taken into account with the 
   * current model rank \f$ r \f$. That is the number of sequence 
   * indices on which a complete tuple can start, \f$ n - r + 1 \f$.
   *
   * @param sequenceLength Length of sequence
   * @param modelRank Number of bases within the current model
   * (mononucleotides = 1, dinucleotides = 2, ...)
   */
  inline size_t getValidSequenceCount(size_t sequenceLength, size_t modelRank)
  { 
    return sequenceLength - modelRank + 1;
  }


  // Reverses the base sequence
  std::string reverse(const std::string& sequence);


  // Returns a position specific wight matrix (PSWM)
  math::Matrix<double> getPositionSpecificWeightMatrix(const StringVector& sequences, 
                                                       const size_t modelRank);

  /**
   * @class ContainsBasetuples
   * @brief Defines an unary operator for sequence lists such that it returns true
   * if and only if at least one sequence of the set contains one of the given
   * basetuples
   * 
   * @todo Test. Source out to cpp and document.
   */
  class ContainsBasetuples : public std::unary_function<StringVector, bool> 
  {
  public:
    explicit ContainsBasetuples(StringVector relevantBasetuples) 
      : mRelevantBasetuples(relevantBasetuples) {}

    bool operator() (const StringVector& probeset) const 
    { 
      // For all sequences in set
      for (StringVector::const_iterator pi = probeset.begin(); pi != probeset.end(); ++pi) {  
        // return 1 if probe contains one of the substrings
        for (std::vector<std::string>::const_iterator tuple = mRelevantBasetuples.begin();
             tuple != mRelevantBasetuples.end(); ++tuple) {
          if (pi->find(*tuple) != std::string::npos) {
            return 1;
          }
        }
      }
      return 0;  // return 0 otherwise
    }    
    
  private:
    /// A list of base tuples 
    std::vector<std::string> mRelevantBasetuples;
  };
 
}

#endif

