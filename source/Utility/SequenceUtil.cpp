/**
 * @file SequenceUtil.cpp Defines functions to work with biological sequences.
 * @author $Author$ 
 * @author Mario Fasold
 * @date $Date$
 */

#include <string>
#include <algorithm>

#include "SequenceUtil.hpp"

using namespace std;

/**
 * Complements the base at position \f[ \ceil {n/2} \f].
 *
 * @param sequence Nucleotide sequence.
 *
 * @return Sequence with complemented middle base
 *
 * @note Due to parameter pass-by-value there is the risk
 * of perfomance issues in some compilers.
 */
std::string sequenceutil::complementMiddlebase(const std::string sequence)
{
  std::string tmpSequence(sequence);
  size_t middlebaseIndex = getMiddlebaseIndex(sequence.size());
  tmpSequence[middlebaseIndex] = complement(tmpSequence[middlebaseIndex]);
  return tmpSequence;
}

/**
 * Returns complementary base, that is A if base is T, T if base is A, 
 * C if base it G and G id base is C.
 *
 * @param base char A,T,G or C
 */
int sequenceutil::complement(int base)
{
  switch (base) {
  case 'A':
    return 'T';
  case 'T':
    return 'A';    
  case 'C':
    return 'G';
  case 'G':
    return 'C';
  }
  return -1;
}

/**
 * Returns complementary base sequence
 *
 * @param sequence String of chars A,T,G or C
 */
std::string sequenceutil::complement(const std::string& sequence)
{
  std::string tmpSequence(sequence);
  transform(tmpSequence.begin(), tmpSequence.end(), tmpSequence.begin(), (int(*)(int))complement);
  return tmpSequence;
}


/**
 * Returns reverse sequence
 *
 * @param sequence String 
 */
std::string sequenceutil::reverse(const std::string& sequence)
{
  std::string tmpSequence(sequence);
  reverse(tmpSequence.begin(), tmpSequence.end());
  return tmpSequence;
}

/**
 * Returns a position specific wight matrix (PSWM) obtained from the relative
 * frequency of the bases (or base tuples) within all sequences of the set.
 *
 * @pre All sequences of the set must have the same length.
 *
 * @param modelRank Size of the base tuples within the current model
 * (mononucleotides = 1, dinucleotides = 2, ...)
 *
 * @return A position specific wight matrix.
 */  
math::Matrix<double> sequenceutil::getPositionSpecificWeightMatrix(const StringVector& sequences, 
                                                                   const size_t modelRank)
{
  // Create Matrix with zeros
  const size_t sequenceLength = sequenceutil::getValidSequenceCount(sequences.front().size(), modelRank);
  math::Matrix<double> pswm(sequenceutil::getBasetupleCount(modelRank), sequenceLength, 0);
  
  // For all sequences and all positions p, add 1 to the matrix_n_p if nucleotide tuple n
  // occurs at position p
  for (StringVector::const_iterator pi = sequences.begin(); pi != sequences.end(); ++pi) {  
    for (size_t sequenceIndex = 0; sequenceIndex < sequenceLength; ++sequenceIndex) {
      pswm[sequenceToIndex(pi->substr(sequenceIndex, modelRank))][sequenceIndex] += 1;
    }
  }

  // Devide each matrix element by the number of probes
  // to obtain relative frequencies
  pswm.flat() /= sequences.size();

  return pswm;
}

