/**
 * @file StringUtil.cpp Provides utility functions for string handling.
 * @author $Author$
 * @author Mario Fasold
 * @date $Date$
 *
 */
#include <string>
#include <vector>
#include <cctype> // for toupper
#include <algorithm>

#include "StringUtil.hpp"

using namespace std;

/**
 * Splits a string up by some delimiters.
 *
 * @param str The string
 * @param delimiters A String with characters used as delimiters
 * @return An array of strings.
 */
vector<string> stringutil::splitString(const string& str, const string& delimiters)
{
  vector<string> tokens;

  // Skip delimiters at beginning
  string::size_type tokenStartPos = str.find_first_not_of(delimiters, 0);
  
  // Find first non-delimiter
  string::size_type tokenEndPos = str.find_first_of(delimiters, tokenStartPos);

  while (string::npos != tokenStartPos) {
    // i do not need:  || string::npos != tokenEndPos

    // Found a token, add it to the vector.
    tokens.push_back(str.substr(tokenStartPos, tokenEndPos - tokenStartPos));

    // Skip delimiters
    tokenStartPos = str.find_first_not_of(delimiters, tokenEndPos);

    // Find next non-delimiter
    tokenEndPos = str.find_first_of(delimiters, tokenStartPos);
  }
  return tokens;
}

/**
 * Converts a string into its uppercase variant
 *
 * @param str Input string.
 *
 * @return String in upper case.
 */
std::string stringutil::getUppercase(const std::string str)
{
  string upcase(str); // Initialize with same size as str

  // Use toupper to convert all chars
  transform(str.begin(), str.end(), upcase.begin(), (int(*)(int)) toupper);
  return upcase;
}


/**
 * Tests if a string ends with a particular substring
 *
 * @param str Input string.
 * @param suffix Suffix to test for.
 *
 * @return True iff string ends with that particular substring.
 */
bool stringutil::endsWith(const std::string& str, const std::string suffix)
{
  // cout << "str.find_last_of(suffix): " << str.find_last_of(suffix) << endl;
  return str.rfind(suffix) == (str.length() - suffix.length());
}



/**
 * Counts the number of occurrences of a substring in a string
 * 
 * @param word query string
 * @param sentence search base
 * 
 * @return Number od occurrences of word in sentence
 */
size_t stringutil::countSubstr(const string& word, const string& sentence)
{
       size_t count = 0;
       string::size_type wordPos( 0 );
       while ( wordPos!=string::npos ){
    	   wordPos = sentence.find(word, wordPos );
               if ( wordPos != string::npos ){
            	   ++count;
                    // To find the other occurrences start search after the last known occurence
                    wordPos += word.length();
               }
       }
       return count;
}

