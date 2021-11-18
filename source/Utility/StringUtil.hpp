/**
 * @file StringUtil.hpp Provides utility functions for string handling.
 * @author $Author$
 * @author Mario Fasold
 * @date $Date$
 */
#include <string>
#include <vector>

#ifndef _STRINGUTIL_
#define _STRINGUTIL_

namespace stringutil {
  std::vector<std::string> splitString(const std::string& str, const std::string& delimiters = " ");
  std::string getUppercase(const std::string str);
  bool endsWith(const std::string& str, const std::string suffix);
  size_t countSubstr(const std::string& word, const std::string& sentence);
}
#endif
