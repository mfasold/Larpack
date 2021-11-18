/**
 * @author $Author$ 
 * @date $Date$
 */
#include <cassert>

#include "StlUtil.hpp"

using namespace std;
namespace larrpack{

/// Tell the compiler which instantitions to use (to avoid template linker errors)
/// @Ref http://www.parashift.com/c++-faq-lite/templates.html#faq-35.13
/// @Ref http://www.parashift.com/c++-faq-lite/templates.html#faq-35.14
/// @ref http://www.comeaucomputing.com/techtalk/templates/#whylinkerror
// template std::vector<int> getVectorSlice(std::vector<int>& entireVector, size_t firstIndex, size_t lastIndex);

}
