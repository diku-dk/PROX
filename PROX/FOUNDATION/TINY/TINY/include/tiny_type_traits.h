#ifndef TINY_TYPE_TRAITS_H
#define TINY_TYPE_TRAITS_H

#include <tiny_scalar_traits.h>
//#include <tiny_sse_traits.h> // 2009-07-13 Kenny: Temporarily disabled to ease debugging

namespace tiny
{

using float_traits = ScalarTraits<float>;
using double_traits = ScalarTraits<double>;
//  typedef SSETraits<float>						   sse_float_traits;  // 2009-07-13 Kenny: Temporarily disabled to ease debugging

} // namespace tiny

//TINY_TYPE_BINDERS_H
#endif

