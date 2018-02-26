#ifndef BIG_IS_FINITE_H
#define BIG_IS_FINITE_H

#ifdef WIN32
#include <float.h>
#endif

namespace big
{
  
#ifdef WIN32
#define is_finite(val) (_finite(val)!=0)  ///< Is finite number test
#else
#define is_finite(val) (finite(val)!=0)  ///< Is finite number test
#endif
  
} // namespace big

//BIG_IS_FINITE_H
#endif
