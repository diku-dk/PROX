#ifndef BIG_PRECISION_H
#define BIG_PRECISION_H

#include <limits>  // for std::numeric_limits<T>::epsilon()

namespace big
{
  template <typename T>
  inline T machine_precision()
  {
    return std::numeric_limits<T>::epsilon();
  }
  
  template <typename T>
  inline T working_precision()
  {
    return std::numeric_limits<T>::epsilon()*10;
  }
  
  template <typename T>
  inline T working_precision(unsigned int scale_factor)
  {
    return std::numeric_limits<T>::epsilon()*scale_factor;
  }
  
}// end of namespace big

//BIG_PRECISION_H
#endif
