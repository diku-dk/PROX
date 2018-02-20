#ifndef BIG_PROJECT_H
#define BIG_PROJECT_H

#include <big_types.h>

namespace big
{
  
  /**
   * Project Operator.
   */
  template < typename T >
  void project(
               boost::numeric::ublas::vector<T> const & x
               , boost::numeric::ublas::vector<T> & new_x
               );

} // namespace big

// BIG_PROJECT_H
#endif
