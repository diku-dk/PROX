#ifndef BIG_COMPUTE_JACOBIAN_H
#define BIG_COMPUTE_JACOBIAN_H

#include <big_types.h>

namespace big
{
  
  /**
   * Compute Jacobian.
   *
   * @param A            The coefficient matrix of the linear equation y = A x + b
   * @param bitmask      A bitmask.
   * @param J            Upon return this argument holds the value of the Jacobian matrix.
   */
  template<typename T>
  void  compute_jacobian(
                         ublas::compressed_matrix<T> const & A
                         , ublas::vector<size_t> const & bitmask
                         , ublas::compressed_matrix<T>  & J
                         );
  
} // namespace big

// BIG_COMPUTE_JACOBIAN_H
#endif
