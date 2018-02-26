#ifndef BIG_PROD_SUB_H
#define BIG_PROD_SUB_H

#include <big_types.h>


namespace big
{
  
  /**
   * Compute y -= prod(A,x)
   *
   * @param A    A compressed matrix.
   * @param x    A vector.
   * @param y    Upon return this argument holds the
   *             result subtracting the value of A times x to the
   *             current value of the argument.
   */
  template<typename T>
  void prod_sub(
                boost::numeric::ublas::compressed_matrix<T> const & A
                , boost::numeric::ublas::vector<T> const & x
                , boost::numeric::ublas::vector<T>       & y
                );
  
} // end namespace big

// BIG_PROD_SUB_H
#endif
