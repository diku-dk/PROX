#ifndef BIG_PROD_H
#define BIG_PROD_H

#include <big_types.h>

namespace big
{
  
  /**
   * Compute y = prod(A,x)
   *
   * @param A    A compressed matrix.
   * @param x    A vector.
   * @param y    Upon return this argument holds the result of A times x.
   */
  template<typename T>
  void prod(
            boost::numeric::ublas::compressed_matrix<T> const & A
            , boost::numeric::ublas::vector<T> const & x
            , boost::numeric::ublas::vector<T>       & y
            );
  
  /**
   * Compute y = prod(A,x)*s
   *
   * @param A    A compressed matrix.
   * @param x    A vector.
   * @param s    A scaling
   * @param y    Upon return this argument holds the result of A times x times s.
   */
  template<typename T>
  void prod(
            boost::numeric::ublas::compressed_matrix<T> const & A
            , boost::numeric::ublas::vector<T> const & x
            , T const & s
            , boost::numeric::ublas::vector<T>       & y
            );
  
} // end namespace big

// BIG_PROD_H
#endif
