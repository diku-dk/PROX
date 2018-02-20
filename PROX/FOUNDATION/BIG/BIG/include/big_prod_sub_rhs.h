#ifndef BIG_PROD_SUB_RHS_H
#define BIG_PROD_SUB_RHS_H

#include <big_types.h>


namespace big
{
  /**
   * Compute y  = prod(A,x) - b
   *
   * @param A        The matrix.
   * @param x        A vector to be multiplied with the matrix.
   * @param b        The right hand side vector.
   * @param y        Upon return this argument holds the value of the computation.
   */
  template<typename T>
  void prod_sub_rhs(
                    boost::numeric::ublas::compressed_matrix<T> const & A
                    , boost::numeric::ublas::vector<T> const & x
                    , boost::numeric::ublas::vector<T> const & b
                    , boost::numeric::ublas::vector<T>       & y
                    );
  
} // end of namespace big

// BIG_PROD_SUB_RHS_H
#endif
