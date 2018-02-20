#ifndef BIG_PROD_TRANS_H
#define BIG_PROD_TRANS_H

#include <big_types.h>

namespace big
{
  
  /**
   * Compute y = prod(trans(A),x)
   *
   * @param A       The matrix.
   * @param x       A vector to be multiplied with the transpose of the specified matrix.
   * @param y       Upon return this argument holds the result of the computations.
   */
  template<typename T>
  void prod_trans(
                  boost::numeric::ublas::compressed_matrix<T> const & A
                  , boost::numeric::ublas::vector<T> const & x
                  , boost::numeric::ublas::vector<T>       & y
                  );
  
  
  /**
   * Compute y = prod(trans(A),x) + b
   *
   * @param A       The matrix.
   * @param x       A vector to be multiplied with the transpose of the specified matrix.
   * @param b       A vector to be added to the matrix vector product
   * @param y       Upon return this argument holds the result of the computations.
   */
  template<typename T>
  void prod_trans(
                  boost::numeric::ublas::compressed_matrix<T> const & A
                  , boost::numeric::ublas::vector<T> const & x
                  , boost::numeric::ublas::vector<T> const & b
                  , boost::numeric::ublas::vector<T>       & y
                  );
  
} // end namespace big

// BIG_PROD_TRANS_H
#endif
