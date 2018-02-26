#ifndef BIG_PROD_ROW_H
#define BIG_PROD_ROW_H

#include <big_types.h>

namespace big
{
  
  /**
   * Compute y_i = prod(row_i(A),x).
   *
   * @param A         A matrix
   * @param x         A vector to be multipled with the matrix.
   * @param i         The index of the i'th row value that should be computed.
   *
   * @return          The value of the product of the i'th row with x-vector.
   */
  template<typename T>
  T prod_row(
             boost::numeric::ublas::compressed_matrix<T> const & A
             , boost::numeric::ublas::vector<T> const & x
             , typename boost::numeric::ublas::vector<T>::size_type i
             );
  
} // end namespace big

// BIG_PROD_ROW_H
#endif
