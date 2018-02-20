#ifndef BIG_MOORE_PENROSE_PSEUDOINVERSE_H
#define BIG_MOORE_PENROSE_PSEUDOINVERSE_H

#include <big_types.h>
#include <big_lu.h>
#include <big_svd.h>

#include <boost/cast.hpp>             // needed for boost::numeric_cast
#include <stdexcept>

namespace big
{
  
  /**
   * Moore-Penrose Pseudoinverse.
   *
   * @param A        The matrix to be inverted.
   * @param invA     Upon return this argument holds the pseudoinverse matrix.
   * @param invert   A functor object or function pointer that specifies the matrix
   *                 inversion method to be used.
   */
  template<class matrix_type, typename invert_function>
  inline void moore_penrose_pseudoinverse(matrix_type const & A, matrix_type& invA, invert_function const & invert)
  {
    typedef typename matrix_type::size_type                size_type;
    
    if(A.size1() <= 0 || A.size2() <= 0)
      throw std::invalid_argument("moore_penrose_pseudoinverse(): A was empty");
    
    size_type m = A.size1();  ///< number of rows
    size_type n = A.size2();  ///< number of columns
    
    invA.resize(n,m,false);
    
    if(m>n)
    {
      // over-determined case, more rows than columns => possible full column rank of A (if not pseudoinverse will not exist)
      //
      //  A x = b
      //
      //  A^T A x = A^T b
      //
      //  x = ((A^T A)^{-1} A^T) b
      //
      
      matrix_type M;
      matrix_type invM;
      
      M.resize(n,n,false);
      invM.resize(n,n,false);
      
      ublas::noalias(M) = ublas::prec_prod(ublas::trans(A),A);
      invert(M,invM);
      ublas::noalias(invA) =ublas::prec_prod(invM,ublas::trans(A)) ;
    }
    else if(m<n)
    {
      // under-determined case, less rows than columns => possible full row rank of A (if not pseudoinverse will not exist)
      //
      //  x = A^T (A A^T)^{-1} b
      //
      //
      
      matrix_type M;
      matrix_type invM;
      
      M.resize(m,m,false);
      invM.resize(m,m,false);
      
      ublas::noalias(M) = ublas::prec_prod(A,ublas::trans(A));
      invert(M,invM);
      ublas::noalias(invA) = ublas::prec_prod(ublas::trans(A),invM);
    }
    else
    {
      // we must have m==n, a square system
      invert(A,invA);
    }
  }
  
  /**
   * Pseudoinverse using LU decompostion.
   *
   * @param A     Matrix to be inverted.
   * @param invA  Upon return holds the pseudoinverse of A.
   */
  template<class matrix_type>
  inline void lu_moore_penrose_pseudoinverse(matrix_type const & A, matrix_type& invA)
  {
    moore_penrose_pseudoinverse( A, invA, &lu_invert<matrix_type> );
  }
  
  /**
   * Pseudoinverse using SVD decompostion.
   *
   * @param A     Matrix to be inverted.
   * @param invA  Upon return holds the pseudoinverse of A.
   */
  template<class matrix_type>
  inline void svd_moore_penrose_pseudoinverse(matrix_type const & A, matrix_type& invA)
  {
    moore_penrose_pseudoinverse( A, invA, &svd_invert<matrix_type> );
  }
  
} // namespace big

// BIG_MOORE_PENROSE_PSEUDOINVERSE_H
#endif
