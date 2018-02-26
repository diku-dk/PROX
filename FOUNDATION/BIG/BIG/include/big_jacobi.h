#ifndef BIG_JACOBI_H
#define BIG_JACOBI_H

#include <big_types.h>
#include <big_prod_row.h>
#include <big_is_number.h>
#include <big_value_traits.h>

#include <boost/cast.hpp> // needed for boost::numeric_cast

#include <cmath>        //needed for std::fabs
#include <stdexcept>
#include <cassert>

namespace big
{
  
  /**
   * Jacobi Iteration Method.
   * This function computes the result of a single iteration using the Jacobi method.
   *
   * @param A       A matrix, A = L + D + U, where D is diagonal, L strict lower triangular, and U strict upper triangular.
   * @param x       Upon return this argument holds the value of x^{k+1} = D^{-1}(b - (L+U) x^k).
   * @param b       The right hand side vector.
   *
   */
  template<typename T>
  inline void jacobi(
                     ublas::compressed_matrix<T>  const & A
                     , ublas::vector<T>                 & x
                     , ublas::vector<T>           const & b
                     )
  {
    using std::fabs;
    
    typedef ublas::vector<T>                vector_type;
    typedef typename vector_type::size_type size_type;
    typedef big::ValueTraits<T>             value_traits;
    
    if(A.size1() <= 0u || A.size2() <= 0u)
      throw std::invalid_argument("jacobi(): A was empty");
    if(b.size() != A.size1())
      throw std::invalid_argument("jacobi(): The size of b must be the same as the number of rows in A");
    if(x.size() != A.size2())
      throw std::invalid_argument("jacobi(): The size of x must be the same as the number of columns in A");
    
    size_type const n = x.size();
    vector_type x_old = x;
    for ( size_type i = 0u; i < n; ++i )
    {
      // This is the straigthforward way of implementing a Jacobi iteration
      //x( i ) = b( i );
      //for ( size_type j = 0; j < i; ++j )
      //{
      //  x( i ) -= A( i, j ) * x_old( j );
      //}
      //for ( size_type j = i + 1; j < n; ++j )
      //{
      //  x( i ) -= A( i, j ) * x_old( j );
      //}
      //if ( A( i, i ) )
      //  x( i ) /= A( i, i );
      //
      // However, knowing we work on a compressed matrix
      // we can do a lot better
      
      T const & Aii  = A(i,i);

      assert( is_number( x(i) )                || !"jacobi(): sum value was not a number?");
      assert( is_number( Aii )                 || !"jacobi(): diag value was not a number?");
      assert( fabs(Aii) > value_traits::zero() || !"jacobi(): Diagonal were zero");
      
      x(i) = ( b(i) - prod_row(A,x_old,i)) / Aii + x_old(i);

      assert( is_number( x(i) )              || !"backward_gauss_seidel(): updated value was not a number?");
    }
  }
  
  /**
   * Jacobi Method.
   * This function is capable of running several iterations of the Jacobi method.
   *
   * @param A                A matrix, A = L + D + U, where D is diagonal, L strict lower triangular, and U strict upper triangular.
   * @param x                At invokation this argument holds the initial value of x^0, Upon return this argument holds the value of x^{k+1} = D^{-1}(b - (L+U) x^k) for some k>0.
   * @param b                The right hand side vector.
   * @param max_iterations   The maximum number of iterations that is allowed.
   * @param epsilon          A stopping threshold (currently not used).
   * @param iterations       Upon return this argument holds the number of used iterations.
   */
  template<typename T>
  inline void jacobi(
                     ublas::compressed_matrix<T>  const & A
                     , ublas::vector<T>                 & x
                     , ublas::vector<T>           const & b
                     , size_t                     const & max_iterations
                     , T                          const & /*epsilon*/
                     , size_t                           & iterations
                     )
  {
    if(max_iterations < 1u)
      throw std::invalid_argument("jacobi(): max_iterations must be a positive number");
    
    iterations = 0u;
    while ( iterations < max_iterations )
    {
      jacobi(A,x,b);
      ++iterations;
    }
  }
  
  /**
   * Jacobi Method.
   * This function is capable of running several iterations of the Jacobi method.
   * It uses a default setting of maximum 15 iterations and a stopping
   * threshold of 10e-4.
   *
   * This function is merely a convenience function with typical default settings.
   *
   * @param A                A matrix, A = L + D + U, where D is diagonal, L strict lower triangular, and U strict upper triangular.
   * @param x                At invokation this argument holds the initial value of x^0, Upon return this argument holds the value of x^{k+1} = D^{-1}(b - (L+U) x^k) for some k>0.
   * @param b                The right hand side vector.
   * @param iterations       Upon return this argument holds the number of used iterations.
   */
  template<typename T>
  inline void jacobi(
                     ublas::compressed_matrix<T>  const & A
                     , ublas::vector<T>            & x
                     , ublas::vector<T>         const  & b
                     , size_t & iterations
                     )
  {
    size_t const max_iterations = 15u;
    T      const epsilon        = boost::numeric_cast<T>( 10e-4 );
    
    jacobi(A,x,b,max_iterations, epsilon, iterations);
  }
    
} // namespace big

// BIG_JACOBI_H
#endif
