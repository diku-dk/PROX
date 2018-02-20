#ifndef BIG_FORWARD_GAUSS_SEIDEL_H
#define BIG_FORWARD_GAUSS_SEIDEL_H

#include <big_types.h>
#include <big_is_number.h>
#include <big_value_traits.h>

#include <stdexcept>
#include <cassert>
#include <cmath>   // needed for std::fabs

namespace big
{
  
  /**
   * Forward Gauss-Seidel Iteration.
   * This function performs a single Gauss-Seidel iteration, that
   * is it solves x^{k+1} = inv(D+L)(b - U x^k).
   *
   * @param A    The matrix.
   * @param x    When called this parameter holds the current
   *             value of the iterate, upon return this arugment
   *             holds the new value of the iterate.
   * @param b    The right hand side vector.
   */
  template<typename T>
  inline void forward_gauss_seidel(
                                   boost::numeric::ublas::compressed_matrix<T> const & A
                                   , boost::numeric::ublas::vector<T>       & x
                                   , boost::numeric::ublas::vector<T> const & b
                                   )
  {
    typedef big::ValueTraits<T>                        value_traits;
    
    using std::fabs;
    
    if(A.size1() <= 0 || A.size2() <= 0)
      throw std::invalid_argument("forward_gauss_seidel(): A was empty");
    
    if(b.size() != A.size1())
      throw std::invalid_argument("forward_gauss_seidel(): The size of b must be the same as the number of rows in A");
    
    if(x.size() != A.size2())
      throw std::invalid_argument("forward_gauss_seidel(): The size of x must be the same as the number of columns in A");
    
    size_t const N = A.filled1() - 1;
    
    for (size_t row = 0u; row < N; ++row)
    {
      size_t const begin = A.index1_data()[row];
      size_t const end   = A.index1_data()[row + 1];
      
      T         sum   = b(row);
      T const & diag  = A(row,row);
      
      for (size_t j = begin; j < end; ++j)
      {
        size_t  const & col  = A.index2_data()[j];
        T       const & A_ij = A.value_data()[j];
        
        assert( ( col >= 0u && col< A.size2() )       || !"forward_gauss_seidel(): column index were out of range");
        assert( is_number( A_ij )                     || !"forward_gauss_seidel(): A_ij value was not a number?");
        
        sum -= A_ij * x( col );
      }
      assert( is_number( sum )                  || !"forward_gauss_seidel(): sum value was not a number?");
      assert( is_number( diag )                 || !"forward_gauss_seidel(): diag value was not a number?");
      
      // 2007-06-10 kenny: Yikes how should we handle a diagonal zero-value?
      assert( fabs(diag) > value_traits::zero() || !"forward_gauss_seidel(): Diagonal were zero");
      x(row) += sum / diag;
      
      assert( is_number( x(row ) )              || !"forward_gauss_seidel(): updated value was not a number?");
    }
  }
  
  /**
   * Forward Gauss-Seidel Iteration.
   * This function performs a several Gauss-Seidel iterations, that
   * is it given x^0 it computes x^{k+1} = inv(D+L)(b - U x^k) for solve k>0.
   *
   * @param A                 The matrix.
   * @param x                 When called this parameter holds the current
   *                          value of the iterate, upon return this arugment
   *                          holds the new value of the iterate.
   * @param b                 The right hand side vector.
   * @param max_iterations    The maximum number of iterations that can be taken.
   * @param iterations        Upon return this argument holds the number of actual iterations that were used.
   */
  template<typename T>
  inline void forward_gauss_seidel(
                                   boost::numeric::ublas::compressed_matrix<T> const & A
                                   , boost::numeric::ublas::vector<T>       & x
                                   , boost::numeric::ublas::vector<T> const & b
                                   , size_t                           const & max_iterations
                                   , size_t                                 & iterations
                                   )
  {
    if(max_iterations < 1u)
      throw std::invalid_argument("forward_gauss_seidel(): max_iterations must be a positive number");
    
    iterations = 0u;
    while(iterations<max_iterations)
    {
      ++iterations;
      forward_gauss_seidel(A,x,b);
    }
  }
    
} // end of namespace big

// BIG_FORWARD_GAUSS_SEIDEL_H
#endif
