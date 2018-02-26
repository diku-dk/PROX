#ifndef BIG_SYMMETRIC_GAUSS_SEIDEL_H
#define BIG_SYMMETRIC_GAUSS_SEIDEL_H

#include <big_forward_gauss_seidel.h>
#include <big_backward_gauss_seidel.h>

#include <stdexcept>

namespace big
{
  
  /**
   * Symmetric Gauss-Seidel Iteration.
   * This function implements a single iteration of a symmetric
   * Gauss Seidel method. That is it performs a forward Gauss-Seidel
   * iteration followed by a backward Gauss-Seidel iteration.
   *
   * @param A       A matrix.
   * @param x       Upon return this argument holds the new value of iterate x.
   * @param b       The right hand side vector.
   *
   */
  template<typename T>
  inline void symmetric_gauss_seidel(
                                     boost::numeric::ublas::compressed_matrix<T> const & A
                                     , boost::numeric::ublas::vector<T>       & x
                                     , boost::numeric::ublas::vector<T> const & b
                                     )
  {
    forward_gauss_seidel( A, x, b );
    backward_gauss_seidel( A, x, b );
  }
  
  /**
   * Symmetric Gauss-Seidel Solver.
   * This function is capable of performing several iterations of the symmetric Gauss Seidel iteration.
   *
   * @param A                A matrix.
   * @param x                At invokation this argument holds the initial value of x^0, Upon return this argument holds the new value of iterate.
   * @param b                The right hand side vector.
   * @param max_iterations   The maximum number of iterations that is allowed.
   * @param epsilon          A stopping threshold (currently not used).
   * @param iterations       Upon return this argument holds the number of used iterations.
   */
  template<typename T>
  inline void symmetric_gauss_seidel(
                                     boost::numeric::ublas::compressed_matrix<T> const & A
                                     , boost::numeric::ublas::vector<T>       & x
                                     , boost::numeric::ublas::vector<T> const & b
                                     , size_t                           const & max_iterations
                                     , size_t                                 & iterations
                                     )
  {
    if(max_iterations < 1u)
      throw std::invalid_argument("symmetric_gauss_seidel(): max_iterations must be a positive number");
    
    iterations = 0u;
    while(iterations<max_iterations)
    {
      ++iterations;
      symmetric_gauss_seidel(A,x,b);
    }
  }
    
} // end of namespace big

// BIG_SYMMETRIC_GAUSS_SEIDEL_H
#endif
