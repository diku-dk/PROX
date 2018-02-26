#ifndef BIG_MIN_MAP_NEWTON_H
#define BIG_MIN_MAP_NEWTON_H

#include <big_types.h>
#include <big_return_codes.h>
#include <big_value_traits.h>
#include <big_precision.h>

#include <big_prod.h>
#include <big_prod_add_rhs.h>
#include <big_prod_trans.h>

#include <big_absolute_convergence.h>
#include <big_stagnation.h>
#include <big_relative_convergence.h>
#include <big_project.h>
#include <big_svd.h>

#include <big_is_number.h>

#include <min_map/big_compute_jacobian.h>
#include <min_map/big_compute_minimum_map.h>
#include <min_map/big_compute_index_sets.h>

#include <cmath>
#include <stdexcept>
#include <cassert>

namespace big
{
  template < typename T>
  inline void min_map_newton(
                             ublas::compressed_matrix<T> const & A
                             , ublas::vector<T> const & b
                             , ublas::vector<T> & x
                             , size_t const & max_iter
                             , T      const & tol_abs
                             , T      const & tol_rel
                             , T      const & tol_stag
                             , size_t       & status
                             , size_t       & iter
                             , T            & err
                             , T      const & alpha
                             , T      const & beta
                             , ublas::vector<T> * profiling = 0
                             )
  {
    using std::fabs;
    using std::min;
    using std::max;
    
    typedef          ublas::compressed_matrix<T>      matrix_type;
    typedef          ublas::vector<T>                 vector_type;
    typedef          ValueTraits<T>                   VT;
    
    if(max_iter <= 0u)
      throw std::invalid_argument("max_iterations must be larger than zero");
    
    if(tol_abs < VT::zero() )
      throw std::invalid_argument("absolute_tolerance must be non-negative");
    
    if(tol_rel < VT::zero() )
      throw std::invalid_argument("relative_tolerance must be non-negative");
    
    if(tol_stag < VT::zero() )
      throw std::invalid_argument("stagnation_tolerance must be non-negative");
    
    if (beta >= VT::one() )
      throw std::invalid_argument("Illegal beta value");
    
    if (alpha <= VT::zero() )
      throw std::invalid_argument("Illegal alpha value");
    
    if(beta<=alpha)
      throw std::invalid_argument("beta must be larger than alpha");
    
    if(profiling == &x)
      throw std::logic_error("profiling must not point to x-vector");
    
    if(profiling == &b)
      throw std::logic_error("profiling must not point to b-vector");
    
    T  const TOO_TINY_STEP_LENGTH = ::boost::numeric_cast<T>(0.00001);  // 2012-09-30 Kenny: This is just a hardwired constant we could pass this an input param?
    
    status = OK;
    
    size_t const N = b.size();
    
    if(N==0u)
      return;
    
    if(x.size()!=N)
      throw std::invalid_argument("size of x-vector were different from b");
    
    //--- Initialization ------------------------------------------
    if(profiling)
    {
      (*profiling).resize( max_iter );
      (*profiling).clear();
    }
    
    err             = VT::infinity();
    iter            = 0u;
    status          = ITERATING;
    big::project( x, x);                          // x = max(0, x ), make sure x is feasible before we start

    
    //--- Start Newton Loop ----------------------------------------
    for (; iter < max_iter; ++iter)
    {
      vector_type y(N);
      big::prod_add_rhs(A, x, b, y);             // y = A x + b
      
      vector_type H(N);
      big::compute_minimum_map(y, x, H);         // H = min(x,y)
      
      T const old_err = err;
      
      err = ublas::inner_prod(H,H)*VT::half();   // err = theta(x) = 1/2 H^T H
      
      if(profiling)
        (*profiling)(iter) = err;
          
      if( big::relative_convergence(old_err, err, tol_rel) )
      {
        status = RELATIVE_CONVERGENCE;
        return;
      }
      
      if( big::absolute_convergence( err, tol_abs) )
      {
        status = ABSOLUTE_CONVERGENCE;
        return;
      }
      
      // Solve Newton Subsystem
      size_t cnt_active   = 0u;    // 2012-09-30 Kenny: I do not really use these two variables for anything in this version of the min-map Newton solver.. they coudld be dropped.
      size_t cnt_inactive = 0u;
      
      // Determine active/free index sets
      //
      //  A = { i | y <  x }
      //  F = { i | y >= x }
      ublas::vector<size_t> bitmask(N);
      big::compute_index_sets( y, x, bitmask, cnt_active, cnt_inactive ); // 2012-09-30 Kenny: To build the J-matrix we really do not need the bitmask vector, we could just pass the y and x vectors to the compute_jacobian matrix in stead. That would simplify things a little. 
      
      // Setup generalized Jacobian matrix
      //
      // J = | A_{AA} A_{AF} |
      //     |   0    I_{FF} |
      //
      // where I_{F,F} is identity matrix for set F
      //
      matrix_type J;
      big::compute_jacobian( A, bitmask, J );    // 2012-09-30 Kenny: Right now I just build the damn matrix... this is wasting
                                                 // memory as J is largely a copy of A just with simpler structure...
      vector_type rhs = -H;
      
      // solve J dx = - H
      vector_type dx(N);
      
      big::svd( J, dx, rhs );   // 2012-09-30 Kenny: Currently this is hardwired to SVD... this is of course insane...
                                // for fluid problems PCG should work fine (probably with a incomplete cholesky preconditioner)
      
      // Test if the search direction is smaller than numerical precision. That is if it is too close to zero.
      T const rho = big::working_precision<T>();   // 2012-09-30 Kenny: This is just a hardwired constant we could pass this an input param?
      if( ublas::norm_inf( dx ) <  rho )
      {
        status = STAGNATION;
        return;
      }
      
      // Test if we have dropped into a local minima if so we are stuck
      vector_type grad(N);
      big::prod_trans(J, H, grad);     // grad = J^T H
      
      if( ublas::norm_2( grad ) <  tol_abs )
      {
        status = LOCAL_MINIMA;
        return;
      }
      
      // Test if our search direction is a 'sufficient' descent direction
      if( ublas::inner_prod( grad, dx)  > - rho*ublas::inner_prod(dx, dx))
      {
        status = NON_DESCENT_DIRECTION;
        return;
      }
      
      //--- Armijo backtracking combined with a projected line-search ---------
      T tau           = VT::one();
      T const f_0     = err;
      vector_type x_k = x;
      
      while (true)
      {
        ublas::noalias( rhs ) = x + tau*dx;
        big::project( rhs, x_k);                         // x = max(0, x0 + tau * dx)
                
        big::prod_add_rhs(A, x_k, b, y);                 // y = A x + b
        big::compute_minimum_map(y, x_k, H);             // H = min(y,x)
        T const f_k = ublas::inner_prod(H,H)*VT::half(); // theta(tau) = 1/2 H^T H

        
        T grad_f_k = alpha * ublas::inner_prod(grad,x_k - x);
        
        //--- Perform Armijo codition to see if we got a sufficient decrease
        if ( f_k <= f_0 + tau*grad_f_k)
            break;
        
        //--- Test if time-step became too small
          if( tau*tau < TOO_TINY_STEP_LENGTH)
            break;
    
        tau = beta * tau;
      }
      
      x = x_k;

    }
    
    if( iter>=max_iter)
    {
      status = MAX_LIMIT;
    }
    
  }
  
} // namespace big

// BIG_MIN_MAP_NEWTON_H
#endif
