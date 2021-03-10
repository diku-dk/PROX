#ifndef PROX_SUBSOLVER_FRICTION_H
#define PROX_SUBSOLVER_FRICTION_H

namespace prox
{
  /**
   * A friction sub solver functor.
   */  
  template<typename T>
  class FrictionSubSolver
  {    
  public:
    
    virtual void operator()(
                    T const & z_s
                    , T const & z_t
                    , T const & z_tau
                    , T const & mu_s
                    , T const & mu_t
                    , T const & mu_tau
                    , T const & lambda_n
                    , T & lambda_s
                    , T & lambda_t
                    , T & lambda_tau
                    ) const = 0;
    
  };
  
} //namespace prox

// PROX_SUBSOLVER_FRICTION_H
#endif
