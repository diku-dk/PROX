#ifndef PROX_SUBSOLVER_NORMAL_H
#define PROX_SUBSOLVER_NORMAL_H

namespace prox
{
  /**
   * Normal sub solver functor.
   */
  template<typename T>
  class NormalSubSolver
  {
  public:
    
    virtual void operator()(
                    T const & z_n
                    , T & lambda_n
                            ) const = 0;
    
  };
    
} //namespace prox

// PROX_SUBSOLVER_NORMAL_H
#endif
