#ifndef HYPER_STRAIN_TENSORS_H
#define HYPER_STRAIN_TENSORS_H

#include <tiny.h>

namespace hyper
{

  template<typename M>
  inline M right_cauchy_green_strain_tensor(M const & F)
  {
    return tiny::prod( trans(F), F);
  }

  template<typename M>
  inline M left_cauchy_green_strain_tensor(M const & F)
  {
    return tiny::prod( F, trans(F));
  }

  template<typename M>
  inline M green_lagrangian_strain_tensor(M const & F)
  {
    typedef typename M::value_traits   VT;

    M const C = right_cauchy_green_strain_tensor(F);
    M const I = M::identity();

    return  VT::half() * ( C - I );
  }

  template<typename M>
  inline M euler_almansi_strain_tensor(M const & F)
  {
    typedef typename M::value_traits   VT;

    M const I = M::identity();
    M const B = left_cauchy_green_strain_tensor(F);

    return VT::half()*(M::identity() - inverse(B));
  }

  template<typename M>
  inline M cauchy_strain_tensor(M const & F)
  {
    typedef typename M::value_traits   VT;

    M const I = M::identity();
    M const varepsilon = VT::half() * (F + tiny::trans(F)) - I;

    return varepsilon;
  }

  template<typename M>
  inline M right_cauchy_strain_tensor(M const & F) { return right_cauchy_green_strain_tensor(F);  }

  template<typename M>
  inline M left_cauchy_strain_tensor(M const & F)  { return left_cauchy_green_strain_tensor(F);   }

  template<typename M>
  inline M green_strain_tensor(M const & F)        { return green_lagrangian_strain_tensor(F);    }

  template<typename M>
  inline M euler_strain_tensor(M const & F)        {  return euler_almansi_strain_tensor(F);      }

}// namespace hyper

// HYPER_STRAIN_TENSORS_H
#endif
