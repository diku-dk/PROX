#ifndef HYPER_NEO_HOOKEAN_H
#define HYPER_NEO_HOOKEAN_H

#include <hyper_strain_tensors.h>
#include <hyper_constitutive_equation.h>

namespace hyper
{

  template<typename MT>
  class NeoHookean
    : public ConstitutiveEquation<MT>
  {
  public:

    typedef          ConstitutiveEquation<MT>    base_type;
    typedef typename base_type::material_type    material_type;
    typedef typename MT::real_type               T;
    typedef typename MT::vector3_type            V;
    typedef typename MT::value_traits            VT;
    typedef typename MT::matrix3x3_type          M;

  public:

    NeoHookean()
      : ConstitutiveEquation<MT>(  )
    {}

    NeoHookean( material_type const & material)
      : ConstitutiveEquation<MT>( material )
    {}

  public:

    /**
     * Compute 1st Piola Kirchhoff Stress tensor
     *
     * @param F    The deformation gradient, F = du/dX where u is displacement
     *             field u = x - X and x is current spatial coordinates and X
     *             is the material coordinates (reference coordinate system)
     *
     * @return     The 1st Piola Kirchhoff stress tensor.
     */
    M P(M const & F) const
    {
      return tiny::prod(F , this->S(F) );
    }

    /**
     * Compute 2nd Piola Kirchhoff Stress tensor
     *
     * @param F    The deformation gradient, F = du/dX where u is displacement
     *             field u = x - X and x is current spatial coordinates and X
     *             is the material coordinates (reference coordinate system)
     *
     * @return     The 2nd Piola Kirchhoff stress tensor.
     */
    M S(M const & F) const
    {
      T const j    = tiny::det( F );
      M const invF = tiny::inverse(F);

      return  j * tiny::prod( invF, tiny::prod( this->sigma(F), tiny::trans(invF)) );
    }

    /**
     * Compute Cauchy Stress tensor
     *
     * @param F    The deformation gradient, F = du/dX where u is displacement
     *             field u = x - X and x is current spatial coordinates and X
     *             is the material coordinates (reference coordinate system)
     *
     * @return     The Cauchy stress tensor.
     */
    virtual M sigma(M const & F) const
    {
      // The consititutive model
      //
      //  sigma = lambda tr(e) I  + 2 mu e
      //
      M const e      =  cauchy_strain_tensor(F);
      T const lambda =  this->m_material.lambda();
      T const mu     =  this->m_material.mu();

      return (lambda * tiny::trace(e)) * M::identity() + (VT::two()*mu) * e;
    }
    
    std::string name() const
    {
      return "Neo Hookean";
    }
    
  };

}// namespace hyper

// HYPER_NEO_HOOKEAN_H
#endif
