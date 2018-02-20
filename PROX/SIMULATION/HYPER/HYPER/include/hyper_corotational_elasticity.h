#ifndef HYPER_COROTATIONAL_ELASTICITY_H
#define HYPER_COROTATIONAL_ELASTICITY_H

#include <hyper_strain_tensors.h>
#include <hyper_constitutive_equation.h>

#include <tiny_polar_decomposition3x3.h>

namespace hyper
{

  template<typename MT>
  class CorotationalElasticity
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

    CorotationalElasticity()
      : ConstitutiveEquation<MT>(  )
    {}

    CorotationalElasticity( material_type const & material)
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
      //
      //  P = lambda tr(S-I) R  + 2 mu (F-R)
      //
      M R;
      M S;

      tiny::polar_decomposition_eigen(F, R, S);

      T const lambda = this->m_material.lambda();
      T const mu     = this->m_material.mu();
      M const I      = M::identity();

      return (lambda * tiny::trace(S-I)) * R  + (VT::two() * mu) * (F-R);
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
      M const invF = tiny::inverse(F);

      return  tiny::prod( invF, this->P(F) );
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
      T const j = tiny::det( F );

      assert( is_number(j) || !"sigma(): j was not a number");
      assert( is_finite(j) || !"sigma(): j was not finite");

      //return (VT::one()/j)* tiny::prod(F, tiny::prod( this->S(F), tiny::trans(F) ) );
      return (VT::one()/j)* tiny::prod( this->P(F), tiny::trans(F) );
    }
    
    std::string name() const
    {
      return "Corotational Elasticity";
    }
    
  };

}// namespace hyper

// HYPER_COROTATIONAL_ELASTICITY_H
#endif
