#ifndef HYPER_SAINT_VERNANT_KIRCHHOFF_H
#define HYPER_SAINT_VERNANT_KIRCHHOFF_H

#include <hyper_strain_tensors.h>
#include <hyper_constitutive_equation.h>

#include <tiny_is_number.h>
#include <tiny_is_finite.h>

namespace hyper
{

  template<typename MT>
  class SaintVernantKirchhoff
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

    SaintVernantKirchhoff()
      : ConstitutiveEquation<MT>(  )
    {}

    SaintVernantKirchhoff( material_type const & material)
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
      // The strain energy is defined as
      //
      //  psi(E) = lambda/2 tr(E)^2 + mu E : E =  lambda/2 tr(E)^2 + mu tr(E^2)
      //
      //  where E is the Green strain tensor
      //
      //    E =  1/2 (C - I)
      //
      // where I is the identify tensor and C is the right Cauchy strain tensor
      //
      //    C = F^T F
      //
      // lambda and mu are the Lamé coefficients.
      //
      // The stress-strain relationship for 2nd Piola Kirchhoff stress tensor is given by
      //
      //  S = d psi / d E
      //
      // Resulting in the consititutive model
      //
      //  S = lambda tr(E) I  + 2 mu E
      //
      M const E      =  green_strain_tensor(F);
      T const lambda = this->m_material.lambda();
      T const mu     = this->m_material.mu();

      return (lambda * tiny::trace(E)) * M::identity() + (VT::two() * mu) * E;
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
    M sigma(M const & F) const
    {
      T const j = tiny::det( F );

      assert( is_number(j) || !"sigma(): j was not a number");
      assert( is_finite(j) || !"sigma(): j was not finite");

      return (VT::one()/j)* tiny::prod(F, tiny::prod( this->S(F), tiny::trans(F) ) );
    }

    std::string name() const
    {
      return "Saint Venant Kirchhoff";
    }
  };

}// namespace hyper

// HYPER_SAINT_VERNANT_KIRCHHOFF_H
#endif
