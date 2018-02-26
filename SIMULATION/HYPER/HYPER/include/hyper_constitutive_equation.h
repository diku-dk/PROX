#ifndef HYPER_CONSTITUTIVE_EQUATION_H
#define HYPER_CONSTITUTIVE_EQUATION_H

#include <hyper_material_parameters.h>

#include <string>

namespace hyper
{

  /**
   * Base class for defining constitutive relations
   * between stress and strain fields.
   */
  template<typename MT>
  class ConstitutiveEquation
  {
  public:

    typedef typename MT::real_type            T;
    typedef typename MT::vector3_type         V;
    typedef typename MT::value_traits         VT;
    typedef typename MT::matrix3x3_type       M;

    typedef          MaterialParameters<T>    material_type;
    
  protected:

    material_type  m_material;

  public:

    material_type const & material() const { return this->m_material;   }
    material_type       & material()       { return this->m_material;   }

  public:

    ConstitutiveEquation()
    : m_material( make_default_material<T>() )
    {}

    ConstitutiveEquation( material_type const & material)
    : m_material( material )
    {}

    virtual ~ConstitutiveEquation(){}

  public:

    /**
     * Compute Cauchy Stress tensor
     *
     * @param F    The deformation gradient, F = du/dX where u is displacement
     *             field u = x - X and x is current spatial coordinates and X
     *             is the material coordinates (reference coordinate system)
     *
     * @return     The Cauchy stress tensor.
     */
    virtual M sigma(M const & F) const = 0;

    /**
     * Get human readable name for constitutive equation. Useful for making out text/error messages.
     *
     * @return    The name of the constitutive model.
     */
    virtual std::string name() const = 0;

  };


}// namespace hyper

// HYPER_CONSTITUTIVE_EQUATION_H
#endif
