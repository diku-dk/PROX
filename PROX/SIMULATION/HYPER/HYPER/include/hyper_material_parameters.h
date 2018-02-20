#ifndef HYPER_MATERIAL_PARAMETERS_H
#define HYPER_MATERIAL_PARAMETERS_H

#include <tiny_value_traits.h>

#include <util_log.h>

#include <string>
#include <cassert>

namespace hyper
{

  template<typename T>
  class MaterialParameters
  {
  protected:

    typedef tiny::ValueTraits<T>  VT;

    T    m_young_modulus;
    T    m_poisson_ratio;
    T    m_mass_density;
    T    m_mass_damping;
    T    m_stiffness_damping;
    T    m_viscous_damping;

    bool m_lumped;

  public:

    T const & E()     const { return this->m_young_modulus;     }
    T       & E()           { return this->m_young_modulus;     }
    T const & nu()    const { return this->m_poisson_ratio;     }
    T       & nu()          { return this->m_poisson_ratio;     }
    T const & rho()   const { return this->m_mass_density;      }
    T       & rho()         { return this->m_mass_density;      }
    T const & alpha() const { return this->m_mass_damping;      }
    T       & alpha()       { return this->m_mass_damping;      }
    T const & beta()  const { return this->m_stiffness_damping; }
    T       & beta()        { return this->m_stiffness_damping; }
    T const & c()     const { return this->m_viscous_damping;   }
    T       & c()           { return this->m_viscous_damping;   }

    /**
     * Converts Young Modulus and Poisson ratio into Lamé's first parameter.
     *
     * @return     Lamé's first parameter
     */

    T lambda() const
    {
      assert( m_poisson_ratio >= VT::zero() || !"lambda(): Illegal Poisson ratio value" );
      assert( m_poisson_ratio > -VT::one()  || !"lambda(): Illegal Poisson ratio value" );
      assert( m_young_modulus >  VT::zero() || !"lambda(): Illegal Young's Modulus value" );

      return (m_young_modulus*m_poisson_ratio) / ( (VT::one() + m_poisson_ratio)*(VT::one() - VT::two()*m_poisson_ratio) );
    }

    /**
     * Converts Young Modulus and Poisson ratio into Lamé's second parameter.
     *
     * @param E    Youngs Modulus
     * @param nu   Poisson Ratio
     *
     * @return     Lamé's second parameter
     */
    T mu() const
    {
      assert( m_poisson_ratio >= VT::zero() || !"mu(): Illegal Poisson ratio value" );
      assert( m_poisson_ratio > -VT::one()  || !"mu(): Illegal Poisson ratio value" );
      assert( m_young_modulus >  VT::zero() || !"mu(): Illegal Young's Modulus value" );

      return m_young_modulus / ( VT::two()*(VT::one()+m_poisson_ratio) );
    }

    bool const & lumped() const { return this->m_lumped; }
    bool       & lumped()       { return this->m_lumped; }

  public:

    MaterialParameters()
    : m_young_modulus(     VT::numeric_cast(1000.0) )
    , m_poisson_ratio(     VT::numeric_cast(0.3   ) )
    , m_mass_density(      VT::numeric_cast(10.0  ) )
    , m_mass_damping(      VT::zero()   )
    , m_stiffness_damping( VT::zero()   )
    , m_viscous_damping(   VT::zero()   )
    , m_lumped( true )
    {}

  };

  template<typename T>
  inline MaterialParameters<T> make_default_material()
  {
    typedef tiny::ValueTraits<T>  VT;

    MaterialParameters<T> MP;

    MP.E()      = VT::numeric_cast(1000.0);
    MP.nu()     = VT::numeric_cast(0.3);
    MP.rho()    = VT::numeric_cast(10.0);
    MP.alpha()  = VT::zero();
    MP.beta()   = VT::zero();
    MP.c()      = VT::numeric_cast(0.0001);
    MP.lumped() = true;

    return MP;
  }

  template<typename T>
  inline MaterialParameters<T> make_cartilage_material()
  {
    typedef tiny::ValueTraits<T>  VT;

    MaterialParameters<T> MP;

    MP.E()      = VT::numeric_cast(0.69e6);
    MP.nu()     = VT::numeric_cast(0.018);
    MP.rho()    = VT::numeric_cast(1000.0);
    MP.alpha()  = VT::zero();
    MP.beta()   = VT::zero();
    MP.c()      = VT::numeric_cast(0.001);
    MP.lumped() = true;

    return MP;
  }

  template<typename T>
  inline MaterialParameters<T> make_cortical_bone_material()
  {
    typedef tiny::ValueTraits<T>  VT;

    MaterialParameters<T> MP;

    MP.E()      = VT::numeric_cast(16.16e9);
    MP.nu()     = VT::numeric_cast(0.33);
    MP.rho()    = VT::numeric_cast(1600.0);
    MP.alpha()  = VT::zero();
    MP.beta()   = VT::zero();
    MP.c()      = VT::numeric_cast(0.001);
    MP.lumped() = true;

    return MP;
  }

  template<typename T>
  inline MaterialParameters<T> make_cancellous_bone_material()
  {
    typedef tiny::ValueTraits<T>  VT;

    MaterialParameters<T> MP;

    MP.E()      = VT::numeric_cast(452e6);
    MP.nu()     = VT::numeric_cast(0.3);
    MP.rho()    = VT::numeric_cast(1600.0);
    MP.alpha()  = VT::zero();
    MP.beta()   = VT::zero();
    MP.c()      = VT::numeric_cast(0.001);
    MP.lumped() = true;

    return MP;
  }

  template<typename T>
  inline MaterialParameters<T> make_rubber_material()
  {
    typedef tiny::ValueTraits<T>  VT;

    MaterialParameters<T> MP;

    MP.E()      = VT::numeric_cast(0.01e9);
    MP.nu()     = VT::numeric_cast(0.48);
    MP.rho()    = VT::numeric_cast(1050.0);
    MP.alpha()  = VT::zero();
    MP.beta()   = VT::zero();
    MP.c()      = VT::numeric_cast(0.001);
    MP.lumped() = true;

    return MP;
  }

  template<typename T>
  inline MaterialParameters<T> make_concrete_material()
  {
    typedef tiny::ValueTraits<T>  VT;

    MaterialParameters<T> MP;

    MP.E()      = VT::numeric_cast(30e9);
    MP.nu()     = VT::numeric_cast(0.20);
    MP.rho()    = VT::numeric_cast(2320.0);
    MP.alpha()  = VT::zero();
    MP.beta()   = VT::zero();
    MP.c()      = VT::numeric_cast(0.001);
    MP.lumped() = true;

    return MP;
  }

  template<typename T>
  inline MaterialParameters<T> make_cobber_material()
  {
    typedef tiny::ValueTraits<T>  VT;

    MaterialParameters<T> MP;

    MP.E()      = VT::numeric_cast(125e9);
    MP.nu()     = VT::numeric_cast(0.35);
    MP.rho()    = VT::numeric_cast(8900.0);
    MP.alpha()  = VT::zero();
    MP.beta()   = VT::zero();
    MP.c()      = VT::numeric_cast(0.001);
    MP.lumped() = true;

    return MP;
  }

  template<typename T>
  inline MaterialParameters<T> make_steel_material()
  {
    typedef tiny::ValueTraits<T>  VT;

    MaterialParameters<T> MP;

    MP.E()      = VT::numeric_cast(210e9);
    MP.nu()     = VT::numeric_cast(0.31);
    MP.rho()    = VT::numeric_cast(7800.0);
    MP.alpha()  = VT::zero();
    MP.beta()   = VT::zero();
    MP.c()      = VT::numeric_cast(0.001);
    MP.lumped() = true;

    return MP;
  }

  template<typename T>
  inline MaterialParameters<T> make_aluminium_material()
  {
    typedef tiny::ValueTraits<T>  VT;

    MaterialParameters<T> MP;

    MP.E()      = VT::numeric_cast(72e9);
    MP.nu()     = VT::numeric_cast(0.34);
    MP.rho()    = VT::numeric_cast(2700.0);
    MP.alpha()  = VT::zero();
    MP.beta()   = VT::zero();
    MP.c()      = VT::numeric_cast(0.001);
    MP.lumped() = true;

    return MP;
  }

  template<typename T>
  inline MaterialParameters<T> make_glass_material()
  {
    typedef tiny::ValueTraits<T>  VT;

    MaterialParameters<T> MP;

    MP.E()      = VT::numeric_cast(50e9);
    MP.nu()     = VT::numeric_cast(0.18);
    MP.rho()    = VT::numeric_cast(2190.0);
    MP.alpha()  = VT::zero();
    MP.beta()   = VT::zero();
    MP.c()      = VT::numeric_cast(0.001);
    MP.lumped() = true;

    return MP;
  }

  template<typename T>
  inline MaterialParameters<T> make_material(std::string const & name)
  {
    if( name.compare("default")==0)
      return make_default_material<T>();
    if( name.compare("cartilage")==0)
      return make_cartilage_material<T>();
    if( name.compare("cortical_bone")==0)
      return make_cortical_bone_material<T>();
    if( name.compare("cancellous_bone")==0)
      return make_cancellous_bone_material<T>();
    if( name.compare("rubber")==0)
      return make_rubber_material<T>();
    if( name.compare("concrete")==0)
      return make_concrete_material<T>();
    if( name.compare("cobber")==0)
      return make_cobber_material<T>();
    if( name.compare("steel")==0)
      return make_steel_material<T>();
    if( name.compare("aluminium")==0)
      return make_aluminium_material<T>();
    if( name.compare("glass")==0)
      return make_glass_material<T>();

    util::Log logging;

    logging << "hyper::make_material(): unrecognized material" << name << " supported names are" << util::Log::newline();
    logging << "\t" << "default"         << util::Log::newline();
    logging << "\t" << "cartilage"       << util::Log::newline();
    logging << "\t" << "cortical_bone"   << util::Log::newline();
    logging << "\t" << "cancellous_bone" << util::Log::newline();
    logging << "\t" << "rubber"          << util::Log::newline();
    logging << "\t" << "concrete"        << util::Log::newline();
    logging << "\t" << "cobber"          << util::Log::newline();
    logging << "\t" << "steel"           << util::Log::newline();
    logging << "\t" << "aluminium"       << util::Log::newline();
    logging << "\t" << "glass"           << util::Log::newline();

    return make_default_material<T>();
  }

}// namespace hyper

// HYPER_MATERIAL_PARAMETERS_H
#endif
