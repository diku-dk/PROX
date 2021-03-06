#ifndef HYPER_PARAMS_H
#define HYPER_PARAMS_H

#include <mesh_array.h>

#include <tiny_is_number.h>
#include <tiny_is_finite.h>

#include <util_log.h>

#include <cassert>

namespace hyper
{    
  /**
   * Paramters for controlling hyper simulation.
   *
   * @tparam   Math types policy
   */
  template< typename MT >
  class Params
  {
  public:

    typedef enum {
      semi_implicit_type
      , adaptive_type
      , implicit_type
    } time_step_method_type;
    
  protected:

    typedef typename MT::real_type    T;
    typedef typename MT::value_traits VT;

  protected:

    T       m_time_step;                   ///< The fixed time stepping to use when time-integrating the dynamics.
    bool    m_use_cfl_condition;           ///< Boolean flag for turning on/off CFL conditions
    T       m_adaptive_min_dt;             ///< The minimum allowed time step size when using adaptive time-stepping.
    T       m_adaptive_max_dt;             ///< The maximum allowed time step size when using adaptive time-stepping.
    T       m_adaptive_halving_tolerance;  ///< The havling tolerance, the maximum allowed error (accuracy).
    size_t  m_adaptive_doubling_count;     ///< The number of consecutive time steps to take before trying ot increase time step-size.
    T       m_envelope;                    ///< Procentage of scale of smallest object size to be used as collision envelope

    mesh_array::TetGenSettings m_tetgen_settings;              ///< Tetgen settings
    time_step_method_type      m_time_step_method_type;

  public:

    T      const & get_envelope()      const { return this->m_envelope;           }

  public:

    void set_envelope(T const & value)              { this->m_envelope         = value; }

  public:

    mesh_array::TetGenSettings const & tetgen_settings() const
    {
      return this->m_tetgen_settings;
    }

    void set_tetgen_settings(mesh_array::TetGenSettings const & value)
    {
      this->m_tetgen_settings = value;
    }

    void set_tetgen_quality_ratio(T const & value)
    {
      assert( value >= VT::zero() || !"set_tetgen_quality_ratio(): illegal value");

      this->m_tetgen_settings.m_quality_ratio = value;
    }

    void set_tetgen_maximum_volume(T const & value)
    {
      assert( value >= VT::zero() || !"set_tetgen_maximum_volume(): illegal value");

      this->m_tetgen_settings.m_maximum_volume = value;
    }

    void set_tetgen_quiet_output(bool const & value)
    {
      this->m_tetgen_settings.m_quiet_output = value;
    }

    void set_tetgen_suppress_splitting(bool const & value)
    {
      this->m_tetgen_settings.m_suppress_splitting = value;
    }

    void set_tetgen_verify_input(bool const & value)
    {
      this->m_tetgen_settings.m_verify_input = value;
    }

    void set_tetgen_filename(std::string const & value)
    {
      this->m_tetgen_settings.m_filename = value;
    }

  public:

    T const & time_step() const
    {
      return this->m_time_step;
    }

    void set_time_step(T const & value)
    {
      assert( value >= VT::zero() || !"set_time_step(): illegal value");

      this->m_time_step = value;
    }

    bool const & use_cfl_condition() const
    {
      return this->m_use_cfl_condition;
    }

    void set_use_cfl_condition(bool const & value)
    {
      assert( value >= VT::zero() || !"set_use_cfl_condition(): illegal value");

      this->m_use_cfl_condition = value;
    }


  public:

    T      const & adaptive_min_dt()            const    {      return this->m_adaptive_min_dt;               }
    T      const & adaptive_max_dt()            const    {      return this->m_adaptive_max_dt;               }
    T      const & adaptive_halving_tolerance() const    {      return this->m_adaptive_halving_tolerance;    }
    size_t const & adaptive_doubling_count()    const    {      return this->m_adaptive_doubling_count;       }

    void set_adaptive_min_dt(T const & value)
    {
      assert( value >= VT::zero() || !"set_adaptive_min_dt(): illegal value");
      this->m_adaptive_min_dt = value;
    }

    void set_adaptive_max_dt(T const & value)
    {
      assert( value >= VT::zero() || !"set_adaptive_max_dt(): illegal value");
      this->m_adaptive_max_dt = value;
    }

    void set_adaptive_halving_tolerance(T const & value)
    {
      assert( value > VT::zero() || !"set_adaptive_halving_tolerance(): illegal value");
      this->m_adaptive_halving_tolerance = value;
    }

    void set_adaptive_doubling_count(size_t const & value)
    {
      assert( value > 0u || !"set_adaptive_doubling_count(): illegal value");
      this->m_adaptive_doubling_count = value;
    }

  public:

    time_step_method_type get_time_step_method_type()const
    {
      return this->m_time_step_method_type;
    }

    void set_time_step_method_type(time_step_method_type const & value )
    {
      this->m_time_step_method_type = value;
    }

    void set_time_step_method_type(std::string const & value )
    {
      if(value.compare("semi_implicit")==0)
      {
        this->m_time_step_method_type = semi_implicit_type;
      }
      else if(value.compare("adaptive")==0)
      {
        this->m_time_step_method_type = adaptive_type;
      }
      else if(value.compare("implicit")==0)
      {
        this->m_time_step_method_type = implicit_type;
      }
      else
      {
        util::Log logging;

        logging << "unknown time step method type =" << value << util::Log::newline();
      }
    }

  public:

    Params()
    : m_time_step(VT::numeric_cast(0.01))
    , m_use_cfl_condition(false)
    , m_adaptive_min_dt(VT::numeric_cast(0.001))
    , m_adaptive_max_dt(VT::numeric_cast(0.01))
    , m_adaptive_halving_tolerance(VT::numeric_cast(0.000001))
    , m_adaptive_doubling_count(5u)
    , m_envelope(VT::numeric_cast(0.01))
    , m_tetgen_settings( mesh_array::tetgen_quality_settings() )
    , m_time_step_method_type(semi_implicit_type)
    {}

    virtual ~Params()
    {
      this->clear();
    }

  public:

    void clear()
    {
      this->m_time_step         = VT::numeric_cast(0.01);
      this->m_use_cfl_condition = false;

      this->m_adaptive_min_dt              = VT::numeric_cast(0.001);
      this->m_adaptive_max_dt              = VT::numeric_cast(0.01);
      this->m_adaptive_halving_tolerance   = VT::numeric_cast(0.000001);
      this->m_adaptive_doubling_count      =  5u;

      this->m_envelope          = VT::numeric_cast(0.01);
      this->m_tetgen_settings   = mesh_array::tetgen_quality_settings();

      this->m_time_step_method_type = semi_implicit_type;
    };

  };

} // namespace hyper

// HYPER_PARAMS_H
#endif 
