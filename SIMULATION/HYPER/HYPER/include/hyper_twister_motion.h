#ifndef HYPER_TWISTER_MOTION_H
#define HYPER_TWISTER_MOTION_H

#include <hyper_scripted_motion.h>
#include <hyper_modifiers.h>

#include <tiny.h>

#include <cassert>

namespace hyper
{

  template<typename MT>
  class TwisterMotion
  : public ScriptedMotion<MT>
  {
  public:

    typedef typename MT::real_type       T;
    typedef typename MT::vector3_type    V;
    typedef typename MT::value_traits    VT;

  protected:

    V m_axis;
    T m_period;

  public:

    V const & axis() const { return this->m_axis; }
    V       & axis()       { return this->m_axis; }

    T const & period() const { return this->m_period; }
    T       & period()       { return this->m_period; }

  public:

    TwisterMotion()
    : m_axis( V::k() )
    , m_period( VT::two() )
    {}

    TwisterMotion(V const & a, T const & p)
    : m_axis( tiny::unit( a ) )
    , m_period( p )
    {
      assert( this->m_period > VT::zero() || !"TwisterMotion(): period must be positive");
    }


  public:

    void compute( Body<MT> & body, typename MT::real_type const & time)
    {
      unsigned int n = time / this->m_period;

      T const radian = ((n%2)==0) ? VT::pi() : -VT::pi();

      twist(
            body
            , this->m_axis
            , radian/VT::numeric_cast(180.0)
            );
    }

  };

} // namespace hyper

// HYPER_TWISTER_MOTION_H
#endif 
