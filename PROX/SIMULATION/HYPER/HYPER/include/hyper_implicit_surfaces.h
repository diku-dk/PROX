#ifndef HYPER_IMPLICIT_SURFACES_H
#define HYPER_IMPLICIT_SURFACES_H

#include <tiny.h>

namespace hyper
{

  template<typename MT>
  class ImplicitSurfaceFunctor
  {
  public:

    typedef typename MT::real_type       T;
    typedef typename MT::vector3_type    V;
    typedef typename MT::value_traits    VT;

    virtual T operator()(V const & p) const = 0;

  };


  template<typename MT>
  class BoxSurface
  : public ImplicitSurfaceFunctor<MT>
  {
  public:

    typedef typename MT::real_type       T;
    typedef typename MT::vector3_type    V;
    typedef typename MT::value_traits    VT;

  protected:

    V m_position;
    V m_dimension;

  public:

    V const & position() const  { return this->m_position;  }
    V       & position()        { return this->m_position;  }

    V const & dimension() const { return this->m_dimension; }
    V       & dimension()       { return this->m_dimension; }

  public:

    BoxSurface()
    : m_position(V::zero())
    , m_dimension(V::make( VT::one(), VT::one(), VT::one()) )
    {}

    BoxSurface(V const & p, V const & d)
    : m_position( p )
    , m_dimension( d )
    {}

  public:

    T operator()(V const & p) const
    {
      V    const & c         = this->m_position;
      V    const & s         = this->m_dimension*VT::half();;
      V    const   d         = tiny::abs( p-c );
      V    const   delta     = d-s;
      //T    const   min_delta = tiny::min(delta);
      T    const   max_delta = tiny::max(delta);
      bool const   inside    = max_delta < VT::zero();

      if (inside)
      {
        return max_delta;
      }

      T const proj_x = delta(0) > VT::zero() ? s(0) : d(0);
      T const proj_y = delta(1) > VT::zero() ? s(1) : d(1);
      T const proj_z = delta(2) > VT::zero() ? s(2) : d(2);

      V const proj   = V::make(proj_x,proj_y, proj_z);

      return norm( d- proj );
    }
  };

  template<typename MT>
  BoxSurface<MT> make_box( typename MT::vector3_type const & position, typename MT::vector3_type const & dimension)
  {
    return BoxSurface<MT>(position, dimension);
  }

} // namespace hyper

// HYPER_IMPLICIT_SURFACES_H
#endif 
