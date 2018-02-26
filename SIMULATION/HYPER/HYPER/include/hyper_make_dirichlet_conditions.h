#ifndef HYPER_MAKE_DIRICHLET_CONDITIONS_H
#define HYPER_MAKE_DIRICHLET_CONDITIONS_H

#include <hyper_dirichlet_info.h>
#include <hyper_implicit_surfaces.h>
#include <hyper_body.h>
#include <hyper_tags.h>

#include <tiny.h>
#include <mesh_array.h>

#include <cassert>
#include <vector>

namespace hyper
{

  template<typename MT>
  inline void make_dirichlet_conditions(
                                  Body<MT> & body
                                  , ImplicitSurfaceFunctor<MT> const & phi
                                  , use_deformed_coordinates const & /*tag*/
                                  )
  {
    typedef typename MT::vector3_type    V;
    typedef typename MT::value_traits    VT;

    body.m_dirichlet_conditions.clear();

    for(size_t k = 0u; k < body.m_X.size(); ++k)
    {
      V const p = V::make( body.m_X(k), body.m_Y(k), body.m_Z(k) );

      if( phi(p) <= VT::zero() )
        body.m_dirichlet_conditions.push_back( make_dirichlet_info<MT>(k, V::zero() ));

    }
  }

  template<typename MT>
  inline void make_dirichlet_conditions(
                                  Body<MT> & body
                                  , ImplicitSurfaceFunctor<MT> const & phi
                                  , use_material_coordinates const & /*tag*/
                                  )
  {
    typedef typename MT::vector3_type    V;
    typedef typename MT::value_traits    VT;

    body.m_dirichlet_conditions.clear();

    for(size_t k = 0u; k < body.m_X0.size(); ++k)
    {
      V const p = V::make( body.m_X0(k), body.m_Y0(k), body.m_Z0(k) );

      if( phi(p) <= VT::zero() )
        body.m_dirichlet_conditions.push_back( make_dirichlet_info<MT>(k, V::zero() ));

    }
  }


  template<typename MT>
  inline void make_dirichlet_conditions(
                                  Body<MT> & body
                                  , ImplicitSurfaceFunctor<MT> const & phi
  )
  {
    make_dirichlet_conditions(body, phi, use_material_coordinates() );
  }


} // namespace hyper

// HYPER_MAKE_DIRICHLET_CONDITIONS_H
#endif 
