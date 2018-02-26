#ifndef HYPER_MAKE_TRACTION_CONDITIONS_H
#define HYPER_MAKE_TRACTION_CONDITIONS_H

#include <hyper_traction_info.h>
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
  inline void make_traction_conditions(
                                  Body<MT> & body
                                  , typename MT::vector3_type const & load
                                  , ImplicitSurfaceFunctor<MT> const & phi
                                  , use_deformed_coordinates const & /*tag*/
                                  )
  {
    typedef typename MT::vector3_type    V;
    typedef typename MT::value_traits    VT;

    body.m_traction_conditions.clear();

    unsigned int const K = body.m_mesh.tetrahedron_size();

    for ( unsigned int idx = 0u; idx < K; ++idx)
    {
      mesh_array::Tetrahedron const & tetrahedron = body.m_mesh.tetrahedron(idx);

      unsigned int const i  = tetrahedron.i();
      unsigned int const j  = tetrahedron.j();
      unsigned int const k  = tetrahedron.k();
      unsigned int const m  = tetrahedron.m();

      V const pi = V::make( body.m_X(i), body.m_Y(i), body.m_Z(i) );
      V const pj = V::make( body.m_X(j), body.m_Y(j), body.m_Z(j) );
      V const pm = V::make( body.m_X(k), body.m_Y(k), body.m_Z(k) );
      V const pk = V::make( body.m_X(m), body.m_Y(m), body.m_Z(m) );

      bool const i_inside = phi(pi) <= VT::zero();
      bool const j_inside = phi(pj) <= VT::zero();
      bool const k_inside = phi(pk) <= VT::zero();
      bool const m_inside = phi(pm) <= VT::zero();

      if( i_inside && k_inside && j_inside )
        body.m_traction_conditions.push_back( make_traction_info<MT>( i, k, j, load ));

      if( i_inside && j_inside && m_inside )
        body.m_traction_conditions.push_back( make_traction_info<MT>( i, j, m, load ));

      if( j_inside && k_inside && m_inside )
        body.m_traction_conditions.push_back( make_traction_info<MT>( j, k, m, load ));

      if( i_inside && m_inside && k_inside )
        body.m_traction_conditions.push_back( make_traction_info<MT>( i, m, k, load ));

    }
  }

  template<typename MT>
  inline void make_traction_conditions(
                                  Body<MT> & body
                                  , typename MT::vector3_type const & load
                                  , ImplicitSurfaceFunctor<MT> const & phi
                                  , use_material_coordinates const & /*tag*/
                                  )
  {
    typedef typename MT::vector3_type    V;
    typedef typename MT::value_traits    VT;

    body.m_traction_conditions.clear();

    unsigned int const K = body.m_mesh.tetrahedron_size();

    for ( unsigned int idx = 0u; idx < K; ++idx)
    {
      mesh_array::Tetrahedron const & tetrahedron = body.m_mesh.tetrahedron(idx);

      unsigned int const i  = tetrahedron.i();
      unsigned int const j  = tetrahedron.j();
      unsigned int const k  = tetrahedron.k();
      unsigned int const m  = tetrahedron.m();

      V const pi = V::make( body.m_X0(i), body.m_Y0(i), body.m_Z0(i) );
      V const pj = V::make( body.m_X0(j), body.m_Y0(j), body.m_Z0(j) );
      V const pk = V::make( body.m_X0(k), body.m_Y0(k), body.m_Z0(k) );
      V const pm = V::make( body.m_X0(m), body.m_Y0(m), body.m_Z0(m) );

      bool const i_inside = phi(pi) <= VT::zero();
      bool const j_inside = phi(pj) <= VT::zero();
      bool const k_inside = phi(pk) <= VT::zero();
      bool const m_inside = phi(pm) <= VT::zero();

      if( i_inside && k_inside && j_inside )
        body.m_traction_conditions.push_back( make_traction_info<MT>( i, k, j, load ));

      if( i_inside && j_inside && m_inside )
        body.m_traction_conditions.push_back( make_traction_info<MT>( i, j, m, load ));

      if( j_inside && k_inside && m_inside )
        body.m_traction_conditions.push_back( make_traction_info<MT>( j, k, m, load ));

      if( i_inside && m_inside && k_inside )
        body.m_traction_conditions.push_back( make_traction_info<MT>( i, m, k, load ));
      
    }

  }


  template<typename MT>
  inline void make_traction_conditions(
                                  Body<MT> & body
                                 , typename MT::vector3_type const & load
                                 , ImplicitSurfaceFunctor<MT> const & phi
                                 )
  {
    make_traction_conditions(body, load, phi, use_material_coordinates() );
  }


} // namespace hyper

// HYPER_MAKE_TRACTION_CONDITIONS_H
#endif 
