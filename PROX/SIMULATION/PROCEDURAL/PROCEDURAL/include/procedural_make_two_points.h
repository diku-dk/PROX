#ifndef PROCEDURAL_MAKE_TWO_POINTS_H
#define PROCEDURAL_MAKE_TWO_POINTS_H

#include <procedural_factory.h>

#include <mesh_array.h>

#include <tiny_math_types.h>

#include <vector>

namespace procedural
{

  template<typename MT>
  inline void make_two_points(
                              content::API *  engine
                              , typename MT::vector3_type const & position
                              , typename MT::quaternion_type const & orientation
                              , typename MT::real_type const & scene_size
                              , MaterialInfo<typename MT::real_type> const & mat_info
                              , bool const & use_spike = true
                              , bool const & use_wedge = true
                              , bool const & use_spike_and_wedge = true
                              , mesh_array::TetGenSettings tetset = mesh_array::tetgen_default_settings()
                              )
  {
    using std::floor;
    using std::ceil;

    typedef typename MT::real_type       T;
    typedef typename MT::vector3_type    V;
    typedef typename MT::quaternion_type Q;
    typedef typename MT::value_traits    VT;

    T      const stone_density	 = get_material_density<MT>(mat_info, "Stone");
    size_t const mid		      	 = get_material_id<MT>(mat_info, "Stone");


    T const H = 6.0*scene_size / 50.0;
    T const B = 2.0*scene_size / 50.0;

    if(use_spike)
    {
      std::vector< V > spike_vertices;
      spike_vertices.resize( 5u );
      spike_vertices[0] = V::make(0.0,H,0.0);
      spike_vertices[1] = V::make( B,0.0, B);
      spike_vertices[2] = V::make(-B,0.0, B);
      spike_vertices[3] = V::make( B,0.0,-B);
      spike_vertices[4] = V::make(-B,0.0,-B);

      GeometryHandle<MT> const spike = create_geometry_handle_convex<MT>(engine, spike_vertices);

      {
        V const T_b2m = spike.Tb2m();
        Q const Q_b2m = spike.Qb2m();

        V const T_m2l = V::make( H, H, 0.0 );
        Q const Q_m2l = Q::Rx( VT::pi()  );

        V const T_l2w = position;
        Q const Q_l2w = orientation;

        V T_b2w;
        Q Q_b2w;

        compute_body_to_world_transform<MT>(
                                            T_b2m
                                            , Q_b2m
                                            , T_m2l
                                            , Q_m2l
                                            , T_l2w
                                            , Q_l2w
                                            , T_b2w
                                            , Q_b2w
                                            );

        size_t const rid = create_rigid_body<MT>(
                                                 engine
                                                 , T_b2w
                                                 , Q_b2w
                                                 , spike
                                                 , mid
                                                 , stone_density
                                                 );

        engine->set_rigid_body_velocity(rid, -H/2.0, -H*2.0, 0.0);
      }
      {
        V const T_b2m = spike.Tb2m();
        Q const Q_b2m = spike.Qb2m();

        V const T_m2l = V::make( H, -H, 0.0 );
        Q const Q_m2l = Q::identity();

        V const T_l2w = position;
        Q const Q_l2w = orientation;

        V T_b2w;
        Q Q_b2w;

        compute_body_to_world_transform<MT>(
                                            T_b2m
                                            , Q_b2m
                                            , T_m2l
                                            , Q_m2l
                                            , T_l2w
                                            , Q_l2w
                                            , T_b2w
                                            , Q_b2w
                                            );

        size_t const rid = create_rigid_body<MT>(
                                                 engine
                                                 , T_b2w
                                                 , Q_b2w
                                                 , spike
                                                 , mid
                                                 , stone_density
                                                 );

        engine->set_rigid_body_fixed(rid, true);
      }
    }

    if(use_spike_and_wedge)
    {
      std::vector< V > spike_vertices;
      spike_vertices.resize( 5u );
      spike_vertices[0] = V::make(0.0,H,0.0);
      spike_vertices[1] = V::make( B,0.0, B);
      spike_vertices[2] = V::make(-B,0.0, B);
      spike_vertices[3] = V::make( B,0.0,-B);
      spike_vertices[4] = V::make(-B,0.0,-B);

      GeometryHandle<MT> const spike = create_geometry_handle_convex<MT>(engine, spike_vertices);

      {
        V const T_b2m = spike.Tb2m();
        Q const Q_b2m = spike.Qb2m();

        V const T_m2l = V::make( 0, H, 0.0 );
        Q const Q_m2l = Q::Rx( VT::pi()  );

        V const T_l2w = position;
        Q const Q_l2w = orientation;

        V T_b2w;
        Q Q_b2w;

        compute_body_to_world_transform<MT>(
                                            T_b2m
                                            , Q_b2m
                                            , T_m2l
                                            , Q_m2l
                                            , T_l2w
                                            , Q_l2w
                                            , T_b2w
                                            , Q_b2w
                                            );

        size_t const rid = create_rigid_body<MT>(
                                                 engine
                                                 , T_b2w
                                                 , Q_b2w
                                                 , spike
                                                 , mid
                                                 , stone_density
                                                 );

        engine->set_rigid_body_velocity(rid, -H/2.0, -H*2.0, 0.0);
      }

      std::vector< V > wedge_vertices;
      wedge_vertices.resize( 6u );

      wedge_vertices[0] = V::make( 0.0,   H, -2.0*B );
      wedge_vertices[1] = V::make( 0.0,   H,  2.0*B );
      wedge_vertices[2] = V::make(   B, 0.0,  2.0*B );
      wedge_vertices[3] = V::make(  -B, 0.0,  2.0*B );
      wedge_vertices[4] = V::make(   B, 0.0, -2.0*B );
      wedge_vertices[5] = V::make(  -B, 0.0, -2.0*B );

      GeometryHandle<MT> const wedge = create_geometry_handle_convex<MT>(engine, wedge_vertices);

      {
        V const T_b2m = wedge.Tb2m();
        Q const Q_b2m = wedge.Qb2m();

        V const T_m2l = V::make( 0, -H, 0.0 );
        Q const Q_m2l = Q::identity();

        V const T_l2w = position;
        Q const Q_l2w = orientation;

        V T_b2w;
        Q Q_b2w;

        compute_body_to_world_transform<MT>(
                                            T_b2m
                                            , Q_b2m
                                            , T_m2l
                                            , Q_m2l
                                            , T_l2w
                                            , Q_l2w
                                            , T_b2w
                                            , Q_b2w
                                            );

        size_t const rid = create_rigid_body<MT>(
                                                 engine
                                                 , T_b2w
                                                 , Q_b2w
                                                 , wedge
                                                 , mid
                                                 , stone_density
                                                 );
        
        engine->set_rigid_body_fixed(rid, true);
      }
      
    }



    if(use_wedge)
    {
      std::vector< V > wedge_vertices;
      wedge_vertices.resize( 6u );

      wedge_vertices[0] = V::make( 0.0,   H, -2.0*B );
      wedge_vertices[1] = V::make( 0.0,   H,  2.0*B );
      wedge_vertices[2] = V::make(   B, 0.0,  2.0*B );
      wedge_vertices[3] = V::make(  -B, 0.0,  2.0*B );
      wedge_vertices[4] = V::make(   B, 0.0, -2.0*B );
      wedge_vertices[5] = V::make(  -B, 0.0, -2.0*B );

      GeometryHandle<MT> const wedge = create_geometry_handle_convex<MT>(engine, wedge_vertices);

      {
        V const T_b2m = wedge.Tb2m();
        Q const Q_b2m = wedge.Qb2m();

        V const T_m2l = V::make( -H, H, 0.0 );
        Q const Q_m2l = Q::Rx( VT::pi()  );

        V const T_l2w = position;
        Q const Q_l2w = orientation;

        V T_b2w;
        Q Q_b2w;

        compute_body_to_world_transform<MT>(
                                            T_b2m
                                            , Q_b2m
                                            , T_m2l
                                            , Q_m2l
                                            , T_l2w
                                            , Q_l2w
                                            , T_b2w
                                            , Q_b2w
                                            );

        size_t const rid = create_rigid_body<MT>(
                                                 engine
                                                 , T_b2w
                                                 , Q_b2w
                                                 , wedge
                                                 , mid
                                                 , stone_density
                                                 );

        engine->set_rigid_body_velocity(rid, H/2.0, -H*2.0, 0.0);
      }
      {
        V const T_b2m = wedge.Tb2m();
        Q const Q_b2m = wedge.Qb2m();

        V const T_m2l = V::make( -H, -H, 0.0 );
        Q const Q_m2l = Q::identity();

        V const T_l2w = position;
        Q const Q_l2w = orientation;

        V T_b2w;
        Q Q_b2w;

        compute_body_to_world_transform<MT>(
                                            T_b2m
                                            , Q_b2m
                                            , T_m2l
                                            , Q_m2l
                                            , T_l2w
                                            , Q_l2w
                                            , T_b2w
                                            , Q_b2w
                                            );

        size_t const rid = create_rigid_body<MT>(
                                                 engine
                                                 , T_b2w
                                                 , Q_b2w
                                                 , wedge
                                                 , mid
                                                 , stone_density
                                                 );
        
        engine->set_rigid_body_fixed(rid, true);
      }
      
    }
    
  }
  
} //namespace procedural

// PROCEDURAL_MAKE_TWO_POINTS_H
#endif
