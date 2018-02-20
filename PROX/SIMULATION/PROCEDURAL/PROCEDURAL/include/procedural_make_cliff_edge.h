#ifndef PROCEDURAL_MAKE_CLIFF_EDGE_H
#define PROCEDURAL_MAKE_CLIFF_EDGE_H

#include <procedural_factory.h>

#include <mesh_array.h>

#include <tiny_math_types.h>

#include <vector>

namespace procedural
{

  template<typename MT>
  inline void make_cliff_edge(
                              content::API *  engine
                              , typename MT::vector3_type const & position
                              , typename MT::quaternion_type const & orientation
                              , typename MT::real_type const & scene_size
                              , MaterialInfo<typename MT::real_type> const & mat_info
                              , mesh_array::TetGenSettings tetset = mesh_array::tetgen_default_settings()
                              )
  {
    using std::floor;
    using std::ceil;

    typedef typename MT::real_type       T;
    typedef typename MT::vector3_type    V;
    typedef typename MT::quaternion_type Q;
    //typedef typename MT::value_traits    VT;

    T      const stone_density	 = get_material_density<MT>(mat_info, "Stone");
    size_t const mid		      	 = get_material_id<MT>(mat_info, "Stone");

    T const H = 0.2*scene_size;
    T const B = 1.5*H;

    std::vector< V > base_vertices;

    //             h6 +-------------------------+ h7
    //
    //
    //   h4 +-------------------------+ h5
    //
    //
    //             h2 +-------------------------+ h3
    //
    //
    //   h0 +-------------------------+ h1

    base_vertices.resize( 8u );
    base_vertices[0] = V::make( -B, -H,  B);
    base_vertices[1] = V::make(  B, -H,  B);
    base_vertices[2] = V::make( -B, -H, -B);
    base_vertices[3] = V::make(  B, -H, -B);
    base_vertices[4] = V::make( -H,  H,  H);
    base_vertices[5] = V::make(  H,  H,  H);
    base_vertices[6] = V::make( -H,  H, -H);
    base_vertices[7] = V::make(  H,  H, -H);

    GeometryHandle<MT> const base = create_geometry_handle_convex<MT>(engine, base_vertices);

    std::vector< V > top_vertices;
    top_vertices.resize( 8u );
    top_vertices[0] = V::make( -H,  H,  H);
    top_vertices[1] = V::make(  H,  H,  H);
    top_vertices[2] = V::make( -H,  H, -H);
    top_vertices[3] = V::make(  H,  H, -H);
    top_vertices[4] = V::make( -H, -H,  H);
    top_vertices[5] = V::make(  H, -H,  H);
    top_vertices[6] = V::make( -H, -H, -H);
    top_vertices[7] = V::make(  H, -H, -H);


    GeometryHandle<MT> const top = create_geometry_handle_convex<MT>(engine, top_vertices);

    {
      V const T_b2m = base.Tb2m();
      Q const Q_b2m = base.Qb2m();

      V const T_m2l = V::make( 0.0, H, 0.0 );
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
                                               , base
                                               , mid
                                               , stone_density
                                               );
      engine->set_rigid_body_fixed(rid, true);
    }

    {
      V const T_b2m = top.Tb2m();
      Q const Q_b2m = top.Qb2m();

      V const T_m2l = V::make( 0.0, 3.0*H, 0.0 );
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
                                               , top
                                               , mid
                                               , stone_density
                                               );

      engine->set_rigid_body_velocity(rid, 0.0, 0.0, 0.0);
    }

  }

} //namespace procedural

// PROCEDURAL_MAKE_CLIFF_EDGE_H
#endif
