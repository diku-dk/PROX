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

    T const stone_density = get_material_density_eigen<T>(mat_info, "Stone");
    size_t const mid = get_material_id_eigen<T>(mat_info, "Stone");

    T const H = 0.2*scene_size;
    T const B = 1.5*H;

    std::vector< EigenVector3<T> > base_vertices;

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
    base_vertices[0] = EigenVector3<T>( -B, -H,  B);
    base_vertices[1] = EigenVector3<T>(  B, -H,  B);
    base_vertices[2] = EigenVector3<T>( -B, -H, -B);
    base_vertices[3] = EigenVector3<T>(  B, -H, -B);
    base_vertices[4] = EigenVector3<T>( -H,  H,  H);
    base_vertices[5] = EigenVector3<T>(  H,  H,  H);
    base_vertices[6] = EigenVector3<T>( -H,  H, -H);
    base_vertices[7] = EigenVector3<T>(  H,  H, -H);

    GeometryHandleEigen<T> const base = create_geometry_handle_convex<T>(engine, base_vertices);

    std::vector<EigenVector3<T>> top_vertices;
    top_vertices.resize( 8u );
    top_vertices[0] = EigenVector3<T>( -H,  H,  H);
    top_vertices[1] = EigenVector3<T>(  H,  H,  H);
    top_vertices[2] = EigenVector3<T>( -H,  H, -H);
    top_vertices[3] = EigenVector3<T>(  H,  H, -H);
    top_vertices[4] = EigenVector3<T>( -H, -H,  H);
    top_vertices[5] = EigenVector3<T>(  H, -H,  H);
    top_vertices[6] = EigenVector3<T>( -H, -H, -H);
    top_vertices[7] = EigenVector3<T>(  H, -H, -H);


    GeometryHandleEigen<T> const top = create_geometry_handle_convex<T>(engine, top_vertices);

    {
        const EigenVector3<T> T_b2m = (base.Tb2m());
        const EigenQuaternion<T> Q_b2m = (base.Qb2m());

        const EigenVector3<T> T_m2l = EigenVector3<T>(0.0, H, 0.0);
        const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

        const EigenVector3<T> T_l2w = toEigen(position);
        const EigenQuaternion<T> Q_l2w = toEigen(orientation);

        EigenVector3<T> T_b2w;
        EigenQuaternion<T> Q_b2w;

        compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w, Q_l2w, T_b2w, Q_b2w);

        size_t const rid = create_rigid_body<T>(engine, T_b2w, Q_b2w, base, mid, stone_density);
        engine->set_rigid_body_fixed(rid, true);
    }

    {
        const EigenVector3<T> T_b2m = (top.Tb2m());
        const EigenQuaternion<T> Q_b2m = (top.Qb2m());

        const EigenVector3<T> T_m2l = EigenVector3<T>(0.0, 3.0 * H, 0.0);
        const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

        const EigenVector3<T> T_l2w = toEigen(position);
        const EigenQuaternion<T> Q_l2w = toEigen(orientation);

        EigenVector3<T> T_b2w;
        EigenQuaternion<T> Q_b2w;

        compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w, Q_l2w, T_b2w, Q_b2w);

        size_t const rid = create_rigid_body<T>(engine, T_b2w, Q_b2w, top, mid, stone_density);

        engine->set_rigid_body_velocity(rid, 0.0, 0.0, 0.0);
    }

  }

} //namespace procedural

// PROCEDURAL_MAKE_CLIFF_EDGE_H
#endif
