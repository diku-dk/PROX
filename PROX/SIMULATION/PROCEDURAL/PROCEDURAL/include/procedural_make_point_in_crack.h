#ifndef PROCEDURAL_POINT_IN_CRACK_H
#define PROCEDURAL_POINT_IN_CRACK_H

#include <procedural_factory.h>

#include <mesh_array.h>

#include <tiny_math_types.h>

#include <vector>

namespace procedural
{

template <typename T>
inline void make_point_in_crack(content::API* engine, const EigenVector3<T>& position,
                                const EigenQuaternion<T>& orientation, const T& scene_size,
                                MaterialInfo<T> const& mat_info, std::string const& obj_path,
                                bool const& use_spike_in_hole = true, bool const& use_spike_in_crack = true,
                                bool const& use_wedge_in_crack = true,
                                mesh_array::TetGenSettings tetset = mesh_array::tetgen_default_settings())
{
    using std::floor;
    using std::ceil;

    T const stone_density = get_material_density_eigen<T>(mat_info, "Stone");
    size_t const mid = get_material_id_eigen<T>(mat_info, "Stone");

    T const H = 6.0 * scene_size / 50.0;
    T const B = 2.0*scene_size / 50.0;

    std::vector<EigenVector3<T>> spike_vertices;

    spike_vertices.resize( 5u );
    spike_vertices[0] = EigenVector3<T>(0.0,H,0.0);
    spike_vertices[1] = EigenVector3<T>( B,0.0, B);
    spike_vertices[2] = EigenVector3<T>(-B,0.0, B);
    spike_vertices[3] = EigenVector3<T>( B,0.0,-B);
    spike_vertices[4] = EigenVector3<T>(-B,0.0,-B);

    GeometryHandleEigen<T> const spike = create_geometry_handle_convex<T>(engine, spike_vertices);


    std::vector<EigenVector3<T>> wedge_vertices;
    wedge_vertices.resize( 6u );

    wedge_vertices[0] = EigenVector3<T>( 0.0,   H, -2.0*B );
    wedge_vertices[1] = EigenVector3<T>( 0.0,   H,  2.0*B );
    wedge_vertices[2] = EigenVector3<T>(   B, 0.0,  2.0*B );
    wedge_vertices[3] = EigenVector3<T>(  -B, 0.0,  2.0*B );
    wedge_vertices[4] = EigenVector3<T>(   B, 0.0, -2.0*B );
    wedge_vertices[5] = EigenVector3<T>(  -B, 0.0, -2.0*B );

    GeometryHandleEigen<T> const wedge = create_geometry_handle_convex<T>(engine, wedge_vertices);

    GeometryHandleEigen<T> const crack = create_geometry_handle_obj<T>(engine, obj_path + "crack.obj", 1.0, 1.0, 4.0,
                                                                       mesh_array::tetgen_cdt_settings());

    GeometryHandleEigen<T> const hole = create_geometry_handle_obj<T>(engine, obj_path + "hole.obj", 1.0, 1.0, 1.0,
                                                                      mesh_array::tetgen_cdt_settings());

    if(use_spike_in_hole)
    {
      {
          const EigenVector3<T> T_b2m = (spike.Tb2m());
          const EigenQuaternion<T> Q_b2m = (spike.Qb2m());

          const EigenVector3<T> T_m2l = EigenVector3<T>(0, H, 0.0);
          const EigenQuaternion<T> Q_m2l = Rotatex(std::numbers::pi_v<T>);

          const EigenVector3<T> T_l2w = (position);
          const EigenQuaternion<T> Q_l2w = (orientation);

          EigenVector3<T> T_b2w;
          EigenQuaternion<T> Q_b2w;

          compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w, Q_l2w, T_b2w, Q_b2w);

          size_t const rid = create_rigid_body<T>(engine, T_b2w, Q_b2w, spike, mid, stone_density);

          engine->set_rigid_body_velocity(rid, B, 0.0, 0.0);
      }
      {
          const EigenVector3<T> T_b2m = hole.Tb2m();
          const EigenQuaternion<T> Q_b2m = hole.Qb2m();

          const EigenVector3<T> T_m2l = EigenVector3<T>(0, 0.0, 0.0);
          const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

          const EigenVector3<T> T_l2w = (position);
          const EigenQuaternion<T> Q_l2w = (orientation);

          EigenVector3<T> T_b2w;
          EigenQuaternion<T> Q_b2w;

          compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w, Q_l2w, T_b2w, Q_b2w);

          size_t const rid = create_rigid_body<T>(engine, T_b2w, Q_b2w, hole, mid, stone_density);

          engine->set_rigid_body_fixed(rid, true);
      }
    }
    if(use_spike_in_crack)
    {
      {
          const EigenVector3<T> T_b2m = (spike.Tb2m());
          const EigenQuaternion<T> Q_b2m = (spike.Qb2m());

          const EigenVector3<T> T_m2l = EigenVector3<T>(H, H, 0.0);
          const EigenQuaternion<T> Q_m2l = Rotatex(std::numbers::pi_v<T>);

          const EigenVector3<T> T_l2w = (position);
          const EigenQuaternion<T> Q_l2w = (orientation);

          EigenVector3<T> T_b2w;
          EigenQuaternion<T> Q_b2w;

          compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w, Q_l2w, T_b2w, Q_b2w);

          size_t const rid = create_rigid_body<T>(engine, T_b2w, Q_b2w, spike, mid, stone_density);

          engine->set_rigid_body_velocity(rid, B, 0.0, 0.0);
      }
      {
          const EigenVector3<T> T_b2m = crack.Tb2m();
          const EigenQuaternion<T> Q_b2m = crack.Qb2m();

          const EigenVector3<T> T_m2l = EigenVector3<T>(H, 0.0, 0.0);
          const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

          const EigenVector3<T> T_l2w = (position);
          const EigenQuaternion<T> Q_l2w = (orientation);

          EigenVector3<T> T_b2w;
          EigenQuaternion<T> Q_b2w;

          compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w, Q_l2w, T_b2w, Q_b2w);

          size_t const rid = create_rigid_body<T>(engine, T_b2w, Q_b2w, crack, mid, stone_density);

          engine->set_rigid_body_fixed(rid, true);
      }
    }
    if(use_wedge_in_crack)
    {
      {
          const EigenVector3<T> T_b2m = (wedge.Tb2m());
          const EigenQuaternion<T> Q_b2m = (wedge.Qb2m());

          const EigenVector3<T> T_m2l = EigenVector3<T>(-H, H, 0.0);
          const EigenQuaternion<T> Q_m2l = Rotatex(std::numbers::pi_v<T>);

          const EigenVector3<T> T_l2w = (position);
          const EigenQuaternion<T> Q_l2w = (orientation);

          EigenVector3<T> T_b2w;
          EigenQuaternion<T> Q_b2w;

          compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w, Q_l2w, T_b2w, Q_b2w);

          size_t const rid = create_rigid_body<T>(engine, T_b2w, Q_b2w, wedge, mid, stone_density);

          engine->set_rigid_body_velocity(rid, B, 0.0, 0.0);
      }
      {
          const EigenVector3<T> T_b2m = crack.Tb2m();
          const EigenQuaternion<T> Q_b2m = crack.Qb2m();

          const EigenVector3<T> T_m2l = EigenVector3<T>(-H, 0.0, 0.0);
          const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

          const EigenVector3<T> T_l2w = (position);
          const EigenQuaternion<T> Q_l2w = (orientation);

          EigenVector3<T> T_b2w;
          EigenQuaternion<T> Q_b2w;

          compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w, Q_l2w, T_b2w, Q_b2w);

          size_t const rid = create_rigid_body<T>(engine, T_b2w, Q_b2w, crack, mid, stone_density);

          engine->set_rigid_body_fixed(rid, true);
      }
    }
}

} //namespace procedural

// PROCEDURAL_POINT_IN_CRACK_H
#endif
