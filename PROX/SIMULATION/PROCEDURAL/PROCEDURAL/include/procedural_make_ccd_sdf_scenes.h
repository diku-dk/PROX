#ifndef PROCEDURAL_MAKE_CCD_SDF_SCENES_H
#define PROCEDURAL_MAKE_CCD_SDF_SCENES_H

#include "procedural.h"
#include <procedural_factory.h>

#include <mesh_array.h>

#include <eigenhelperall.h>

#include <vector>

namespace procedural
{
template <typename T>
void make_sdf_thin_object_drop(content::API* engine,
                               const EigenVector3<T>& position,
                               const EigenQuaternion<T>& orientation,
                               const T& scene_size,
                               MaterialInfo<T> const& mat_info,
                               mesh_array::TetGenSettings tetset
                               = mesh_array::tetgen_default_settings())
{
    using std::ceil;
    using std::floor;

    T const stone_density = get_material_density_eigen<T>(mat_info, "Stone");
    size_t const mid = get_material_id_eigen<T>(mat_info, "Stone");

    T const H = 6.0 * scene_size / 50.0;
    T const B = 2.0 * scene_size / 50.0;

    if (true)
    {
        std::vector<EigenVector3<T>> spike_vertices;
        spike_vertices.resize(5u);
        spike_vertices[0] = EigenVector3<T>(0.0, H, 0.0);
        spike_vertices[1] = EigenVector3<T>(B, 0.0, B);
        spike_vertices[2] = EigenVector3<T>(-B, 0.0, B);
        spike_vertices[3] = EigenVector3<T>(B, 0.0, -B);
        spike_vertices[4] = EigenVector3<T>(-B, 0.0, -B);

        GeometryHandleEigen<T> const spike
            = create_geometry_handle_convex<T>(engine, spike_vertices);

        const EigenVector3<T> T_b2m = (spike.Tb2m());
        const EigenQuaternion<T> Q_b2m = (spike.Qb2m());

        const EigenVector3<T> T_m2l = EigenVector3<T>(-H * 0.5, H, H * 0.15);
        const EigenQuaternion<T> Q_m2l
            = Rotatex(std::numbers::pi_v<T> * T(1.5));

        const EigenVector3<T> T_l2w = (position);
        const EigenQuaternion<T> Q_l2w = (orientation);

        EigenVector3<T> T_b2w;
        EigenQuaternion<T> Q_b2w;

        compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w,
                                           Q_l2w, T_b2w, Q_b2w);

        size_t const rid = create_rigid_body<T>(engine, T_b2w, Q_b2w, spike,
                                                mid, stone_density);

        engine->set_rigid_body_velocity(rid, H * 0.25 * 4.0, 0.0, 0.0);
    }

    if (false)
    {
        std::vector<EigenVector3<T>> wedge_vertices;
        wedge_vertices.resize(6u);

        wedge_vertices[0] = EigenVector3<T>(0.0, H, -2.0 * B);
        wedge_vertices[1] = EigenVector3<T>(0.0, H, 2.0 * B);
        wedge_vertices[2] = EigenVector3<T>(B, 0.0, 2.0 * B);
        wedge_vertices[3] = EigenVector3<T>(-B, 0.0, 2.0 * B);
        wedge_vertices[4] = EigenVector3<T>(B, 0.0, -2.0 * B);
        wedge_vertices[5] = EigenVector3<T>(-B, 0.0, -2.0 * B);

        GeometryHandleEigen<T> const wedge
            = create_geometry_handle_convex<T>(engine, wedge_vertices);

        const EigenVector3<T> T_b2m = (wedge.Tb2m());
        const EigenQuaternion<T> Q_b2m = (wedge.Qb2m());

        const EigenVector3<T> T_m2l = EigenVector3<T>(0.0, H * 3.0, 0.0);
        const EigenQuaternion<T> Q_m2l = Rotatex(std::numbers::pi_v<T>);

        const EigenVector3<T> T_l2w = (position);
        const EigenQuaternion<T> Q_l2w = (orientation);

        EigenVector3<T> T_b2w;
        EigenQuaternion<T> Q_b2w;

        compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w,
                                           Q_l2w, T_b2w, Q_b2w);

        size_t const rid = create_rigid_body<T>(engine, T_b2w, Q_b2w, wedge,
                                                mid, stone_density);

        engine->set_rigid_body_velocity(rid, 0.0, 0.0, 0.0);
    }

    if (true)
    {
        std::vector<EigenVector3<T>> spike_vertices;
        spike_vertices.resize(5u);
        spike_vertices[0] = EigenVector3<T>(0.0, H, 0.0);
        spike_vertices[1] = EigenVector3<T>(B, 0.0, B);
        spike_vertices[2] = EigenVector3<T>(-B, 0.0, B);
        spike_vertices[3] = EigenVector3<T>(B, 0.0, -B);
        spike_vertices[4] = EigenVector3<T>(-B, 0.0, -B);

        GeometryHandleEigen<T> const spike
            = create_geometry_handle_convex<T>(engine, spike_vertices);

        const EigenVector3<T> T_b2m = (spike.Tb2m());
        const EigenQuaternion<T> Q_b2m = (spike.Qb2m());

        const EigenVector3<T> T_m2l = EigenVector3<T>(0.0, H * 3.0, 0.0);
        const EigenQuaternion<T> Q_m2l = Rotatex(std::numbers::pi_v<T>);

        const EigenVector3<T> T_l2w = (position);
        const EigenQuaternion<T> Q_l2w = (orientation);

        EigenVector3<T> T_b2w;
        EigenQuaternion<T> Q_b2w;

        compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w,
                                           Q_l2w, T_b2w, Q_b2w);

        size_t const rid = create_rigid_body<T>(engine, T_b2w, Q_b2w, spike,
                                                mid, stone_density);

        engine->set_rigid_body_velocity(rid, H * 0.25, 0.0, 0.0);
    }

    if (true)
    {
        std::vector<EigenVector3<T>> spike_vertices;
        spike_vertices.resize(5u);
        spike_vertices[0] = EigenVector3<T>(0.0, H, 0.0);
        spike_vertices[1] = EigenVector3<T>(B, 0.0, B);
        spike_vertices[2] = EigenVector3<T>(-B, 0.0, B);
        spike_vertices[3] = EigenVector3<T>(B, 0.0, -B);
        spike_vertices[4] = EigenVector3<T>(-B, 0.0, -B);

        GeometryHandleEigen<T> const spike
            = create_geometry_handle_convex<T>(engine, spike_vertices);

        const EigenVector3<T> T_b2m = (spike.Tb2m());
        const EigenQuaternion<T> Q_b2m = (spike.Qb2m());

        const EigenVector3<T> T_m2l = EigenVector3<T>(-H * 2.0, H, 0.0);
        const EigenQuaternion<T> Q_m2l = Rotatex(std::numbers::pi_v<T>);

        const EigenVector3<T> T_l2w = (position);
        const EigenQuaternion<T> Q_l2w = (orientation);

        EigenVector3<T> T_b2w;
        EigenQuaternion<T> Q_b2w;

        compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w,
                                           Q_l2w, T_b2w, Q_b2w);

        size_t const rid = create_rigid_body<T>(engine, T_b2w, Q_b2w, spike,
                                                mid, stone_density);

        engine->set_rigid_body_velocity(rid, H * 0.15, 0.0, 0.0);
    }

    if (true)
    {
        std::vector<EigenVector3<T>> spike_vertices;
        spike_vertices.resize(5u);
        spike_vertices[0] = EigenVector3<T>(0.0, H, 0.0);
        spike_vertices[1] = EigenVector3<T>(B, 0.0, B);
        spike_vertices[2] = EigenVector3<T>(-B, 0.0, B);
        spike_vertices[3] = EigenVector3<T>(B, 0.0, -B);
        spike_vertices[4] = EigenVector3<T>(-B, 0.0, -B);

        GeometryHandleEigen<T> const spike
            = create_geometry_handle_convex<T>(engine, spike_vertices);

        const EigenVector3<T> T_b2m = (spike.Tb2m());
        const EigenQuaternion<T> Q_b2m = (spike.Qb2m());

        const EigenVector3<T> T_m2l = EigenVector3<T>(-H * 2.0, H * 3.0, 0.75);
        const EigenQuaternion<T> Q_m2l = Rotatex(std::numbers::pi_v<T>);

        const EigenVector3<T> T_l2w = (position);
        const EigenQuaternion<T> Q_l2w = (orientation);

        EigenVector3<T> T_b2w;
        EigenQuaternion<T> Q_b2w;

        compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w,
                                           Q_l2w, T_b2w, Q_b2w);

        size_t const rid = create_rigid_body<T>(engine, T_b2w, Q_b2w, spike,
                                                mid, stone_density);

        engine->set_rigid_body_velocity(rid, H * 0.0, 0.0, 0.0);
    }
}

template <typename T>
void make_sdf_thin_thin_object_drop(content::API* engine,
                                    const EigenVector3<T>& position,
                                    const EigenQuaternion<T>& orientation,
                                    const T& scene_size,
                                    MaterialInfo<T> const& mat_info,
                                    mesh_array::TetGenSettings tetset
                                    = mesh_array::tetgen_default_settings())
{
    using std::ceil;
    using std::floor;

    T const stone_density = get_material_density_eigen<T>(mat_info, "Stone");
    size_t const mid = get_material_id_eigen<T>(mat_info, "Stone");

    T const H = 6.0 * scene_size / 50.0;
    T const B = 2.0 * scene_size / 50.0;

    if (true)
    {
        std::vector<EigenVector3<T>> spike_vertices;
        spike_vertices.resize(5u);
        spike_vertices[0] = EigenVector3<T>(0.0, H, 0.0);
        spike_vertices[1] = EigenVector3<T>(B, 0.0, B);
        spike_vertices[2] = EigenVector3<T>(-B, 0.0, B);
        spike_vertices[3] = EigenVector3<T>(B, 0.0, -B);
        spike_vertices[4] = EigenVector3<T>(-B, 0.0, -B);

        GeometryHandleEigen<T> const spike
            = create_geometry_handle_convex<T>(engine, spike_vertices);

        const EigenVector3<T> T_b2m = (spike.Tb2m());
        const EigenQuaternion<T> Q_b2m = (spike.Qb2m());

        const EigenVector3<T> T_m2l = EigenVector3<T>(0.0, 0.0, 0.0);
        const EigenQuaternion<T> Q_m2l
            = Rotatey(std::numbers::pi_v<T> * T(1.0));

        const EigenVector3<T> T_l2w = (position);
        const EigenQuaternion<T> Q_l2w = (orientation);

        EigenVector3<T> T_b2w;
        EigenQuaternion<T> Q_b2w;

        compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w,
                                           Q_l2w, T_b2w, Q_b2w);

        size_t const rid = create_rigid_body<T>(engine, T_b2w, Q_b2w, spike,
                                                mid, stone_density);

        engine->set_rigid_body_velocity(rid, 0.0, 0.0, 0.0);
    }
    if (true)
    {
        std::vector<EigenVector3<T>> spike_vertices;
        spike_vertices.resize(5u);
        spike_vertices[0] = EigenVector3<T>(0.0, H, 0.0);
        spike_vertices[1] = EigenVector3<T>(B, 0.0, B);
        spike_vertices[2] = EigenVector3<T>(-B, 0.0, B);
        spike_vertices[3] = EigenVector3<T>(B, 0.0, -B);
        spike_vertices[4] = EigenVector3<T>(-B, 0.0, -B);

        GeometryHandleEigen<T> const spike
            = create_geometry_handle_convex<T>(engine, spike_vertices);

        const EigenVector3<T> T_b2m = (spike.Tb2m());
        const EigenQuaternion<T> Q_b2m = (spike.Qb2m());

        const EigenVector3<T> T_m2l = EigenVector3<T>(0.0, H * 2.0, 0.0);
        const EigenQuaternion<T> Q_m2l
            = Rotatex(std::numbers::pi_v<T> * T(1.0));

        const EigenVector3<T> T_l2w = (position);
        const EigenQuaternion<T> Q_l2w = (orientation);

        EigenVector3<T> T_b2w;
        EigenQuaternion<T> Q_b2w;

        compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w,
                                           Q_l2w, T_b2w, Q_b2w);

        size_t const rid = create_rigid_body<T>(engine, T_b2w, Q_b2w, spike,
                                                mid, stone_density);

        engine->set_rigid_body_velocity(rid, 0.0, 0.0, 0.0);
    }
}

template <typename T>
void make_ccd_thin_thin_object_shoot(content::API* engine,
                                     const EigenVector3<T>& position,
                                     const EigenQuaternion<T>& orientation,
                                     const T& scene_size,
                                     MaterialInfo<T> const& mat_info,
                                     mesh_array::TetGenSettings tetset
                                     = mesh_array::tetgen_default_settings())
{

    T width = 5.0f;
    T height = 0.1f;
    T depth = 5.0f;
    T const H = 6.0 * scene_size / 50.0;

    T const stone_density = get_material_density_eigen<T>(mat_info, "Stone");
    size_t const mid = get_material_id_eigen<T>(mat_info, "Stone");

    GeometryHandleEigen<T> stone_handle
        = create_geometry_handle_box<T>(engine, width, height, depth);

    {
        T y = H;

        const EigenVector3<T> T_b2m = stone_handle.Tb2m();
        const EigenQuaternion<T> Q_b2m = stone_handle.Qb2m();

        const EigenVector3<T> T_m2l = EigenVector3<T>(0, y, 0);
        const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

        const EigenVector3<T> T_l2w = position;
        const EigenQuaternion<T> Q_l2w = orientation;

        EigenVector3<T> T_b2w;
        EigenQuaternion<T> Q_b2w;

        compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w,
                                           Q_l2w, T_b2w, Q_b2w);

        const size_t rid = create_rigid_body<T>(
            engine, T_b2w, Q_b2w, stone_handle, mid, stone_density);
        engine->set_rigid_body_velocity(rid, 0.0, 0.0, 0.0);
    }

    {
        make_cannonball<T>(engine, 2.0f, EigenVector3<T>(0.0, H * 2.0, 0.0),
                           EigenQuaternion<T>::Identity(),
                           EigenVector3<T>(0.0, -2.0, 0.0), mat_info);
    }
}

template <typename T>
void make_ccd_cannonball_fast_shoot(content::API* engine,
                                    const EigenVector3<T>& position,
                                    const EigenQuaternion<T>& orientation,
                                    const T& scene_size,
                                    MaterialInfo<T> const& mat_info,
                                    mesh_array::TetGenSettings tetset
                                    = mesh_array::tetgen_default_settings())
{
    T width = 5.0f;
    T height = 0.1f;
    T depth = 5.0f;
    T const H = 6.0 * scene_size / 50.0;

    T const stone_density = get_material_density_eigen<T>(mat_info, "Stone");
    size_t const mid = get_material_id_eigen<T>(mat_info, "Stone");

    GeometryHandleEigen<T> stone_handle
        = create_geometry_handle_box<T>(engine, width, height, depth);

    {
        T y = H;

        const EigenVector3<T> T_b2m = stone_handle.Tb2m();
        const EigenQuaternion<T> Q_b2m = stone_handle.Qb2m();

        const EigenVector3<T> T_m2l = EigenVector3<T>(0, y, 0);
        const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

        const EigenVector3<T> T_l2w = position;
        const EigenQuaternion<T> Q_l2w = orientation;

        EigenVector3<T> T_b2w;
        EigenQuaternion<T> Q_b2w;

        compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w,
                                           Q_l2w, T_b2w, Q_b2w);

        const size_t rid = create_rigid_body<T>(
            engine, T_b2w, Q_b2w, stone_handle, mid, stone_density);
        engine->set_rigid_body_velocity(rid, 0.0, 0.0, 0.0);
    }

    {
        make_cannonball<T>(engine, 0.5f, EigenVector3<T>(0.0, H * 2.0, 0.0),
                           EigenQuaternion<T>::Identity(),
                           EigenVector3<T>(0.0, -200.0, 0.0), mat_info);
    }
}

template <typename T>
void make_ccd_drop_fast_thin(content::API* engine,
                             const EigenVector3<T>& position,
                             const EigenQuaternion<T>& orientation,
                             const T& scene_size,
                             MaterialInfo<T> const& mat_info,
                             mesh_array::TetGenSettings tetset
                             = mesh_array::tetgen_default_settings())
{
    T width = 0.1f;
    T height = 5.0f;
    T depth = 0.1f;
    T const H = 6.0 * scene_size / 50.0;

    T const stone_density = get_material_density_eigen<T>(mat_info, "Stone");
    size_t const mid = get_material_id_eigen<T>(mat_info, "Stone");

    GeometryHandleEigen<T> stone_handle
        = create_geometry_handle_box<T>(engine, width, height, depth);

    {
        T y = H;

        const EigenVector3<T> T_b2m = stone_handle.Tb2m();
        const EigenQuaternion<T> Q_b2m = stone_handle.Qb2m();

        const EigenVector3<T> T_m2l = EigenVector3<T>(0, y * 3.0, 0);
        const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

        const EigenVector3<T> T_l2w = position;
        const EigenQuaternion<T> Q_l2w = orientation;

        EigenVector3<T> T_b2w;
        EigenQuaternion<T> Q_b2w;

        compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w,
                                           Q_l2w, T_b2w, Q_b2w);

        const size_t rid = create_rigid_body<T>(
            engine, T_b2w, Q_b2w, stone_handle, mid, stone_density);
        engine->set_rigid_body_velocity(rid, 0.0, -50.0, 0.0);
    }
}

template <typename T>
void make_ccd_cones_fast_towards_eachother(
    content::API* engine, const EigenVector3<T>& position,
    const EigenQuaternion<T>& orientation, const T& scene_size,
    MaterialInfo<T> const& mat_info,
    mesh_array::TetGenSettings tetset = mesh_array::tetgen_default_settings())
{
    ;
}

template <typename T>
void make_ccd_multiple_bounces_in_one_iteration(
    content::API* engine, const EigenVector3<T>& position,
    const EigenQuaternion<T>& orientation, const T& scene_size,
    MaterialInfo<T> const& mat_info,
    mesh_array::TetGenSettings tetset = mesh_array::tetgen_default_settings())
{
    T width = 5.0f;
    T height = 5.0f;
    T depth = 0.1f;
    T const H = 6.0 * scene_size / 50.0;

    T const stone_density = get_material_density_eigen<T>(mat_info, "Stone");
    size_t const mid = get_material_id_eigen<T>(mat_info, "Stone");

    GeometryHandleEigen<T> stone_handle
        = create_geometry_handle_box<T>(engine, width, height, depth);

    {
        T y = H;

        const EigenVector3<T> T_b2m = stone_handle.Tb2m();
        const EigenQuaternion<T> Q_b2m = stone_handle.Qb2m();

        const EigenVector3<T> T_m2l = EigenVector3<T>(0.0, 2.5, 0.0);
        const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

        const EigenVector3<T> T_l2w = position;
        const EigenQuaternion<T> Q_l2w = orientation;

        EigenVector3<T> T_b2w;
        EigenQuaternion<T> Q_b2w;

        compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w,
                                           Q_l2w, T_b2w, Q_b2w);

        const size_t rid = create_rigid_body<T>(
            engine, T_b2w, Q_b2w, stone_handle, mid, stone_density);
        engine->set_rigid_body_velocity(rid, 0.0, 0.0, 0.0);
    }
    {
        //We may have to increase speed when doing CCD, o.w. it may not be in same iteration we have multiple bounces!
        make_cannonball<T>(engine, 0.5f, EigenVector3<T>(0.0, H * 2.0, -20.0),
                           EigenQuaternion<T>::Identity(),
                           EigenVector3<T>(0.0, -50.0, 200.0), mat_info);
    }
}

template <typename T>
void make_ccd_small_space(content::API* engine, const EigenVector3<T>& position,
                          const EigenQuaternion<T>& orientation,
                          const T& scene_size, MaterialInfo<T> const& mat_info,
                          mesh_array::TetGenSettings tetset
                          = mesh_array::tetgen_default_settings())
{
    T width = 5.0f;
    T height = 0.1f;
    T depth = 5.0f;
    T H = 6.0 * scene_size / 50.0;

    T const stone_density = get_material_density_eigen<T>(mat_info, "Stone");
    size_t const mid = get_material_id_eigen<T>(mat_info, "Stone");

    GeometryHandleEigen<T> stone_handle
        = create_geometry_handle_box<T>(engine, width, height, depth);

    {
        T y = H;

        const EigenVector3<T> T_b2m = stone_handle.Tb2m();
        const EigenQuaternion<T> Q_b2m = stone_handle.Qb2m();

        const EigenVector3<T> T_m2l = EigenVector3<T>(0, y, 0);
        const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

        const EigenVector3<T> T_l2w = position;
        const EigenQuaternion<T> Q_l2w = orientation;

        EigenVector3<T> T_b2w;
        EigenQuaternion<T> Q_b2w;

        compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w,
                                           Q_l2w, T_b2w, Q_b2w);

        const size_t rid = create_rigid_body<T>(
            engine, T_b2w, Q_b2w, stone_handle, mid, stone_density);
        engine->set_rigid_body_velocity(rid, 0.0, 0.0, 0.0);
    }
    {
        make_cannonball<T>(engine, 0.5f,
                           EigenVector3<T>(0.0, H + 0.01f + 0.5f, 0.0),
                           EigenQuaternion<T>::Identity(),
                           EigenVector3<T>(0.0, 50.0, 200.0), mat_info);
    }

    {
        T y = H;

        const EigenVector3<T> T_b2m = stone_handle.Tb2m();
        const EigenQuaternion<T> Q_b2m = stone_handle.Qb2m();

        const EigenVector3<T> T_m2l
            = EigenVector3<T>(0, y + 0.1 + 1.0f + 0.0f, 0);
        const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

        const EigenVector3<T> T_l2w = position;
        const EigenQuaternion<T> Q_l2w = orientation;

        EigenVector3<T> T_b2w;
        EigenQuaternion<T> Q_b2w;

        compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w,
                                           Q_l2w, T_b2w, Q_b2w);

        const size_t rid = create_rigid_body<T>(
            engine, T_b2w, Q_b2w, stone_handle, mid, stone_density);
        engine->set_rigid_body_velocity(rid, 0.0, 0.0, 0.0);
    }
}

template <typename T>
void make_ccd_large_large_fast(content::API* engine,
                               const EigenVector3<T>& position,
                               const EigenQuaternion<T>& orientation,
                               const T& scene_size,
                               MaterialInfo<T> const& mat_info,
                               mesh_array::TetGenSettings tetset
                               = mesh_array::tetgen_default_settings())
{
    T width = 5.0f;
    T height = 0.1f;
    T depth = 5.0f;
    T const H = 6.0 * scene_size / 50.0;

    T const stone_density = get_material_density_eigen<T>(mat_info, "Stone");
    size_t const mid = get_material_id_eigen<T>(mat_info, "Stone");

    GeometryHandleEigen<T> stone_handle
        = create_geometry_handle_box<T>(engine, width, height, depth);

    {
        T y = H;

        const EigenVector3<T> T_b2m = stone_handle.Tb2m();
        const EigenQuaternion<T> Q_b2m = stone_handle.Qb2m();

        const EigenVector3<T> T_m2l = EigenVector3<T>(0, y, 0);
        const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

        const EigenVector3<T> T_l2w = position;
        const EigenQuaternion<T> Q_l2w = orientation;

        EigenVector3<T> T_b2w;
        EigenQuaternion<T> Q_b2w;

        compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w,
                                           Q_l2w, T_b2w, Q_b2w);

        const size_t rid = create_rigid_body<T>(
            engine, T_b2w, Q_b2w, stone_handle, mid, stone_density);
        engine->set_rigid_body_velocity(rid, 0.0, 0.0, 0.0);
    }

    {
        T y = H;

        const EigenVector3<T> T_b2m = stone_handle.Tb2m();
        const EigenQuaternion<T> Q_b2m = stone_handle.Qb2m();

        const EigenVector3<T> T_m2l = EigenVector3<T>(0, y * 2.0, 0);
        const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

        const EigenVector3<T> T_l2w = position;
        const EigenQuaternion<T> Q_l2w = orientation;

        EigenVector3<T> T_b2w;
        EigenQuaternion<T> Q_b2w;

        compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w,
                                           Q_l2w, T_b2w, Q_b2w);

        const size_t rid = create_rigid_body<T>(
            engine, T_b2w, Q_b2w, stone_handle, mid, stone_density);
        engine->set_rigid_body_velocity(rid, 0.0, -200.0, 0.0);
    }
}

} // namespace procedural

#endif // PROCEDURAL_MAKE_CCD_SDF_SCENES_H
