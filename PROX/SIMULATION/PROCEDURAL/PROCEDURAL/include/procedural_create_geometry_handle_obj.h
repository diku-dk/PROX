#ifndef PROCEDURAL_CREATE_GEOMETRY_HANDLE_OBJ_H
#define PROCEDURAL_CREATE_GEOMETRY_HANDLE_OBJ_H

#include "grid_helpers.h"
#include <content.h>

#include <procedural_factory_types.h>
#include <procedural_types.h>
#include <procedural_noise.h>

#include <eigenhelperall.h>
#include <mesh_array.h>
#include <mass.h>
#include <util_string_helper.h>

#include <string>

namespace procedural
{
template <typename T>
inline GeometryHandleEigen<T>
create_geometry_handle_obj(content::API* engine, std::string const& rel_file_name, const T& scale_x, const T& scale_y,
                           const T& scale_z, mesh_array::TetGenSettings tetset = mesh_array::tetgen_default_settings())
{
    using std::min;

    static size_t counter = 0u;

    std::string const geom_name = "nonuniform_obj_" + util::to_string(counter++);

    typename detail::MeshData<T> data;

    mesh_array::read_obj(rel_file_name, data.m_mesh, data.m_X, data.m_Y, data.m_Z);

    mesh_array::scale<T>(scale_x, scale_y, scale_z, data.m_mesh, data.m_X, data.m_Y, data.m_Z);

    mass::Properties<T> props_mf = mass::compute_mesh<T>(1, data.m_mesh.triangle_size(), &data);
    mass::Properties<T> props_bf = mass::translate_to_body_frame(props_mf);
    mass::Properties<T> props = mass::rotate_to_body_frame(props_mf);

      // Change geometry from model space to body space
    const EigenVector3<T> d = -EigenVector3<T>(props_mf.m_x, props_mf.m_y, props_mf.m_z);
    const EigenQuaternion<T> R = EigenQuaternion<T>(props.m_Qs, props.m_Qx, props.m_Qy, props.m_Qz);

    mesh_array::translate<T>(d, data.m_mesh, data.m_X, data.m_Y, data.m_Z);
    mesh_array::rotate<T>((R).conjugate(), data.m_mesh, data.m_X, data.m_Y, data.m_Z);

    mesh_array::T4Mesh volume;
    mesh_array::VertexAttribute<T, mesh_array::T4Mesh> volX;
    mesh_array::VertexAttribute<T, mesh_array::T4Mesh> volY;
    mesh_array::VertexAttribute<T, mesh_array::T4Mesh> volZ;

    mesh_array::tetgen(data.m_mesh, data.m_X, data.m_Y, data.m_Z, volume, volX, volY, volZ, tetset);

    size_t const gid = engine->create_collision_geometry(geom_name);
    size_t const N = volume.vertex_size();
    size_t const K = volume.tetrahedron_size();

    std::vector<size_t> verts;
    verts.resize(N);

    for (size_t i = 0; i < N; ++i) { verts[i] = i; }

    std::vector<T> coords;
    coords.resize(3u * N);

    for (size_t i = 0; i < N; ++i)
    {
        mesh_array::Vertex const v = volume.vertex(i);
        coords[3 * i] = volX(v);
        coords[3 * i + 1] = volY(v);
        coords[3 * i + 2] = volZ(v);
    }

    std::vector<size_t> tets;
    tets.resize(4u * K);

    for (size_t i = 0; i < K; ++i)
    {
        mesh_array::Tetrahedron const t = volume.tetrahedron(i);
        tets[4 * i] = t.i();
        tets[4 * i + 1] = t.j();
        tets[4 * i + 2] = t.k();
        tets[4 * i + 3] = t.m();
    }

    engine->set_tetramesh_shape(gid, N, K, &verts[0], &tets[0], &coords[0]);

    grid::Grid<T, T> grid;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    grid::build_VF_from_T3Mesh(data.m_mesh, data.m_X, data.m_Y, data.m_Z, V, F);
    //grid::build_VF_from_tris_dedup_exact(tris, V, F);
    grid = grid::projectGridToSDF<T, T>(
        V, F, Eigen::Matrix<size_t, 3, 1>(64, 64, 64));
    grid::extractIsosurfaceFromGrid<T, T>(grid);
    grid.m_temporaryGridStructure = grid::build_triangle_list_from_T3Mesh(
        data.m_mesh, data.m_X, data.m_Y, data.m_Z);
    //ProxData::geometry_type& geometry = engine->get_geometry_type(gid);
    //geometry.m_signedDistanceMap.setSignedDistanceGrid(grid);
    engine->set_geometry_type(gid, grid);

    return GeometryHandleEigen<T>(props.m_m, props.m_Ixx, props.m_Iyy, props.m_Izz,
                                  EigenVector3<T>(props_mf.m_x, props_mf.m_y, props_mf.m_z),
                                  EigenQuaternion<T>(props.m_Qs, props.m_Qx, props.m_Qy, props.m_Qz), gid);
}

  template <typename T>
  inline GeometryHandleEigen<T> create_geometry_handle_obj(content::API* engine, std::string const& rel_file_name,
                                                           const T& scale = 1, bool const& blind_copy = false,
                                                           mesh_array::TetGenSettings tetset
                                                           = mesh_array::tetgen_default_settings())
  {
      using std::min;

      static size_t counter = 0u;

      std::string const geom_name = "obj_" + util::to_string(counter++);

      detail::MeshData<T> data;

      mesh_array::read_obj(rel_file_name, data.m_mesh, data.m_X, data.m_Y, data.m_Z);

      if (!blind_copy)
      {
          T max_x = 0;
          T min_x = 0;

          T max_y = 0;
          T min_y = 0;

          T max_z = 0;
          T min_z = 0;

          mesh_array::compute_extents<T>(data.m_mesh, data.m_X, min_x, max_x);
          mesh_array::compute_extents<T>(data.m_mesh, data.m_Y, min_y, max_y);
          mesh_array::compute_extents<T>(data.m_mesh, data.m_Z, min_z, max_z);

          T const x_scale = 1 / (max_x - min_x);
          T const y_scale = 1 / (max_y - min_y);
          T const z_scale = 1 / (max_z - min_z);

          T const uniform_scale = min(min(x_scale, y_scale), z_scale) * scale;

          mesh_array::scale<T>(uniform_scale, uniform_scale, uniform_scale, data.m_mesh, data.m_X, data.m_Y, data.m_Z);
      }

      mass::Properties<T> props_mf = mass::compute_mesh<T>(1, data.m_mesh.triangle_size(), &data);
      mass::Properties<T> props_bf = mass::translate_to_body_frame(props_mf);
      mass::Properties<T> props = mass::rotate_to_body_frame(props_mf);

      // Change geometry from model space to body space
      if (!blind_copy)
      {

          const EigenVector3<T> d = -EigenVector3<T>(props_mf.m_x, props_mf.m_y, props_mf.m_z);
          const EigenQuaternion<T> R = EigenQuaternion<T>(props.m_Qs, props.m_Qx, props.m_Qy, props.m_Qz);

          mesh_array::translate<T>(d, data.m_mesh, data.m_X, data.m_Y, data.m_Z);
          mesh_array::rotate<T>((R).conjugate(), data.m_mesh, data.m_X, data.m_Y, data.m_Z);
      }

      if (blind_copy)
      {
          tetset.m_maximum_volume = 0;
          tetset.m_quality_ratio = 0;
          tetset.m_suppress_splitting = true;
      }

      mesh_array::T4Mesh volume;
      mesh_array::VertexAttribute<T, mesh_array::T4Mesh> volX;
      mesh_array::VertexAttribute<T, mesh_array::T4Mesh> volY;
      mesh_array::VertexAttribute<T, mesh_array::T4Mesh> volZ;

      mesh_array::tetgen(data.m_mesh, data.m_X, data.m_Y, data.m_Z, volume, volX, volY, volZ, tetset);

      size_t const gid = engine->create_collision_geometry(geom_name);
      size_t const N = volume.vertex_size();
      size_t const K = volume.tetrahedron_size();

      std::vector<size_t> verts;
      verts.resize(N);

      for (size_t i = 0; i < N; ++i) { verts[i] = i; }

      std::vector<T> coords;
      coords.resize(3u * N);

      for (size_t i = 0; i < N; ++i)
      {
          mesh_array::Vertex const v = volume.vertex(i);
          coords[3 * i] = volX(v);
          coords[3 * i + 1] = volY(v);
          coords[3 * i + 2] = volZ(v);
      }

      std::vector<size_t> tets;
      tets.resize(4u * K);

      for (size_t i = 0; i < K; ++i)
      {
          mesh_array::Tetrahedron const t = volume.tetrahedron(i);
          tets[4 * i] = t.i();
          tets[4 * i + 1] = t.j();
          tets[4 * i + 2] = t.k();
          tets[4 * i + 3] = t.m();
      }

      engine->set_tetramesh_shape(gid, N, K, &verts[0], &tets[0], &coords[0]);

      grid::Grid<T, T> grid;
      Eigen::MatrixXd V;
      Eigen::MatrixXi F;

      EigenVector3<T> centerTranslation
          = EigenVector3<T>(props_mf.m_x, props_mf.m_y, props_mf.m_z);
      grid::build_VF_from_T3Mesh(data.m_mesh, data.m_X, data.m_Y, data.m_Z, V,
                                 F);
      //grid::build_VF_from_tris_dedup_exact(tris, V, F);
      grid = grid::projectGridToSDF<T, T>(
          V, F, Eigen::Matrix<size_t, 3, 1>(64, 64, 64));
      grid::extractIsosurfaceFromGrid<T, T>(grid);
      grid.m_temporaryGridStructure = grid::build_triangle_list_from_T3Mesh(
          data.m_mesh, data.m_X, data.m_Y, data.m_Z);
      //ProxData::geometry_type& geometry = engine->get_geometry_type(gid);
      //geometry.m_signedDistanceMap.setSignedDistanceGrid(grid);
      engine->set_geometry_type(gid, grid);

      return GeometryHandleEigen<T>(
          props.m_m, props.m_Ixx, props.m_Iyy, props.m_Izz, centerTranslation,
          EigenQuaternion<T>(props.m_Qs, props.m_Qx, props.m_Qy, props.m_Qz),
          gid);
  }
  } // namespace procedural

// PROCEDURAL_CREATE_GEOMETRY_HANDLE_OBJ_H
#endif
