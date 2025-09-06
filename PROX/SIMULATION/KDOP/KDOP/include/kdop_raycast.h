#ifndef KDOP_RAYCAST_H
#define KDOP_RAYCAST_H

#include <kdop_tree.h>

#include <raycast/geometry_raycast_tetrahedron.h>
#include <raycast/geometry_raycast_dop.h>

#include <mesh_array.h>
#include <mesh_array_t4mesh.h>
#include <mesh_array_vertex_attribute.h>
#include <mesh_array_compute_surface_map.h>

namespace kdop
{
  namespace details
  {

  template <typename T, size_t K>
  inline void raycast(
      geometry::RayEigen<T> const& ray, size_t const& node_idx,
      SubTree<T, K> const& branch, mesh_array::T4Mesh const& mesh,
      mesh_array::VertexAttribute<T, mesh_array::T4Mesh> const& X,
      mesh_array::VertexAttribute<T, mesh_array::T4Mesh> const& Y,
      mesh_array::VertexAttribute<T, mesh_array::T4Mesh> const& Z,
      mesh_array::TetrahedronAttribute<mesh_array::TetrahedronSurfaceInfo,
                                       mesh_array::T4Mesh> const& surface_map,
      EigenVector3<T>& hit_point, T& length)
  {
      T const dop_threshold = 0.01;

      Node<T,K> const & node = branch.m_nodes[node_idx];

      EigenVector3<T> p = EigenVector3<T>(0,0,0);
      T t = std::numeric_limits<T>::max();

      if (!geometry::compute_raycast_dop<T, K>(ray, node.m_volume, p, t,
                                               dop_threshold))
          return;

      if(t >= length)
        return;

      bool const is_leaf = node.is_leaf();

      if(is_leaf)
      {
        mesh_array::Tetrahedron const mT = mesh.tetrahedron( node.m_start );

        EigenVector3<T> p0(X(mT.i()), Y(mT.i()), Z(mT.i()));
        EigenVector3<T> p1(X(mT.j()), Y(mT.j()), Z(mT.j()));
        EigenVector3<T> p2(X(mT.k()), Y(mT.k()), Z(mT.k()));
        EigenVector3<T> p3(X(mT.m()), Y(mT.m()), Z(mT.m()));

        std::vector<bool> surf(4u,false);

        surf[0] = surface_map(mT).m_i;
        surf[1] = surface_map(mT).m_j;
        surf[2] = surface_map(mT).m_k;
        surf[3] = surface_map(mT).m_m;

        EigenVector3<T> p = EigenVector3<T>(0,0,0);
        T s = std::numeric_limits<T>::max();

        geometry::TetrahedronEigen<T> const gT
            = geometry::make_tetrahedron((p0), (p1), (p2), (p3));

        bool const did_hit = geometry::compute_raycast_tetrahedron<T>(geometry::convertRayToEigen( ray), gT, p, s, surf);

        if(did_hit)
        {
            hit_point    = s < length ? fromEigen(p) : hit_point;
          length       = s < length ? s : length;
        }

      }
      else
      {

        for(size_t idx = node.m_start; idx <= node.m_end; ++idx)
        {
            raycast<T, K>(ray, idx, branch, mesh, X, Y, Z, surface_map,
                          hit_point, length);
        }

      }
  }

  } // end of namespace details

  template <typename T, size_t K>
  inline bool raycast(
      geometry::RayEigen<T> const& ray, Tree<T, K> const& tree,
      mesh_array::T4Mesh const& mesh,
      mesh_array::VertexAttribute<T, mesh_array::T4Mesh> const& X,
      mesh_array::VertexAttribute<T, mesh_array::T4Mesh> const& Y,
      mesh_array::VertexAttribute<T, mesh_array::T4Mesh> const& Z,
      mesh_array::TetrahedronAttribute<mesh_array::TetrahedronSurfaceInfo,
                                       mesh_array::T4Mesh> const& surface_map,
      EigenVector3<T>& hit_point, T& length)
  {
      T const dop_threshold = (0.01f);

      hit_point = {0, 0, 0};
      length = std::numeric_limits<T>::max();

      EigenVector3<T> p = EigenVector3<T>(0, 0, 0);
      T t = std::numeric_limits<T>::max();

      if (!geometry::compute_raycast_dop(ray, tree.m_root, p, t, dop_threshold))
          return false;

      size_t const C = tree.branches().size();

      for (size_t c = 0u; c < C; ++c)
      {
          SubTree<T, K> const& branch = tree.branches()[c];

          details::raycast<T, K>(ray, 0, branch, mesh, X, Y, Z, surface_map,
                                 hit_point, length);
      }

      return length < std::numeric_limits<T>::max();
  }

}// namespace kdop

// KDOP_RAYCAST_H
#endif
