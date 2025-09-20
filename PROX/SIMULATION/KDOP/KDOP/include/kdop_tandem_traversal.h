#ifndef KDOP_TANDEM_TRAVERSAL_H
#define KDOP_TANDEM_TRAVERSAL_H

#include <grid.h>
#include <kdop_tags.h>
#include <kdop_test_pair.h>
#include <kdop_tree.h>
#include <kdop_select_contact_point_algorithm.h>
#include <types/geometry_dop.h>

#include <mesh_array.h>

#include <util_profiling.h>
#include <eigenhelperall.h>

namespace kdop
{
  namespace details
  {
  template <size_t K, typename T>
  inline void traversal(
      size_t const& node_idx_A, SubTree<T, K> const& branch_A,
      mesh_array::T4Mesh const& mesh_A,
      mesh_array::VertexAttribute<T, mesh_array::T4Mesh> const& X_A,
      mesh_array::VertexAttribute<T, mesh_array::T4Mesh> const& Y_A,
      mesh_array::VertexAttribute<T, mesh_array::T4Mesh> const& Z_A,
      mesh_array::TetrahedronAttribute<mesh_array::TetrahedronSurfaceInfo,
                                       mesh_array::T4Mesh> const& surface_map_A,
      size_t const& node_idx_B, SubTree<T, K> const& branch_B,
      mesh_array::T4Mesh const& mesh_B,
      mesh_array::VertexAttribute<T, mesh_array::T4Mesh> const& X_B,
      mesh_array::VertexAttribute<T, mesh_array::T4Mesh> const& Y_B,
      mesh_array::VertexAttribute<T, mesh_array::T4Mesh> const& Z_B,
      mesh_array::TetrahedronAttribute<mesh_array::TetrahedronSurfaceInfo,
                                       mesh_array::T4Mesh> const& surface_map_B,
      geometry::ContactsCallback<T>& callback)
  {
      using namespace mesh_array;

      Node<T,K> const & node_A = branch_A.m_nodes[node_idx_A];
      Node<T,K> const & node_B = branch_B.m_nodes[node_idx_B];

      if (!geometry::overlap_dop_dop(node_A.m_volume, node_B.m_volume)) return;

      bool const A_is_leaf = node_A.is_leaf();
      bool const B_is_leaf = node_B.is_leaf();

      if(A_is_leaf && B_is_leaf)
      {
        PAUSE_TIMER("tandem_traversal");
        RESUME_TIMER("exact_test");

        Tetrahedron const & tet_A = mesh_A.tetrahedron( node_A.m_start );
        Tetrahedron const & tet_B = mesh_B.tetrahedron( node_B.m_start );

        bool const & surface_Ai = surface_map_A( tet_A ).m_i;
        bool const & surface_Aj = surface_map_A( tet_A ).m_j;
        bool const & surface_Ak = surface_map_A( tet_A ).m_k;
        bool const & surface_Am = surface_map_A( tet_A ).m_m;

        if (
            !surface_Ai &&
            !surface_Aj &&
            !surface_Ak &&
            !surface_Am
            )
        {
          PAUSE_TIMER("exact_test");
          RESUME_TIMER("tandem_traversal");
          return;
        }

        bool const & surface_Bi = surface_map_B( tet_B ).m_i;
        bool const & surface_Bj = surface_map_B( tet_B ).m_j;
        bool const & surface_Bk = surface_map_B( tet_B ).m_k;
        bool const & surface_Bm = surface_map_B( tet_B ).m_m;

        if (
            !surface_Bi &&
            !surface_Bj &&
            !surface_Bk &&
            !surface_Bm
            )
        {
          PAUSE_TIMER("exact_test");
          RESUME_TIMER("tandem_traversal");
          return; // all faces of B are internal
        }

        const EigenVector3<T> a0
            = EigenVector3<T>(X_A(tet_A.i()), Y_A(tet_A.i()), Z_A(tet_A.i()));
        const EigenVector3<T> a1
            = EigenVector3<T>(X_A(tet_A.j()), Y_A(tet_A.j()), Z_A(tet_A.j()));
        const EigenVector3<T> a2
            = EigenVector3<T>(X_A(tet_A.k()), Y_A(tet_A.k()), Z_A(tet_A.k()));
        const EigenVector3<T> a3
            = EigenVector3<T>(X_A(tet_A.m()), Y_A(tet_A.m()), Z_A(tet_A.m()));

        const EigenVector3<T> b0
            = EigenVector3<T>(X_B(tet_B.i()), Y_B(tet_B.i()), Z_B(tet_B.i()));
        const EigenVector3<T> b1
            = EigenVector3<T>(X_B(tet_B.j()), Y_B(tet_B.j()), Z_B(tet_B.j()));
        const EigenVector3<T> b2
            = EigenVector3<T>(X_B(tet_B.k()), Y_B(tet_B.k()), Z_B(tet_B.k()));
        const EigenVector3<T> b3
            = EigenVector3<T>(X_B(tet_B.m()), Y_B(tet_B.m()), Z_B(tet_B.m()));

        std::vector<bool> surface_A( 4u, false );
        std::vector<bool> surface_B( 4u, false );

        geometry::TetrahedronEigen<T> const gtet_A
            = geometry::make_tetrahedron((a0), (a1), (a2), (a3));
        geometry::TetrahedronEigen<T> const gtet_B
            = geometry::make_tetrahedron((b0), (b1), (b2), (b3));

        surface_A[0] = surface_Ai;
        surface_A[1] = surface_Aj;
        surface_A[2] = surface_Ak;
        surface_A[3] = surface_Am;

        surface_B[0] = surface_Bi;
        surface_B[1] = surface_Bj;
        surface_B[2] = surface_Bk;
        surface_B[3] = surface_Bm;

        SelectContactPointAlgorithm::call_algorithm(gtet_A, gtet_B, callback, surface_A, surface_B );

        PAUSE_TIMER("exact_test");
        RESUME_TIMER("tandem_traversal");

      }
      else if(!A_is_leaf && !B_is_leaf)
      {
        for(size_t a = node_A.m_start; a <= node_A.m_end; ++a)
        {
          for(size_t b = node_B.m_start; b <= node_B.m_end; ++b)
          {
              traversal<K, T>(a, branch_A, mesh_A, X_A, Y_A, Z_A, surface_map_A,
                              b, branch_B, mesh_B, X_B, Y_B, Z_B, surface_map_B,
                              callback);
          }
        }
      }
      else if(!A_is_leaf && B_is_leaf)
      {
        for(size_t a = node_A.m_start; a <= node_A.m_end; ++a)
        {
            traversal<K, T>(a, branch_A, mesh_A, X_A, Y_A, Z_A, surface_map_A,
                            node_idx_B, branch_B, mesh_B, X_B, Y_B, Z_B,
                            surface_map_B, callback);
        }
      }
      else if( A_is_leaf && !B_is_leaf)
      {
        for(size_t b = node_B.m_start; b <= node_B.m_end; ++b)
        {
            traversal<K, T>(node_idx_A, branch_A, mesh_A, X_A, Y_A, Z_A,
                            surface_map_A, b, branch_B, mesh_B, X_B, Y_B, Z_B,
                            surface_map_B, callback);
        }
      }
  }

  }// namespace details

  template <size_t K, typename T>
  inline void tandem_traversal(TestPair<K, T>& work_item)
  {

    if(!geometry::overlap_dop_dop(work_item.m_tree_a->m_root, work_item.m_tree_b->m_root))
      return;

    size_t const C_A = work_item.m_tree_a->branches().size();
    size_t const C_B = work_item.m_tree_b->branches().size();

    for( size_t a = 0u; a < C_A; ++a)
    {
      SubTree<T,K> const & branch_A = work_item.m_tree_a->branches()[a];

      for( size_t b = 0u; b < C_B; ++b)
      {
        SubTree<T,K> const & branch_B = work_item.m_tree_b->branches()[b];

        details::traversal<K, T>(
            0, branch_A, *(work_item.m_mesh_a), *(work_item.m_x_a),
            *(work_item.m_y_a), *(work_item.m_z_a),
            *(work_item.m_surface_map_a), 0, branch_B, *(work_item.m_mesh_b),
            *(work_item.m_x_b), *(work_item.m_y_b), *(work_item.m_z_b),
            *(work_item.m_surface_map_b), *(work_item.m_callback));
      }
    }
  }

  template <size_t K, typename T>
  inline void tandem_traversal(std::vector<TestPair<K, T>>& work_pool,
                               sequential const& /*tag*/
  )
  {
    if( work_pool.empty() )
      return;

    typedef TestPair<K, T> work_item_type;
    typedef          std::vector< work_item_type >    work_pool_type;
    typedef typename work_pool_type::iterator         work_item_iterator;

    START_TIMER("tandem_traversal");
    START_TIMER("exact_test");
    PAUSE_TIMER("exact_test");

    work_item_iterator end     = work_pool.end();
    work_item_iterator current = work_pool.begin();

    for (; current != end; ++current) { tandem_traversal<K, T>(*current); }

    RESUME_TIMER("exact_test");
    STOP_TIMER("exact_test");
    STOP_TIMER("tandem_traversal");
  }

  }// namespace kdop

  template <typename T>
  Eigen::AlignedBox<T, 3>
  transformAABB(const Eigen::AlignedBox<T, 3>& aabb,
                const Eigen::Quaternion<T>& rotation,
                const Eigen::Matrix<T, 3, 1>& translation)
  {
      // Create transformation matrix from rotation and translation
      Eigen::Transform<T, 3, Eigen::Affine> transform
          = Eigen::Transform<T, 3, Eigen::Affine>::Identity();
      transform.rotate(rotation);
      transform.translate(translation);

      // Transform all 8 corners of the AABB and create a new AABB that contains them
      Eigen::AlignedBox<T, 3> transformed_aabb;

      for (int i = 0; i < 8; ++i)
      {
          // Get the corner using the enum values from Eigen::AlignedBox
          auto corner_type
              = static_cast<typename Eigen::AlignedBox<T, 3>::CornerType>(i);
          Eigen::Vector3<T> corner = aabb.corner(corner_type);
          transformed_aabb.extend(transform * corner);
      }

      return transformed_aabb;
  }

  namespace kdop
  {
  template <size_t K, typename T>
  inline bool
  overlap_dop_aabb(geometry::DOP<T, K> const& dop,
                   Eigen::AlignedBox<T, 3> const& aabb,
                   geometry::DirectionTable<T, K / 2> const& directions)
  {
      size_t const N = K / 2;

      // Get AABB center and half-extents
      Eigen::Vector3<T> center = (aabb.min() + aabb.max()) * 0.5;
      Eigen::Vector3<T> half_extents = (aabb.max() - aabb.min()) * 0.5;

      for (size_t k = 0u; k < N; ++k)
      {
          // Get the axis direction
          Eigen::Vector3<T> axis = directions(k);

          // Project AABB onto the axis
          T projection_center = center.dot(axis);
          T projection_radius = half_extents.x() * std::abs(axis.x())
                              + half_extents.y() * std::abs(axis.y())
                              + half_extents.z() * std::abs(axis.z());

          // Create interval for AABB projection
          geometry::Interval<T> aabb_interval;
          aabb_interval.lower() = (projection_center - projection_radius);

          aabb_interval.upper() = (projection_center + projection_radius);

          // Get DOP interval for this axis
          geometry::Interval<T> dop_interval = dop(k);

          // Check for separation along this axis
          if (!overlap_interval_interval(aabb_interval, dop_interval))
              return false;
      }

      return true;
  }

  template <size_t K, typename T>
  inline void traversal_sdf(
      size_t const& node_idx_A, SubTree<T, K> const& branch_A,
      mesh_array::T4Mesh const& mesh_A,
      mesh_array::VertexAttribute<T, mesh_array::T4Mesh> const& X_A,
      mesh_array::VertexAttribute<T, mesh_array::T4Mesh> const& Y_A,
      mesh_array::VertexAttribute<T, mesh_array::T4Mesh> const& Z_A,
      mesh_array::TetrahedronAttribute<mesh_array::TetrahedronSurfaceInfo,
                                       mesh_array::T4Mesh> const& surface_map_A,
      const grid::Grid<T, T>& sdf, const EigenVector3<T>& transformTranslation,
      const EigenQuaternion<T>& transformRotation,
      geometry::DirectionTable<T, K / 2> const& directions,
      geometry::ContactsCallback<T>& callback)
  {
      using namespace mesh_array;

      Node<T, K> const& node_A = branch_A.m_nodes[node_idx_A];

      //      if (!geometry::overlap_dop_dop(node_A.m_volume, node_B.m_volume)) return;

      bool const A_is_leaf = node_A.is_leaf();

      if (A_is_leaf)
      {
          PAUSE_TIMER("tandem_traversal");
          RESUME_TIMER("exact_test");

          Tetrahedron const& tet_A = mesh_A.tetrahedron(node_A.m_start);

          bool const& surface_Ai = surface_map_A(tet_A).m_i;
          bool const& surface_Aj = surface_map_A(tet_A).m_j;
          bool const& surface_Ak = surface_map_A(tet_A).m_k;
          bool const& surface_Am = surface_map_A(tet_A).m_m;

          if (!surface_Ai && !surface_Aj && !surface_Ak && !surface_Am)
          {
              PAUSE_TIMER("exact_test");
              RESUME_TIMER("tandem_traversal");
              return;
          }

          const EigenVector3<T> a0
              = EigenVector3<T>(X_A(tet_A.i()), Y_A(tet_A.i()), Z_A(tet_A.i()));
          const EigenVector3<T> a1
              = EigenVector3<T>(X_A(tet_A.j()), Y_A(tet_A.j()), Z_A(tet_A.j()));
          const EigenVector3<T> a2
              = EigenVector3<T>(X_A(tet_A.k()), Y_A(tet_A.k()), Z_A(tet_A.k()));
          const EigenVector3<T> a3
              = EigenVector3<T>(X_A(tet_A.m()), Y_A(tet_A.m()), Z_A(tet_A.m()));

          std::vector<bool> surface_A(4u, false);

          geometry::TetrahedronEigen<T> const gtet_A
              = geometry::make_tetrahedron((a0), (a1), (a2), (a3));

          surface_A[0] = surface_Ai;
          surface_A[1] = surface_Aj;
          surface_A[2] = surface_Ak;
          surface_A[3] = surface_Am;

          /*          SelectContactPointAlgorithm::call_algorithm(gtet_A, gtet_B, callback,
                                                      surface_A, surface_B);*/
          // Check each surface triangle against SDF

          geometry::Triangle<T> tri0 = geometry::get_opposite_face(0, gtet_A);
          geometry::Triangle<T> tri1 = geometry::get_opposite_face(1, gtet_A);
          geometry::Triangle<T> tri2 = geometry::get_opposite_face(2, gtet_A);
          geometry::Triangle<T> tri3 = geometry::get_opposite_face(3, gtet_A);

          if (surface_A[0] || true)
          { // Face opposite vertex i (vertices j,k,m)
              EigenVector3<T> contactPoint;
              EigenVector3<T> normal;
              T penetration;
              bool isPenetrating = grid::optimizeTriangleFWTransform<T, T>(
                  tri0.p(0), tri0.p(1), tri0.p(2), transformTranslation,
                  transformRotation, sdf, contactPoint, normal, penetration);
              if (isPenetrating) callback(contactPoint, normal, penetration);
          }
          if (surface_A[1] || true)
          { // Face opposite vertex j (vertices i,k,m)

              EigenVector3<T> contactPoint;
              EigenVector3<T> normal;
              T penetration;
              bool isPenetrating = grid::optimizeTriangleFWTransform<T, T>(
                  tri0.p(0), tri0.p(1), tri0.p(2), transformTranslation,
                  transformRotation, sdf, contactPoint, normal, penetration);
              if (isPenetrating) callback(contactPoint, normal, penetration);
          }
          if (surface_A[2] || true)
          { // Face opposite vertex k (vertices i,j,m)

              EigenVector3<T> contactPoint;
              EigenVector3<T> normal;
              T penetration;
              bool isPenetrating = grid::optimizeTriangleFWTransform<T, T>(
                  tri0.p(0), tri0.p(1), tri0.p(2), transformTranslation,
                  transformRotation, sdf, contactPoint, normal, penetration);
              if (isPenetrating) callback(contactPoint, normal, penetration);
          }
          if (surface_A[3] || true)
          { // Face opposite vertex m (vertices i,j,k)
              EigenVector3<T> contactPoint;
              EigenVector3<T> normal;
              T penetration;
              bool isPenetrating = grid::optimizeTriangleFWTransform<T, T>(
                  tri0.p(0), tri0.p(1), tri0.p(2), transformTranslation,
                  transformRotation, sdf, contactPoint, normal, penetration);
              if (isPenetrating) callback(contactPoint, normal, penetration);
          }

          PAUSE_TIMER("exact_test");
          RESUME_TIMER("tandem_traversal");
      }
      else if (!A_is_leaf)
      {
          for (size_t a = node_A.m_start; a <= node_A.m_end; ++a)
          {
              traversal_sdf<K, T>(a, branch_A, mesh_A, X_A, Y_A, Z_A,
                                  surface_map_A, sdf, transformTranslation,
                                  transformRotation, directions, callback);
          }
      }
  }

  // Updated traversal function for SDF vs TetraMesh
  template <size_t K, typename T>
  inline void traversal_sdf_old(
      size_t const& node_idx, SubTree<T, K> const& branch,
      mesh_array::T4Mesh const& mesh,
      mesh_array::VertexAttribute<T, mesh_array::T4Mesh> const& X,
      mesh_array::VertexAttribute<T, mesh_array::T4Mesh> const& Y,
      mesh_array::VertexAttribute<T, mesh_array::T4Mesh> const& Z,
      mesh_array::TetrahedronAttribute<mesh_array::TetrahedronSurfaceInfo,
                                       mesh_array::T4Mesh> const& surface_map,
      const grid::Grid<T, T>& sdf, const EigenVector3<T>& transformTranslation,
      const EigenQuaternion<T>& transformRotation,
      geometry::DirectionTable<T, K / 2> const& directions,
      geometry::ContactsCallback<T>& callback)
  {
      using namespace mesh_array;

      Node<T, K> const& node = branch.m_nodes[node_idx];

      // Get SDF bounding box (min and max corners)
/*      auto sdf_bbox = sdf.getBoundingBox();
      Eigen::AlignedBox<T, 3> sdf_aabb(
          Eigen::Vector3<T>(sdf_bbox.min_x, sdf_bbox.min_y, sdf_bbox.min_z),
          Eigen::Vector3<T>(sdf_bbox.max_x, sdf_bbox.max_y, sdf_bbox.max_z));*/
      Eigen::Matrix<T, 3, 1> sdfMin = sdf.min();
      Eigen::Matrix<T, 3, 1> sdfMax = sdf.max();
      Eigen::AlignedBox<T, 3> sdf_aabb(
          Eigen::Vector3<T>(sdfMin.x(), sdfMin.y(), sdfMin.z()),
          Eigen::Vector3<T>(sdfMax.x(), sdfMax.y(), sdfMax.z()));

      // Apply SDF transform to the AABB
      //auto transform = work_item.m_sdf->getTransform();
      auto translation = transformTranslation;
      auto rotation = transformRotation;
      //sdf_aabb = transform * sdf_aabb;
      sdf_aabb = transformAABB(sdf_aabb, rotation, translation);
      //      sdf_aabb = transform * sdf_aabb;

      // Check if node's DOP overlaps with SDF's AABB
      // if (!overlap_dop_aabb(node.m_volume, sdf_aabb, directions)) return;

      if (node.is_leaf())
      {
          PAUSE_TIMER("tandem_traversal");
          RESUME_TIMER("exact_test");

          Tetrahedron const& tet = mesh.tetrahedron(node.m_start);

          // Check if any face is a surface face
          bool const& surface_i = surface_map(tet).m_i;
          bool const& surface_j = surface_map(tet).m_j;
          bool const& surface_k = surface_map(tet).m_k;
          bool const& surface_m = surface_map(tet).m_m;

          if (!surface_i && !surface_j && !surface_k && !surface_m)
          {
              PAUSE_TIMER("exact_test");
              RESUME_TIMER("tandem_traversal");
              return; // all faces are internal
          }

          // Get vertex positions
          const EigenVector3<T> v0
              = EigenVector3<T>(X(tet.i()), Y(tet.i()), Z(tet.i()));
          const EigenVector3<T> v1
              = EigenVector3<T>(X(tet.j()), Y(tet.j()), Z(tet.j()));
          const EigenVector3<T> v2
              = EigenVector3<T>(X(tet.k()), Y(tet.k()), Z(tet.k()));
          const EigenVector3<T> v3
              = EigenVector3<T>(X(tet.m()), Y(tet.m()), Z(tet.m()));

          std::cerr << "V0: " << v0 << "\n";
          std::cerr << "V1: " << v1 << "\n";
          std::cerr << "V2: " << v2 << "\n";
          std::cerr << "V3: " << v3 << "\n";
          Eigen::Vector3<T> vRes = v0 + v1 + v2 + v3;

          // Check each surface triangle against SDF
          if (surface_i)
          { // Face opposite vertex i (vertices j,k,m)
              EigenVector3<T> contactPoint;
              EigenVector3<T> normal;
              T penetration;
              bool isPenetrating = grid::optimizeTriangleFWTransform<T, T>(
                  v1, v2, v3, translation, rotation, sdf, contactPoint, normal,
                  penetration);
              if (isPenetrating) callback(contactPoint, normal, penetration);
          }
          if (surface_j)
          { // Face opposite vertex j (vertices i,k,m)

              EigenVector3<T> contactPoint;
              EigenVector3<T> normal;
              T penetration;
              bool isPenetrating = grid::optimizeTriangleFWTransform<T, T>(
                  v0, v1, v2, translation, rotation, sdf, contactPoint, normal,
                  penetration);
              if (isPenetrating) callback(contactPoint, normal, penetration);
          }
          if (surface_k)
          { // Face opposite vertex k (vertices i,j,m)

              EigenVector3<T> contactPoint;
              EigenVector3<T> normal;
              T penetration;
              bool isPenetrating = grid::optimizeTriangleFWTransform<T, T>(
                  v0, v1, v3, translation, rotation, sdf, contactPoint, normal,
                  penetration);
              if (isPenetrating) callback(contactPoint, normal, penetration);
          }
          if (surface_m)
          { // Face opposite vertex m (vertices i,j,k)
              EigenVector3<T> contactPoint;
              EigenVector3<T> normal;
              T penetration;
              bool isPenetrating = grid::optimizeTriangleFWTransform<T, T>(
                  v0, v1, v2, translation, rotation, sdf, contactPoint, normal,
                  penetration);
              if (isPenetrating) callback(contactPoint, normal, penetration);
          }

          PAUSE_TIMER("exact_test");
          RESUME_TIMER("tandem_traversal");
      }
      else
      {
          // Recursively process child nodes
          for (size_t i = node.m_start; i <= node.m_end; ++i)
          {
              traversal_sdf<K, T>(i, branch, mesh, X, Y, Z, surface_map, sdf,
                                  transformTranslation, transformRotation,
                                  directions, callback);
          }
      }
  }

  // Top-level function for SDF vs TetraMesh collision
  template <size_t K, typename T>
  inline void tandem_traversal_sdf(kdop::TestPairSDFStruct<K, T>& work_item)
  {
      if (!work_item.m_tree_a || !work_item.m_grid_b) return;

      // Get SDF bounding box and transform it
      Eigen::Matrix<T, 3, 1> sdfMin = work_item.m_grid_b->min();
      Eigen::Matrix<T, 3, 1> sdfMax = work_item.m_grid_b->max();
      Eigen::AlignedBox<T, 3> sdf_aabb(
          Eigen::Vector3<T>(sdfMin.x(), sdfMin.y(), sdfMin.z()),
          Eigen::Vector3<T>(sdfMax.x(), sdfMax.y(), sdfMax.z()));

      //auto transform = work_item.m_sdf->getTransform();
      const Eigen::Matrix<T, 3, 1>* translation
          = work_item.m_transformTranslation_b;
      const Eigen::Quaternion<T>* rotation = work_item.m_transformRotation_b;
      //sdf_aabb = transform * sdf_aabb;
      sdf_aabb = transformAABB(sdf_aabb, *rotation, *translation);

      geometry::DirectionTable<T, K / 2> directions
          = geometry::DirectionTableHelper<T, K / 2>::make();
      // Check root-level overlap
      //if (!overlap_dop_aabb(work_item.m_tree_a->m_root, sdf_aabb, directions))
      //return;

      // Process all branches
      for (auto const& branch : work_item.m_tree_a->branches())
      {
          traversal_sdf<K, T>(0, branch, *(work_item.m_mesh_a),
                              *(work_item.m_x_a), *(work_item.m_y_a),
                              *(work_item.m_z_a), *(work_item.m_surface_map_a),
                              *(work_item.m_grid_b),
                              *(work_item.m_transformTranslation_b),
                              *(work_item.m_transformRotation_b), directions,
                              *(work_item.m_callback));
      }
  }

  template <size_t K, typename T>
  inline void
  tandem_traversal_sdf(std::vector<kdop::TestPairSDFStruct<K, T>>& work_pool,
                       sequential const& /*tag*/
  )
  {
      if (work_pool.empty()) return;

      START_TIMER("tandem_traversal");
      START_TIMER("exact_test");
      PAUSE_TIMER("exact_test");

      for (auto& item : work_pool) { tandem_traversal_sdf<K, T>(item); }

      RESUME_TIMER("exact_test");
      STOP_TIMER("exact_test");
      STOP_TIMER("tandem_traversal");
  }
  } // namespace kdop

  // KDOP_TANDEM_TRAVERSAL_H
#endif
