#ifndef NARROW_TETRAMESH_TETRAMESH_H
#define NARROW_TETRAMESH_TETRAMESH_H

#include "narrow_object.h"
#include "narrow_geometry.h"
#include "narrow_system.h"
#include "narrow_test_pair.h"

#ifdef HAS_DIKUCL
#include <cl/gproximity/kdop_cl_gproximity_tandem_traversal.h>
#include <cl/kdop_cl_tandem_traversal.h>
#endif // HAS_DIKUCL

#include <kdop_tandem_traversal.h>

namespace narrow
{

  namespace details
  {

  template <typename T>
  inline void dispatch_tetramesh_tetramesh(System<T> const& system,
                                           std::vector<TestPair<T>>& test_pairs)
  {
      assert( ! test_pairs.empty() || !"dispatch_tetramesh_tetramesh : test_pairs are empty" );
      typedef typename kdop::TestPair<8, T> kdop_pair_type;
//      typedef typename kdop::TestPairSDF<8, T> kdop_sdf_pair_type;

      std::vector<kdop_pair_type> kdop_test_pairs;

      for (auto& current : test_pairs)
      {
          const auto& objA = current.obj_a();
          const auto& objB = current.obj_b();

          const auto& geoA = system.get_geometry(objA.get_geometry_idx());
          const auto& geoB = system.get_geometry(objB.get_geometry_idx());

          //Callback is a reference to our pairs, so we go from pairs to be tested to
          // making a kdop-pair-type
          kdop_pair_type const test_pair = kdop_pair_type(
              objA.m_tree, objB.m_tree, geoA.m_tetramesh.m_mesh,
              geoB.m_tetramesh.m_mesh, objA.m_X, objB.m_X, objA.m_Y, objB.m_Y,
              objA.m_Z, objB.m_Z, geoA.m_tetramesh.m_surface_map,
              geoB.m_tetramesh.m_surface_map, current.callback());

          kdop_test_pairs.push_back(test_pair);
      }
#ifdef HAS_DIKUCL

      if(system.params().use_open_cl())
      {
        if(system.params().use_gproximity())
        {

            kdop::tandem_traversal<8, T>(kdop_test_pairs,
                                         kdop::dikucl::gproximity(),
                                         system.params().open_cl_platform(),
                                         system.params().open_cl_device());
        }
        else
        {

            kdop::tandem_traversal<8, T>(kdop_test_pairs, kdop::dikucl(),
                                         system.params().open_cl_platform(),
                                         system.params().open_cl_device());
        }
      } else {
#endif // HAS_DIKUCL

        // use regular tandem traversal if DIKUCL is not available or should not be used
        kdop::tandem_traversal<8, T>(kdop_test_pairs, kdop::sequential());

#ifdef HAS_DIKUCL
      }
#endif // HAS_DIKUCL
  }

  } // namespace details

  namespace details
  {

  template <typename T>
  inline void dispatch_tetramesh_sdf(System<T> const& system,
                                     std::vector<TestPair<T>>& test_pairs)
  {
      assert(!test_pairs.empty()
             || !"dispatch_tetramesh_tetramesh : test_pairs are empty");
      //typedef typename kdop::TestPairSDF<8, T> kdop_sdf_pair_type;

      std::vector<kdop::TestPairSDFStruct<8, T>> kdop_test_sdf_pairs;
      kdop::TestPairSDFStruct<8, T> sdf_pair_type;

      for (auto& current : test_pairs)
      {
          const auto& objA = current.obj_a();
          const auto& objB = current.obj_b();

          const auto& geoA = system.get_geometry(objA.get_geometry_idx());
          const auto& geoB = system.get_geometry(objB.get_geometry_idx());

          //Callback is a reference to our pairs, so we go from pairs to be tested to
          // making a kdop-pair-type
          /*          const kdop_sdf_pair_type test_pair = kdop_sdf_pair_type(
              objA.m_tree, geoA.m_tetramesh.m_mesh, objA.m_X, objA.m_Y,
              objA.m_Z, geoA.m_tetramesh.m_surface_map,
              geoB.m_signedDistanceMap, current.callback());*/
          /*2
          sdf_pair_type.m_mesh_a = &geoA.m_tetramesh.m_mesh;
          sdf_pair_type.m_tree_a = &objA.m_tree;
          sdf_pair_type.m_x_a = &objA.m_X;
          sdf_pair_type.m_y_a = &objA.m_Y;
          sdf_pair_type.m_z_a = &objA.m_Z;
          sdf_pair_type.m_grid_b
              = &geoB.m_signedDistanceMap.getSignedDistanceGrid();
          sdf_pair_type.m_surface_map_a = &geoA.m_tetramesh.m_surface_map;
          sdf_pair_type.m_triangles_a
              = &geoA.m_signedDistanceMap.getSignedDistanceGrid()
                     .m_temporaryGridStructure;
          sdf_pair_type.m_transformTranslation_a = &current.t_a();
          sdf_pair_type.m_transformRotation_a = &current.Q_a();
          sdf_pair_type.m_transformTranslation_b = &current.t_b();
          sdf_pair_type.m_transformRotation_b = &current.Q_b();
          sdf_pair_type.m_callback = &current.callback();*/

          //kdop_test_sdf_pairs.push_back(sdf_pair_type);
          sdf_pair_type.m_mesh_a = &geoB.m_tetramesh.m_mesh;
          sdf_pair_type.m_tree_a = &objB.m_tree;
          sdf_pair_type.m_x_a = &objB.m_X;
          sdf_pair_type.m_y_a = &objB.m_Y;
          sdf_pair_type.m_z_a = &objB.m_Z;
          sdf_pair_type.m_grid_b
              = &geoA.m_signedDistanceMap.getSignedDistanceGrid();
          sdf_pair_type.m_surface_map_a = &geoB.m_tetramesh.m_surface_map;
          sdf_pair_type.m_triangles_a
              = &geoB.m_signedDistanceMap.getSignedDistanceGrid()
                     .m_temporaryGridStructure;
          sdf_pair_type.m_transformTranslation_a = &current.t_b();
          sdf_pair_type.m_transformRotation_a = &current.Q_b();
          sdf_pair_type.m_transformTranslation_b = &current.t_a();
          sdf_pair_type.m_transformRotation_b = &current.Q_a();
          sdf_pair_type.m_callback = &current.callback();
          kdop_test_sdf_pairs.push_back(sdf_pair_type);
      }
#ifdef HAS_DIKUCL

      if (system.params().use_open_cl())
      {
          if (system.params().use_gproximity())
          {

              kdop::tandem_traversal<8, T>(kdop_test_pairs,
                                           kdop::dikucl::gproximity(),
                                           system.params().open_cl_platform(),
                                           system.params().open_cl_device());
          }
          else
          {

              kdop::tandem_traversal<8, T>(kdop_test_pairs, kdop::dikucl(),
                                           system.params().open_cl_platform(),
                                           system.params().open_cl_device());
          }
      }
      else
      {
#endif // HAS_DIKUCL

          // use regular tandem traversal if DIKUCL is not available or should not be used
          kdop::tandem_traversal_sdf<8, T>(kdop_test_sdf_pairs,
                                           kdop::sequential());

#ifdef HAS_DIKUCL
      }
#endif // HAS_DIKUCL
  }

  template <typename T>
  inline T dispatch_tetramesh_sdf_CCD(
      System<T> const& system, std::vector<TestPairCCD<T>>& test_pairs,
      T startTime, T endTime,
      std::vector<kdop::BodyVelocities<T>>& bodyContacts, bool& onlyZEROTOI)
  {
      assert(!test_pairs.empty()
             || !"dispatch_tetramesh_tetramesh : test_pairs are empty");
      //typedef typename kdop::TestPairSDF<8, T> kdop_sdf_pair_type;

      std::vector<kdop::TestPairSDFStruct<8, T>> kdop_test_sdf_pairs;
      kdop::TestPairSDFStruct<8, T> sdf_pair_type;

      for (auto& current : test_pairs)
      {
          const auto& objA = current.obj_a();
          const auto& objB = current.obj_b();

          const auto& geoA = system.get_geometry(objA.get_geometry_idx());
          const auto& geoB = system.get_geometry(objB.get_geometry_idx());

          //Callback is a reference to our pairs, so we go from pairs to be tested to
          // making a kdop-pair-type
          /*          const kdop_sdf_pair_type test_pair = kdop_sdf_pair_type(
              objA.m_tree, geoA.m_tetramesh.m_mesh, objA.m_X, objA.m_Y,
              objA.m_Z, geoA.m_tetramesh.m_surface_map,
              geoB.m_signedDistanceMap, current.callback());*/
          /*2
          sdf_pair_type.m_mesh_a = &geoA.m_tetramesh.m_mesh;
          sdf_pair_type.m_tree_a = &objA.m_tree;
          sdf_pair_type.m_x_a = &objA.m_X;
          sdf_pair_type.m_y_a = &objA.m_Y;
          sdf_pair_type.m_z_a = &objA.m_Z;
          sdf_pair_type.m_grid_b
              = &geoB.m_signedDistanceMap.getSignedDistanceGrid();
          sdf_pair_type.m_surface_map_a = &geoA.m_tetramesh.m_surface_map;
          sdf_pair_type.m_triangles_a
              = &geoA.m_signedDistanceMap.getSignedDistanceGrid()
                     .m_temporaryGridStructure;
          sdf_pair_type.m_transformTranslation_a = &current.t_a();
          sdf_pair_type.m_transformRotation_a = &current.Q_a();
          sdf_pair_type.m_transformTranslation_b = &current.t_b();
          sdf_pair_type.m_transformRotation_b = &current.Q_b();
          sdf_pair_type.m_callback = &current.callback();*/

          //kdop_test_sdf_pairs.push_back(sdf_pair_type);
          //Below is swapped -- WHy did I do that?!
          /*          sdf_pair_type.m_mesh_a = &geoB.m_tetramesh.m_mesh;
          sdf_pair_type.m_tree_a = &objB.m_tree;
          sdf_pair_type.m_x_a = &objB.m_X;
          sdf_pair_type.m_y_a = &objB.m_Y;
          sdf_pair_type.m_z_a = &objB.m_Z;
          sdf_pair_type.m_grid_b
              = &geoA.m_signedDistanceMap.getSignedDistanceGrid();
          sdf_pair_type.m_surface_map_a = &geoB.m_tetramesh.m_surface_map;
          sdf_pair_type.m_triangles_a
              = &geoB.m_signedDistanceMap.getSignedDistanceGrid()
                     .m_temporaryGridStructure;
          sdf_pair_type.m_transformTranslation_a = &current.t_b();
          sdf_pair_type.m_transformRotation_a = &current.Q_b();
          sdf_pair_type.m_transformTranslation_b = &current.t_a();
          sdf_pair_type.m_transformRotation_b = &current.Q_a();
          sdf_pair_type.m_callback = nullptr;
          kdop_test_sdf_pairs.push_back(sdf_pair_type);*/
          sdf_pair_type.m_mesh_a = &geoA.m_tetramesh.m_mesh;
          sdf_pair_type.m_tree_a = &objA.m_tree;
          sdf_pair_type.m_x_a = &objA.m_X;
          sdf_pair_type.m_y_a = &objA.m_Y;
          sdf_pair_type.m_z_a = &objA.m_Z;
          sdf_pair_type.m_grid_b
              = &geoB.m_signedDistanceMap.getSignedDistanceGrid();
          sdf_pair_type.m_surface_map_a = &geoA.m_tetramesh.m_surface_map;
          sdf_pair_type.m_triangles_a
              = &geoA.m_signedDistanceMap.getSignedDistanceGrid()
                     .m_temporaryGridStructure;
          sdf_pair_type.m_transformTranslation_a = &current.t_a();
          sdf_pair_type.m_transformRotation_a = &current.Q_a();
          sdf_pair_type.m_transformTranslation_b = &current.t_b();
          sdf_pair_type.m_transformRotation_b = &current.Q_b();
          sdf_pair_type.m_callback = nullptr;
          kdop_test_sdf_pairs.push_back(sdf_pair_type);
      }

      std::vector<T> TOIs;
      kdop::tandem_traversal_sdf_CCD<8, T>(kdop_test_sdf_pairs,
                                           kdop::sequential(), TOIs, startTime,
                                           endTime, bodyContacts);
      onlyZEROTOI = false;
      if (TOIs.size() > 0)
      {
          T earliestTOI = endTime;
          for (size_t i = 0; i < TOIs.size(); ++i)
          {
              T currTOI = earliestTOI;
              if (std::abs<T>(TOIs[i]) > 0.000001) { currTOI = TOIs[i]; }
              else { onlyZEROTOI = true; }

              earliestTOI = std::min<T>(earliestTOI, currTOI);
          }
          return earliestTOI;
      }
      else { return std::numeric_limits<T>::max(); }
  }

  template <typename T>
  inline T dispatch_tetramesh_sdf_CCD_WARM_START(
      System<T> const& system, std::vector<TestPairCCD<T>>& test_pairs,
      T startTime, T endTime,
      std::vector<kdop::BodyVelocities<T>>& bodyContacts, bool& onlyZEROTOI,
      std::vector<std::vector<T>>& warmBodiesStart)
  {
      assert(!test_pairs.empty()
             || !"dispatch_tetramesh_tetramesh : test_pairs are empty");
      //typedef typename kdop::TestPairSDF<8, T> kdop_sdf_pair_type;

      std::vector<kdop::TestPairSDFStruct<8, T>> kdop_test_sdf_pairs;
      kdop::TestPairSDFStruct<8, T> sdf_pair_type;

      for (auto& current : test_pairs)
      {
          const auto& objA = current.obj_a();
          const auto& objB = current.obj_b();

          const auto& geoA = system.get_geometry(objA.get_geometry_idx());
          const auto& geoB = system.get_geometry(objB.get_geometry_idx());

          sdf_pair_type.m_mesh_a = &geoA.m_tetramesh.m_mesh;
          sdf_pair_type.m_tree_a = &objA.m_tree;
          sdf_pair_type.m_x_a = &objA.m_X;
          sdf_pair_type.m_y_a = &objA.m_Y;
          sdf_pair_type.m_z_a = &objA.m_Z;
          sdf_pair_type.m_grid_b
              = &geoB.m_signedDistanceMap.getSignedDistanceGrid();
          sdf_pair_type.m_surface_map_a = &geoA.m_tetramesh.m_surface_map;
          sdf_pair_type.m_triangles_a
              = &geoA.m_signedDistanceMap.getSignedDistanceGrid()
                     .m_temporaryGridStructure;
          sdf_pair_type.m_transformTranslation_a = &current.t_a();
          sdf_pair_type.m_transformRotation_a = &current.Q_a();
          sdf_pair_type.m_transformTranslation_b = &current.t_b();
          sdf_pair_type.m_transformRotation_b = &current.Q_b();
          sdf_pair_type.m_callback = nullptr;
          kdop_test_sdf_pairs.push_back(sdf_pair_type);
      }

      std::vector<T> TOIs;
      kdop::tandem_traversal_sdf_CCD_WARM_START<8, T>(
          kdop_test_sdf_pairs, kdop::sequential(), TOIs, startTime, endTime,
          bodyContacts, warmBodiesStart);
      onlyZEROTOI = false;
      if (TOIs.size() > 0)
      {
          T earliestTOI = endTime;
          for (size_t i = 0; i < TOIs.size(); ++i)
          {
              T currTOI = earliestTOI;
              if (std::abs<T>(TOIs[i]) > 0.000001) { currTOI = TOIs[i]; }
              else { onlyZEROTOI = true; }

              earliestTOI = std::min<T>(earliestTOI, currTOI);
          }
          return earliestTOI;
      }
      else { return std::numeric_limits<T>::max(); }
  }

  } // namespace details

} //namespace narrow

// NARROW_TETRAMESH_TETRAMESH_H
#endif
