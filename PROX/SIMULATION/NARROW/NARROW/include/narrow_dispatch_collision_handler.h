#ifndef NARROW_DISPATCH_COLLISION_HANDLER_H
#define NARROW_DISPATCH_COLLISION_HANDLER_H

#include <narrow_test_pair.h>
#include <narrow_box_box.h>
#include <narrow_box_sphere.h>
#include <narrow_sphere_box.h>
#include <narrow_sphere_sphere.h>
#include <narrow_tetramesh_tetramesh.h>
#include <narrow_spheres_tetramesh.h>

#include <cassert>
#include <vector>

namespace narrow
{

  namespace details
  {
  template <typename T>
  inline void dispatch_primitives(System<T> const& system,
                                  std::vector<TestPair<T>>& test_pairs)
  {
      assert( ! test_pairs.empty() || !"dispatch_primitives : test_pairs are empty" );
      using M = tiny::MathTypes<T>;
      for (auto& elem : test_pairs)
      {
          const auto& geoA
              = system.get_geometry(elem.obj_a().get_geometry_idx());
          const auto& geoB
              = system.get_geometry(elem.obj_b().get_geometry_idx());
          auto const& boxesA = geoA.m_boxes;
          auto const& spheresA = geoA.m_spheres;
          //convex_container   const & hullsA     = geoA.m_hulls;

          auto const& boxesB = geoB.m_boxes;
          auto const& spheresB = geoB.m_spheres;
          //convex_container   const & hullsB     = geoB.m_hulls;

          detail::box_box<T>(boxesA, boxesB, elem.t_a(), elem.Q_a(), elem.t_b(),
                             elem.Q_b(), system.params().get_envelope(),
                             elem.callback());

          detail::box_sphere<T>(
              boxesA, spheresB, elem.t_a(), elem.Q_a(), elem.t_b(), elem.Q_b(),
              system.params().get_envelope(), elem.callback());

          detail::sphere_box<T>(
              spheresA, boxesB, elem.t_a(), elem.Q_a(), elem.t_b(), elem.Q_b(),
              system.params().get_envelope(), elem.callback());

          detail::sphere_sphere<T>(
              spheresA, spheresB, elem.t_a(), elem.Q_a(), elem.t_b(),
              elem.Q_b(), system.params().get_envelope(), elem.callback());
          /*
        detail::convex_convex<M>(
                           hullsA
                           , hullsB
                           , elem.t_a()
                           , elem.Q_a()
                           , elem.t_b()
                           , elem.Q_b()
                           , system.params().get_envelope()
                           , elem.callback()
                           );
 */
      }
  }

  template <typename T>
  inline void dispatch_mixed(System<T> const& system,
                             std::vector<TestPair<T>>& test_pairs)
  {
      assert( ! test_pairs.empty() || !"dispatch_mixed : test_pairs are empty" );
      using M = tiny::MathTypes<T>;
      for (auto& elem : test_pairs)
      {
          auto const& geoA
              = system.get_geometry(elem.obj_a().get_geometry_idx());
          auto const& geoB
              = system.get_geometry(elem.obj_b().get_geometry_idx());

          //box_container      const & boxesA     = geoA.m_boxes;
          auto const& spheresA = geoA.m_spheres;

          //box_container      const & boxesB     = geoB.m_boxes;
          auto const& spheresB = geoB.m_spheres;

          if (geoA.m_tetramesh.has_data())
          {
              detail::spheres_tetramesh<M>(spheresB, elem.t_b(), elem.Q_b(),
                                           elem.obj_a(), geoA, elem.callback(),
                                           true);
          }

        if ( geoB.m_tetramesh.has_data() )
        {
            detail::spheres_tetramesh<M>(spheresA, elem.t_a(), elem.Q_a(),
                                         elem.obj_b(), geoB, elem.callback(),
                                         false);
        }

      }
  }

  } //namespace details

  template <typename T>
  inline void
  dispatch_collision_handlers(System<T> const& system,
                              std::vector<TestPair<T>> const& test_pairs)
  {
      assert(!test_pairs.empty()
             || !"dispatch_collision_handlers : test_pairs are empty");

      std::vector<TestPair<T>> tetramesh_pairs;
      std::vector<TestPair<T>> primitive_pairs;
      std::vector<TestPair<T>> mixed_pairs;

      for (const auto& elem : test_pairs)
      {
          bool const A_is_tetramesh
              = system.get_geometry(elem.obj_a().get_geometry_idx())
                    .m_tetramesh.has_data();
          bool const B_is_tetramesh
              = system.get_geometry(elem.obj_b().get_geometry_idx())
                    .m_tetramesh.has_data();

          if (A_is_tetramesh && B_is_tetramesh)
          {
              tetramesh_pairs.push_back(elem);
          }
          else if (A_is_tetramesh && !B_is_tetramesh)
          {
              mixed_pairs.push_back(elem);
          }
          else if (!A_is_tetramesh && B_is_tetramesh)
          {
              mixed_pairs.push_back(elem);
          }
          else { primitive_pairs.push_back(elem); }
      }

      if (!primitive_pairs.empty())
          details::dispatch_primitives(system, primitive_pairs);
      if (!tetramesh_pairs.empty())
          details::dispatch_tetramesh_tetramesh(system, tetramesh_pairs);
      if (!mixed_pairs.empty()) details::dispatch_mixed(system, mixed_pairs);
  }

} //namespace narrow

// NARROW_DISPATCH_COLLISION_HANDLER_H
#endif
