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

      using M = tiny::MathTypes<T>;
      typedef typename M::vector3_type V;
      typedef typename kdop::TestPair<8, T> kdop_pair_type;

      std::vector< kdop_pair_type > kdop_test_pairs;

      for (auto& current : test_pairs)
      {
          const auto& objA = current.obj_a();
          const auto& objB = current.obj_b();

          const auto& geoA = system.get_geometry(objA.get_geometry_idx());
          const auto& geoB = system.get_geometry(objB.get_geometry_idx());

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
        } else {

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

} //namespace narrow

// NARROW_TETRAMESH_TETRAMESH_H
#endif
