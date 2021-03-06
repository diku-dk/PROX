#ifndef NARROW_TETRAMESH_TETRAMESH_H
#define NARROW_TETRAMESH_TETRAMESH_H

#include "narrow_object.h"
#include "narrow_geometry.h"

#include <kdop_tandem_traversal.h>

namespace narrow
{

  namespace details
  {

    template< typename M>
    inline void dispatch_tetramesh_tetramesh( System<M> const & system, std::vector< TestPair<M> > & test_pairs )
    {
      assert( ! test_pairs.empty() || !"dispatch_tetramesh_tetramesh : test_pairs are empty" );
             
      typedef typename M::vector3_type                     V;
      typedef typename M::real_type                        T;
      typedef typename std::vector<TestPair<M> >::iterator pair_iterator;
      typedef typename kdop::TestPair<V, 8, T>             kdop_pair_type;

      pair_iterator current = test_pairs.begin();
      pair_iterator end     = test_pairs.end();

      std::vector< kdop_pair_type > kdop_test_pairs;

      for(;current!=end;++current)
      {
        Object<M> const & objA = current->obj_a();
        Object<M> const & objB = current->obj_b();

        Geometry<M> const & geoA = system.get_geometry( objA.get_geometry_idx() );
        Geometry<M> const & geoB = system.get_geometry( objB.get_geometry_idx() );

        kdop_pair_type const test_pair = kdop_pair_type(
                                                        objA.m_tree
                                                        , objB.m_tree
                                                        , geoA.m_tetramesh.m_mesh
                                                        , geoB.m_tetramesh.m_mesh
                                                        , objA.m_X
                                                        , objB.m_X
                                                        , objA.m_Y
                                                        , objB.m_Y
                                                        , objA.m_Z
                                                        , objB.m_Z
                                                        , geoA.m_tetramesh.m_surface_map
                                                        , geoB.m_tetramesh.m_surface_map
                                                        , current->callback()
                                                        );

        kdop_test_pairs.push_back( test_pair );
      }

      kdop::tandem_traversal<V, 8, T>(  kdop_test_pairs
                                        , kdop::sequential()
                                        );
              
    }

  } // namespace details
  
} //namespace narrow

// NARROW_TETRAMESH_TETRAMESH_H
#endif
