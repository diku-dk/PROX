#ifndef HYPER_UPDATE_KDOP_BVH_H
#define HYPER_UPDATE_KDOP_BVH_H

#include <hyper_engine.h>

#include <kdop_refit_tree.h>

namespace hyper
{

  template<typename MT>
  inline void update_kdop_bvh( Engine<MT> & engine )
  {
    typedef typename MT::real_type    T;
    typedef typename MT::vector3_type V;

    typename Engine<MT>::body_iterator current = engine.body_begin();
    typename Engine<MT>::body_iterator end     = engine.body_end();

    for (; current != end; ++current)
    {
      Body<MT> & body    = *(current);

      if(body.empty())
        continue;

      kdop::refit_tree<V,8,T>(body.m_tree, body.m_mesh, body.m_X, body.m_Y, body.m_Z );
    }
    
  }
  
} // namespace hyper

// HYPER_UPDATE_KDOP_BVH_H
#endif
