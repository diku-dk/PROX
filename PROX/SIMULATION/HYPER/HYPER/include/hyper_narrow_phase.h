#ifndef HYPER_NARROW_PHASE_H
#define HYPER_NARROW_PHASE_H

#include <hyper_engine.h>
#include <kdop_tandem_traversal.h>

namespace hyper
{

  template< typename MT>
  inline void narrow_phase(
                           Engine<MT> & engine
                           , std::vector< kdop::TestPair<typename MT::vector3_type, 8, typename MT::real_type> > & work_pool
                           )
  {
    typedef typename MT::real_type    T;
    typedef typename MT::vector3_type V;

    kdop::tandem_traversal<V, 8, T>(
                                    work_pool
                                    , kdop::sequential()
                                    );

  }

} //namespace narrow

// HYPER_NARROW_PHASE_H
#endif
