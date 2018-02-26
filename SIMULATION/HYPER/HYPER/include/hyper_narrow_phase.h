#ifndef HYPER_NARROW_PHASE_H
#define HYPER_NARROW_PHASE_H

#include <hyper_engine.h>

#ifdef HAS_DIKUCL
#include <cl/gproximity/kdop_cl_gproximity_tandem_traversal.h>
#include <cl/kdop_cl_tandem_traversal.h>
#endif // HAS_DIKUCL

#include <kdop_tandem_traversal.h>

#include <util_profiling.h>

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

    START_TIMER("hyper_tandem_traversal");
#ifdef HAS_DIKUCL

    if(engine.params().use_open_cl())
    {
      if(engine.params().use_gproximity())
      {

        kdop::tandem_traversal<V, 8, T>(
                                        work_pool
                                        , kdop::dikucl::gproximity()
                                        , engine.params().open_cl_platform()
                                        , engine.params().open_cl_device());

      } else {

        kdop::tandem_traversal<V, 8, T>(
                                        work_pool
                                        , kdop::dikucl()
                                        , engine.params().open_cl_platform()
                                        , engine.params().open_cl_device());

      }
    } else {
#endif // HAS_DIKUCL

      kdop::tandem_traversal<V, 8, T>(
                                      work_pool
                                      , kdop::sequential()
                                      );

#ifdef HAS_DIKUCL
    }
#endif // HAS_DIKUCL
    STOP_TIMER("hyper_tandem_traversal");
  }

} //namespace narrow

// HYPER_NARROW_PHASE_H
#endif
