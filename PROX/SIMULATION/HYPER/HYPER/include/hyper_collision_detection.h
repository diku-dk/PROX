#ifndef HYPER_COLLISION_DETECTION_H
#define HYPER_COLLISION_DETECTION_H

#include <hyper_engine.h>
#include <hyper_contact_callback.h>

#include <hyper_update_kdop_bvh.h>
#include <hyper_update_bounding_volumes.h>

#include <hyper_broad_phase.h>
#include <hyper_narrow_phase.h>

#include <util_profiling.h>

#include <cassert>
#include <vector>

namespace hyper
{

  /**
   * Inovoke collision detection system.
   * This function basically wraps the entire collision detection system into one functional ``unit''.
   *
   * @tparam contact_container   The type of contact container
   * @tparam M                   The math types used
   *
   * @param data                 Collision detection system data that are reused from invokation to invokation.
   * @param contacts             Upon return this container will hold all current contact points in the configuration.
   * @param tag                  Tag dispatcher for transfering math types to collision detection system.
   */
  template< typename MT>
  inline void collision_detection( Engine<MT> & engine )
  {
    typedef typename MT::real_type                         T;
    typedef typename MT::value_traits                      VT;
    typedef typename MT::vector3_type                      V;

    typedef          Engine<MT>                            engine_type;
    typedef typename engine_type::body_type                body_type;
    typedef typename engine_type::body_iterator            body_iterator;
//    typedef typename engine_type::const_body_iterator      const_body_iterator;
//    typedef typename engine_type::contact_type             contact_type;
//    typedef typename engine_type::contact_iterator         contact_iterator;
//    typedef typename engine_type::const_contact_iterator   const_contact_iterator;
//    typedef typename engine_type::params_type              params_type;

    typedef          detail::ContactCallbackFunctor<MT>     callback_type;

    START_TIMER("collision_detection");

    //--- First update collision detection data structures
    START_TIMER("collision_detection_updating_kdop");

    update_kdop_bvh( engine );

    STOP_TIMER("collision_detection_updating_kdop");
    START_TIMER("updating_bounding_sphere");

    for (body_iterator body = engine.body_begin();body != engine.body_end(); ++body)
    {
      update_bounding_volumes( *body );
    }

    STOP_TIMER("updating_bounding_sphere");


    //--- Second perform broad phase collision detection --------------------
    typedef          std::vector< std::pair<body_type*,body_type*> >  overlaps_container;
    typedef typename overlaps_container::iterator                     overlap_iterator;

    START_TIMER("broad_phase");

    std::vector< std::pair<body_type*,body_type*> > overlaps;
    T efficiency = VT::zero();

    broad_phase( engine, overlaps, efficiency);

    RECORD("broad_efficiency", efficiency);

    STOP_TIMER("broad_phase");

    //--- Third we perform narrow phase collision detction ---------------------
    START_TIMER("narrow_phase");

    engine.contacts().clear();

    std::vector< callback_type > callbacks;  // 2014-10-19 Kenny: Argh, I hate this design choice.... really ugly
    callbacks.resize( overlaps.size() );
    typename std::vector< callback_type >::iterator callback = callbacks.begin();


    std::vector< kdop::TestPair<V, 8, T> > narrow_work_pool;
    narrow_work_pool.resize( overlaps.size() );

    typename std::vector< kdop::TestPair<V, 8, T>  >::iterator work_item = narrow_work_pool.begin();

    for(overlap_iterator o = overlaps.begin(); o != overlaps.end(); ++o, ++callback, ++work_item)
    {
      body_type * bodyA = static_cast<body_type *>(o->first);    // 2009-11-25 Kenny: hmm can we not get rid of static casts?
      body_type * bodyB = static_cast<body_type *>(o->second);   // 2009-11-25 Kenny: hmm can we not get rid of static casts?

      *callback = callback_type( bodyA, bodyB, engine.contacts() );

      *work_item = kdop::TestPair<V, 8, T>(
                                  bodyA->m_tree
                                  , bodyB->m_tree
                                  , bodyA->m_mesh
                                  , bodyB->m_mesh
                                  , bodyA->m_X
                                  , bodyB->m_X
                                  , bodyA->m_Y
                                  , bodyB->m_Y
                                  , bodyA->m_Z
                                  , bodyB->m_Z
                                  , bodyA->m_opposite_face_on_surf
                                  , bodyB->m_opposite_face_on_surf
                                  , (*callback)
                                  );
    }

    narrow_phase( engine, narrow_work_pool );

    STOP_TIMER("narrow_phase");

    RECORD("contacts", engine.contacts().size() );
    
    STOP_TIMER("collision_detection");    
  }
} // namespace hyper

// HYPER_COLLISION_DETECTION_H
#endif
