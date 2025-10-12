#ifndef PROX_COLLISION_DETECTION_H
#define PROX_COLLISION_DETECTION_H

#include <broad.h>
#include <broad_statistics.h>

#ifdef HAS_DIKUCL
#include <cl/narrow_cl_update_kdop_bvh.h>
#endif // HAS_DIKUCL

#include <geometry.h>

#include <narrow.h>
#include <narrow_tags.h>
#include <narrow_update_kdop_bvh.h>

#include <prox_contact_point.h>
#include <prox_params.h>
#include <prox_rigid_body.h>

#include <util_profiling.h>

#include <cassert>
#include <vector>
#include <type_traits>

namespace prox
{

namespace detail
{

/**
     * Contact Point Generation Callback Functor.
     * This functor acts as the glue between the contact point generation algorithm and the Prox contact point data types and containers.
     *
     * The functor hides the nasty details and data types from the contact point generation library.
     */
template <typename T>
class ContactCallbackFunctor : public geometry::ContactsCallback<T>
{
public:
    typedef prox::RigidBody<T> body_type;
    typedef prox::ContactPoint<T> contact_type;

protected:
    body_type* m_body_i;   ///< A pointer to body i of the contact.
    body_type* m_body_j;   ///< A pointer to body j of the contact.
    std::vector< contact_type >*
        m_results; ///< A pointer to a contact point container where all generated contacts should be added to.

public:
    ContactCallbackFunctor()
        : m_body_i(0)
        , m_body_j(0)
        , m_results(0)
    {
    }

    ContactCallbackFunctor(body_type* A, body_type* B,
                           std::vector<prox::ContactPoint<T>>& results)
        : m_body_i(A)
        , m_body_j(B)
        , m_results(&results)
    {
        assert(A || !"ContactCallbackFunctor(...) body A was null");
        assert(B || !"ContactCallbackFunctor(...) body B was null");
    }

    ~ContactCallbackFunctor() {}

    ContactCallbackFunctor& operator=(ContactCallbackFunctor const& callback)
    {
        if (this != &callback)
        {
            this->m_body_i = callback.m_body_i;
            this->m_body_j = callback.m_body_j;
            this->m_results = callback.m_results;
        }
        return *this;
    }

    /*
       * Add new contact callback function.
       * The contact point generation library invokes this function call whenever it has discovered a new contact that should be stored.
       *
       * @param p    The position of the new contact point.
       * @param n    The unit normal of the new contact point.
       * @param d    The penetration distance meassure of the contact point.
       */
    void tempParenthesisOperatorImpl(const EigenVector3<T>& point, const EigenVector3<T>& normal,
                                     const T& distance)
    {
        assert(this->m_body_i || !"operator(): body i is null");
        assert(this->m_body_j || !"operator(): body j is null");
        assert(this->m_results || !"operator(): results is null");

        contact_type contact;
        contact.position = point;
        contact.depth = distance;
        contact.normal = normal;
        contact.bodyI = m_body_i;
        contact.bodyJ = m_body_j;

        this->m_results->push_back(contact);
    }
};

} // namespace detail

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
template < typename T>
requires std::is_floating_point_v<T>
inline void collision_detection(std::vector< RigidBody<T>>& bodies, broad::System<T>& broad_system,
                                narrow::System<T>& narrow_system,
                                std::vector< ContactPoint<T> >& contacts, Params<T> const& params)
{
    typedef typename broad::System<T>::overlap_type overlap_type;
    typedef std::vector<overlap_type> overlap_container;

    START_TIMER("collision_detection");

    //--- First we preprocess data structures for doing collision detection ----
    {
        START_TIMER("collision_detection_preprocessing");

        //--- Update kDOP BVHs to reflect changes in tetramesh geometry ----------
        START_TIMER("collision_detection_creating_kdop_work_pool");

        std::vector<narrow::KDopBvhUpdateWorkItem<T>> kdop_bvh_update_work_pool;

        kdop_bvh_update_work_pool.reserve(
            bodies.size()); // Make sure all space we may need is pre-allocated.

        for (auto body = bodies.begin(); body != bodies.end(); ++body)
        {
            const auto& geometry = narrow_system.get_geometry(body->get_geometry_idx());

            if (geometry.m_tetramesh.has_data())
            {
                auto work_item = narrow::KDopBvhUpdateWorkItem<T>(
                    *body, geometry, body->get_position(),
                    body->get_orientation());
                kdop_bvh_update_work_pool.push_back(work_item);
            }
            /*if (geometry.m_signedDistanceMap.hasData())
            {
                CoordSysEigen<T> data(body->get_position(),
                                      body->get_orientation());
                geometry.m_signedDistanceMap.setSignedDistanceTransform(data);
            }*/
        }
        STOP_TIMER("collision_detection_creating_kdop_work_pool");

        START_TIMER("collision_detection_updating_kdop");
        if (!kdop_bvh_update_work_pool.empty())
        {
#ifdef HAS_DIKUCL
            if (narrow_system.params().use_open_cl())
            {
                narrow::update_kdop_bvh(kdop_bvh_update_work_pool, narrow::dikucl(),
                                        narrow_system.params().open_cl_platform(),
                                        narrow_system.params().open_cl_device());
            }
            else
            {
#endif // HAS_DIKUCL

                // use regular updating of KDOP BVHs if DIKUCL is not available or should not be used
                narrow::update_kdop_bvh(kdop_bvh_update_work_pool, narrow::sequential());

#ifdef HAS_DIKUCL
            }
#endif // HAS_DIKUCL
        }
        STOP_TIMER("collision_detection_updating_kdop");

        //--- Update bounding spheres (radius) of all geometries in the system -----
        for (auto geometry = narrow_system.begin(); geometry != narrow_system.end(); ++geometry)
        {
            // 2013-06-07 Kenny code review: Ideally only geometry
            //            that has changed should have the radius updated.
            //            Note that tetrameshes are special as only their
            //            material space mesh are shared between the
            //            associated bodies.

            geometry->update_radius();
        }

        // 2013-07-06 Kenny code review: This will clear all data in
        //            the broad phase collision detection system. This
        //            may not be the most efficient approach, as all
        //            configuration data needs to be build up from stracth
        //            again. Ideally only newly added bodies or geometry
        //            changed bodies should be cleared/added to the broad
        //            phase system.
        broad_system.clear();

        //--- Update body radius and connect their AABBs to the broad phase system -
        for (auto body = bodies.begin(); body != bodies.end(); ++body)
        {
            auto const& geometry = narrow_system.get_geometry(body->get_geometry_idx());

            body->set_radius(geometry.get_radius());

            // 2013-07-06 Kenny code review: Here we re-connect all bodies
            //            to the broad phase collision detection system.
            //            See my review comment above about efficiency.
            broad_system.connect(&(*body));
        }

        // 2015-03-03 Kenny code review: This optimal spacing requires
        // sorting of all objects, so it runs O(n lg n). However, grid
        // algorithm is trying to run in O(n). Hence, one could argue that
        // a sweep-line algorithm would be better as its performance do
        // not depend on the obejct sizes.
        broad_system.compute_optimal_cell_spacing();

        STOP_TIMER("collision_detection_preprocessing");
    }

    //--- Second we perform broad phase collision detection --------------------
    overlap_container overlaps;
    {
        START_TIMER("broad_phase");

        // We use the broad phase collision detection system to quickly determine pairs of overlapping bodies

        float efficiency = 0.0f;

        if (params.use_all_pair())
        {
            broad::find_overlaps(broad_system, overlaps, efficiency, broad::all_pair_algorithm());
        }
        else { broad::find_overlaps(broad_system, overlaps, efficiency, broad::grid_algorithm()); }

        STOP_TIMER("broad_phase");

#ifdef USE_PROFILING
        RECORD("broad_efficiency", efficiency);
        broad::Statistics stats = broad::make_statistics(broad_system);
        RECORD("broad_fill_ratio", stats.m_fill_ratio);
        RECORD("min_obj", stats.m_min_obj);
        RECORD("max_obj", stats.m_max_obj);
        RECORD("avg_obj", stats.m_avg_obj);
        RECORD("std_obj", stats.m_std_obj);
#endif
    }

    //--- Third we perform narrow phase collision detction ---------------------
    {
        START_TIMER("narrow_phase");

        //--- Make sure we do not carry any old contact information around. ------
        contacts.clear();

        std::vector<narrow::TestPair<T>> narrow_test_pairs;

        typedef detail::ContactCallbackFunctor<T> callback_type;

        std::vector<callback_type>
            callbacks; // 2014-10-19 Kenny: Argh, I hate this design choice.... really ugly
        callbacks.resize(overlaps.size());
        auto callback = callbacks.begin();

        for (auto o = overlaps.begin(); o != overlaps.end(); ++o, ++callback)
        {
            auto* bodyA = dynamic_cast<RigidBody<T>*>(o->first);
            auto* bodyB = dynamic_cast<RigidBody<T>*>(o->second);

            //--- Verify if we need to test the two bodies or if we can skip them --
            if (bodyA->is_fixed() && bodyB->is_fixed()) continue;
            if (bodyA->is_fixed() && bodyB->is_scripted()) continue;
            if (bodyA->is_scripted() && bodyB->is_fixed()) continue;
            if (bodyA->is_scripted() && bodyB->is_scripted()) continue;

            *callback = callback_type(bodyA, bodyB, contacts);

            narrow::TestPair<T> narrow_pair(
                *bodyA, *bodyB, (bodyA->get_position()),
                (bodyA->get_orientation()), (bodyB->get_position()),
                (bodyB->get_orientation()), (*callback));

            narrow_test_pairs.push_back(narrow_pair);

            // 2015-02-01 Kenny code review: Disptaching one-by-one is old-style
            // fashion way of doing this. The batch dispatch is nicer as it allows
            // narrow phase system to exploit parallism. The one-by-one behavior
            // forces sequential behavior.
            if (!narrow_system.params().use_batching())
            {
                narrow::dispatch_collision_handlers(narrow_system, narrow_test_pairs);
                narrow_test_pairs.clear();
            }
        }

        if (!narrow_test_pairs.empty())
        {
            narrow::dispatch_collision_handlers(narrow_system, narrow_test_pairs);
        }

        STOP_TIMER("narrow_phase");
    }

    //--- Forth phase we perform contact point reduction -----------------------
    if (params.stepper_params().contact_reduction())
    {
        START_TIMER("contact_reduction");

        std::sort(contacts.begin(), contacts.end());
        auto last = std::unique(contacts.begin(), contacts.end());
        contacts.erase(last, contacts.end());

        STOP_TIMER("contact_reduction");
    }

    //--- Fifth phase we collect statistics on contact points ------------------
    {
        RECORD("contacts", contacts.size());

        util::Profiling::Monitor* penetration_monitor = util::Profiling::get_monitor("penetration");

        penetration_monitor->clear();

        for (const auto& contact : contacts) { RECORD("penetration", contact.depth); }

        RECORD("average_penetration", penetration_monitor->get_avg());
        RECORD("minimum_penetration", penetration_monitor->get_min());
        RECORD("maximum_penetration", penetration_monitor->get_max());

        penetration_monitor->clear();
    }

    STOP_TIMER("collision_detection");
}
} // namespace prox

namespace detail
{

/**
     * Contact Point Generation Callback Functor.
     * This functor acts as the glue between the contact point generation algorithm and the Prox contact point data types and containers.
     *
     * The functor hides the nasty details and data types from the contact point generation library.
     */
template <typename T>
class ContactCallbackFunctor : public geometry::ContactsCallback<T>
{
public:
    typedef prox::RigidBody<T> body_type;
    typedef prox::ContactPoint<T> contact_type;

protected:
    body_type* m_body_i; ///< A pointer to body i of the contact.
    body_type* m_body_j; ///< A pointer to body j of the contact.
    std::vector<contact_type>*
        m_results; ///< A pointer to a contact point container where all generated contacts should be added to.

public:
    ContactCallbackFunctor()
        : m_body_i(0)
        , m_body_j(0)
        , m_results(0)
    {
    }

    ContactCallbackFunctor(body_type* A, body_type* B,
                           std::vector<prox::ContactPoint<T>>& results)
        : m_body_i(A)
        , m_body_j(B)
        , m_results(&results)
    {
        assert(A || !"ContactCallbackFunctor(...) body A was null");
        assert(B || !"ContactCallbackFunctor(...) body B was null");
    }

    ~ContactCallbackFunctor() {}

    ContactCallbackFunctor& operator=(ContactCallbackFunctor const& callback)
    {
        if (this != &callback)
        {
            this->m_body_i = callback.m_body_i;
            this->m_body_j = callback.m_body_j;
            this->m_results = callback.m_results;
        }
        return *this;
    }

    /*
       * Add new contact callback function.
       * The contact point generation library invokes this function call whenever it has discovered a new contact that should be stored.
       *
       * @param p    The position of the new contact point.
       * @param n    The unit normal of the new contact point.
       * @param d    The penetration distance meassure of the contact point.
       */
    void tempParenthesisOperatorImpl(const EigenVector3<T>& point,
                                     const EigenVector3<T>& normal,
                                     const T& distance)
    {
        assert(this->m_body_i || !"operator(): body i is null");
        assert(this->m_body_j || !"operator(): body j is null");
        assert(this->m_results || !"operator(): results is null");

        contact_type contact;
        contact.position = point;
        contact.depth = distance;
        contact.normal = normal;
        contact.bodyI = m_body_i;
        contact.bodyJ = m_body_j;

        this->m_results->push_back(contact);
    }
};

} // namespace detail

namespace prox
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
template <typename T>
requires std::is_floating_point_v<T>
inline T collision_detection_CCD(std::vector<RigidBody<T>>& bodies,
                                 narrow::System<T>& narrow_system,
                                 std::vector<ContactPoint<T>>& contacts,
                                 T startTime, T endTime)
{
    typedef typename broad::System<T>::overlap_type overlap_type;
    typedef std::vector<overlap_type> overlap_container;

    START_TIMER("continuous_collision_detection");

    const std::size_t n = bodies.size();
    std::vector<narrow::TestPairCCD<T>> narrow_test_pairs;

    //Note all this callback mess is not needed: I only use it for quick compatibility, this will be changed!
    typedef detail::ContactCallbackFunctor<T> callback_type;

    std::vector<callback_type> callbacks;
    callbacks.resize(n * n);
    auto callback = callbacks.begin();
    ;
    contacts.clear();
    std::vector<kdop::BodyVelocities<T>> bodyVels;
    T earliestTOI = std::numeric_limits<T>::max();
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            if (i == j) { continue; }
            RigidBody<T>* pair_body_A = &bodies[i];
            RigidBody<T>* pair_body_B = &bodies[j];

            auto* bodyA = dynamic_cast<RigidBody<T>*>(pair_body_A);
            auto* bodyB = dynamic_cast<RigidBody<T>*>(pair_body_B);

            //--- Verify if we need to test the two bodies or if we can skip them --
            if (bodyA->is_fixed() && bodyB->is_fixed()) continue;
            if (bodyA->is_fixed() && bodyB->is_scripted()) continue;
            if (bodyA->is_scripted() && bodyB->is_fixed()) continue;
            if (bodyA->is_scripted() && bodyB->is_scripted()) continue;

            //*callback = callback_type(bodyA, bodyB, contacts);

            narrow::TestPairCCD<T> narrow_pair(
                *bodyA, *bodyB, (bodyA->get_position()),
                (bodyA->get_orientation()), (bodyB->get_position()),
                (bodyB->get_orientation()));

            kdop::BodyVelocities<T> body{
                .bodyALinVel = &pair_body_A->get_velocity(),
                .bodyAAngVel = &pair_body_A->get_spin(),
                .bodyACenterTranslation = &pair_body_A->get_position(),
                .bodyACenterRotation = &pair_body_A->get_orientation(),
                .bodyBLinVel = &pair_body_B->get_velocity(),
                .bodyBAngVel = &pair_body_B->get_spin(),
                .bodyBCenterTranslation = &pair_body_B->get_position(),
                .bodyBCenterRotation = &pair_body_B->get_orientation(),
            };

            narrow_test_pairs.push_back(narrow_pair);
            bodyVels.push_back(body);
            if (!narrow_system.params().use_batching())
            {
                T toi = narrow::dispatch_collision_handlers_CCD(
                    narrow_system, narrow_test_pairs, startTime, endTime,
                    bodyVels);
                earliestTOI = std::min<T>(toi, earliestTOI);
                narrow_test_pairs.clear();
                bodyVels.clear();
            }
        }
    }
    if (narrow_system.params().use_batching())
    {
        T toi = narrow::dispatch_collision_handlers_CCD(
            narrow_system, narrow_test_pairs, startTime, endTime, bodyVels);
        earliestTOI = std::min<T>(toi, earliestTOI);
        narrow_test_pairs.clear();
        bodyVels.clear();
    }
    std::cerr << "WE HAVE A TOI OF " << earliestTOI << "\n";
    STOP_TIMER("continuous_collision_detection");
    return earliestTOI;
}

/*
    //--- First we preprocess data structures for doing collision detection ----
    {
        START_TIMER("continuous_collision_detection_preprocessing");

        //--- Update kDOP BVHs to reflect changes in tetramesh geometry ----------
        START_TIMER("continuous_collision_detection_creating_kdop_work_pool");

        std::vector<narrow::KDopBvhUpdateWorkItem<T>> kdop_bvh_update_work_pool;

        kdop_bvh_update_work_pool.reserve(
            bodies.size()); // Make sure all space we may need is pre-allocated.

        for (auto body = bodies.begin(); body != bodies.end(); ++body)
        {
            const auto& geometry
                = narrow_system.get_geometry(body->get_geometry_idx());

            if (geometry.m_tetramesh.has_data())
            {
                auto work_item = narrow::KDopBvhUpdateWorkItem<T>(
                    *body, geometry, body->get_position(),
                    body->get_orientation());
                kdop_bvh_update_work_pool.push_back(work_item);
            }

        }
        STOP_TIMER("continuous_collision_detection_creating_kdop_work_pool");

        START_TIMER("continuous_collision_detection_updating_kdop");
        if (!kdop_bvh_update_work_pool.empty())
        {
#ifdef HAS_DIKUCL
            if (narrow_system.params().use_open_cl())
            {
                narrow::update_kdop_bvh(
                    kdop_bvh_update_work_pool, narrow::dikucl(),
                    narrow_system.params().open_cl_platform(),
                    narrow_system.params().open_cl_device());
            }
            else
            {
#endif // HAS_DIKUCL

                // use regular updating of KDOP BVHs if DIKUCL is not available or should not be used
                narrow::update_kdop_bvh(kdop_bvh_update_work_pool,
                                        narrow::sequential());

#ifdef HAS_DIKUCL
            }
#endif // HAS_DIKUCL
        }
        STOP_TIMER("continuous_collision_detection_updating_kdop");

        //--- Update bounding spheres (radius) of all geometries in the system -----
        for (auto geometry = narrow_system.begin();
             geometry != narrow_system.end(); ++geometry)
        {
            // 2013-06-07 Kenny code review: Ideally only geometry
            //            that has changed should have the radius updated.
            //            Note that tetrameshes are special as only their
            //            material space mesh are shared between the
            //            associated bodies.

            geometry->update_radius();
        }

        // 2013-07-06 Kenny code review: This will clear all data in
        //            the broad phase collision detection system. This
        //            may not be the most efficient approach, as all
        //            configuration data needs to be build up from stracth
        //            again. Ideally only newly added bodies or geometry
        //            changed bodies should be cleared/added to the broad
        //            phase system.
        broad_system.clear();

        //--- Update body
 and connect their AABBs to the broad phase system -
        for (auto body = bodies.begin(); body != bodies.end(); ++body)
        {
            auto const& geometry
                = narrow_system.get_geometry(body->get_geometry_idx());

            body->set_radius(geometry.get_radius());

            // 2013-07-06 Kenny code review: Here we re-connect all bodies
            //            to the broad phase collision detection system.
            //            See my review comment above about efficiency.
            broad_system.connect(&(*body));
        }

        // 2015-03-03 Kenny code review: This optimal spacing requires
        // sorting of all objects, so it runs O(n lg n). However, grid
        // algorithm is trying to run in O(n). Hence, one could argue that
        // a sweep-line algorithm would be better as its performance do
        // not depend on the obejct sizes.
        broad_system.compute_optimal_cell_spacing();

        STOP_TIMER("continuous_collision_detection_preprocessing");
    }

    //--- Second we perform broad phase collision detection --------------------
    overlap_container overlaps;
    {
        START_TIMER("broad_phase");

        // We use the broad phase collision detection system to quickly determine pairs of overlapping bodies

        float efficiency = 0.0f;

        if (params.use_all_pair())
        {
            broad::find_overlaps(broad_system, overlaps, efficiency,
                                 broad::all_pair_algorithm());
        }
        else
        {
            broad::find_overlaps(broad_system, overlaps, efficiency,
                                 broad::grid_algorithm());
        }

        STOP_TIMER("broad_phase");

#ifdef USE_PROFILING
        RECORD("broad_efficiency", efficiency);
        broad::Statistics stats = broad::make_statistics(broad_system);
        RECORD("broad_fill_ratio", stats.m_fill_ratio);
        RECORD("min_obj", stats.m_min_obj);
        RECORD("max_obj", stats.m_max_obj);
        RECORD("avg_obj", stats.m_avg_obj);
        RECORD("std_obj", stats.m_std_obj);
#endif
    }

    //--- Third we perform narrow phase collision detction ---------------------
    {
        START_TIMER("narrow_phase");

        //--- Make sure we do not carry any old contact information around. ------
        contacts.clear();

        std::vector<narrow::TestPair<T>> narrow_test_pairs;
        typedef detail::ContactCallbackFunctor<T> callback_type;

        std::vector<callback_type>
            callbacks; // 2014-10-19 Kenny: Argh, I hate this design choice.... really ugly
        callbacks.resize(overlaps.size());
        auto callback = callbacks.begin();

        for (auto o = overlaps.begin(); o != overlaps.end(); ++o, ++callback)
        {
            auto* bodyA = dynamic_cast<RigidBody<T>*>(o->first);
            auto* bodyB = dynamic_cast<RigidBody<T>*>(o->second);

            //--- Verify if we need to test the two bodies or if we can skip them --
            if (bodyA->is_fixed() && bodyB->is_fixed()) continue;
            if (bodyA->is_fixed() && bodyB->is_scripted()) continue;
            if (bodyA->is_scripted() && bodyB->is_fixed()) continue;
            if (bodyA->is_scripted() && bodyB->is_scripted()) continue;

            *callback = callback_type(bodyA, bodyB, contacts);

            narrow::TestPair<T> narrow_pair(
                *bodyA, *bodyB, (bodyA->get_position()),
                (bodyA->get_orientation()), (bodyB->get_position()),
                (bodyB->get_orientation()), (*callback));

            narrow_test_pairs.push_back(narrow_pair);

            // 2015-02-01 Kenny code review: Disptaching one-by-one is old-style
            // fashion way of doing this. The batch dispatch is nicer as it allows
            // narrow phase system to exploit parallism. The one-by-one behavior
            // forces sequential behavior.
            if (!narrow_system.params().use_batching())
            {
                narrow::dispatch_collision_handlers(narrow_system,
                                                    narrow_test_pairs);
                narrow_test_pairs.clear();
            }
        }

        if (!narrow_test_pairs.empty())
        {
            narrow::dispatch_collision_handlers(narrow_system,
                                                narrow_test_pairs);
        }

        STOP_TIMER("narrow_phase");
    }

    //--- Forth phase we perform contact point reduction -----------------------
    if (params.stepper_params().contact_reduction())
    {
        START_TIMER("contact_reduction");

        std::sort(contacts.begin(), contacts.end());
        auto last = std::unique(contacts.begin(), contacts.end());
        contacts.erase(last, contacts.end());

        STOP_TIMER("contact_reduction");
    }

    //--- Fifth phase we collect statistics on contact points ------------------
    {
        RECORD("contacts", contacts.size());

        util::Profiling::Monitor* penetration_monitor
            = util::Profiling::get_monitor("penetration");

        penetration_monitor->clear();

        for (const auto& contact : contacts)
        {
            RECORD("penetration", contact.depth);
        }

        RECORD("average_penetration", penetration_monitor->get_avg());
        RECORD("minimum_penetration", penetration_monitor->get_min());
        RECORD("maximum_penetration", penetration_monitor->get_max());

        penetration_monitor->clear();
    }*/

} // namespace prox

// PROX_COLLISION_DETECTION_H
#endif
