#ifndef PROX_COMPUTE_RAYCAST_H
#define PROX_COMPUTE_RAYCAST_H

#include <prox_rigid_body.h>
#include <geometry.h>
#include <narrow.h>

#include <vector>
#include <cassert>

namespace prox
{

template < typename M>
inline bool compute_raycast(geometry::Ray<typename M::vector3_type> const& ray,
                            std::vector< RigidBody< typename M::real_type > >& bodies,
                            narrow::System<typename M::real_type>& narrow_system,
                            size_t& body_index, typename M::vector3_type& point,
                            typename M::real_type& distance)
{
    typedef typename M::tiny_types                                  tiny_types;
    typedef typename M::vector3_type                                V;
    typedef typename M::value_traits                                VT;
    typedef typename M::real_type                                   T;

    //--- Initialization -------------------------------------------------------
    point        = V::zero();
    distance     = std::numeric_limits<T>::max();
    bool did_hit = false;

    //--- Preprocessing geometries so they reflect current state ---------------
    std::vector< narrow::KDopBvhUpdateWorkItem< tiny_types > > kdop_bvh_update_work_pool;

    kdop_bvh_update_work_pool.reserve( bodies.size() ); // Make sure all space we may need is pre-allocated.

    for (auto body = bodies.begin(); body != bodies.end(); ++body)
    {
        auto const& geometry = narrow_system.get_geometry(body->get_geometry_idx());

        if (geometry.m_tetramesh.has_data())
        {
            narrow::KDopBvhUpdateWorkItem<tiny_types> work_item
                = narrow::KDopBvhUpdateWorkItem<tiny_types>(*body, geometry, body->get_position(),
                                                            body->get_orientation());
            kdop_bvh_update_work_pool.push_back(work_item);
        }
    }

    if( ! kdop_bvh_update_work_pool.empty() )
    {
      narrow::update_kdop_bvh(
                              kdop_bvh_update_work_pool
                              , narrow::sequential()
                              );
    }

    for (auto geometry = narrow_system.begin(); geometry != narrow_system.end(); ++geometry)
    {
      geometry->update_radius();
    }

    for (auto body = bodies.begin(); body != bodies.end(); ++body)
    {
        auto const& geometry = narrow_system.get_geometry(body->get_geometry_idx());

        body->set_radius(geometry.get_radius());
    }

    //--- Test each body for intersection and find the "closest" one -----------
    size_t idx = 0u;
    for (auto body = bodies.begin(); body != bodies.end(); ++body, ++idx)
    {
      T min_x = 0;
      T min_y = 0;
      T min_z = 0;
      T max_x = 0;
      T max_y = 0;
      T max_z = 0;

      body->get_box( min_x, min_y, min_z, max_x, max_y, max_z );

      const EigenVector3<T> min_coord = EigenVector3<T>( min_x, min_y, min_z );
      const EigenVector3<T> max_coord = EigenVector3<T>( max_x, max_y, max_z );
      geometry::AABBEigen<T> const aabb      = geometry::make_aabb( min_coord, max_coord );

      T       aabb_distance  = std::numeric_limits<T>::max();
      EigenVector3<T> aabb_point     = EigenVector3<T>(0,0,0);

      if ( geometry::compute_raycast_aabb<T>( geometry::convertRayToEigen<V>(ray), (aabb), aabb_point, aabb_distance) )
      {
        T       body_distance  = std::numeric_limits<T>::max();
        V       body_point     = V::zero();

        auto const& geometry = narrow_system.get_geometry(body->get_geometry_idx());

        bool body_hit = narrow::raycast(
                                   ray
                                   , *body
                                   , geometry
                                   , body->get_position()
                                   , body->get_orientation()
                                   , body_point
                                   , body_distance
                                   );

        //if( aabb_distance < distance )
        if( body_hit && (body_distance < distance) )
        {
          body_index    = idx;
          //point    = aabb_point;
          //distance = aabb_distance;
          point    = body_point;
          distance = body_distance;
          did_hit  = true;
        }

      }

    }

    return did_hit;
}

} // namespace prox

// PROX_COMPUTE_RAYCAST_H
#endif
