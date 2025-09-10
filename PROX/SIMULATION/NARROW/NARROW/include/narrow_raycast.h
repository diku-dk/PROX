#ifndef NARROW_RAYCAST_H
#define NARROW_RAYCAST_H

#include <narrow_object.h>
#include <narrow_geometry.h>

#include <eigenhelperall.h>
#include <geometry.h>

#include <kdop_raycast.h>

namespace narrow
{

template <typename T>
inline bool raycast(geometry::RayEigen<T> const& ray, Object<T> const& objA,
                    Geometry<T> const& geoA, EigenVector3<T> const& tA,
                    EigenQuaternion<T> const& qA, EigenVector3<T>& point,
                    T& distance)
{
    point = EigenVector3<T>{0, 0, 0};
    distance = std::numeric_limits<T>::max();

    if (geoA.m_tetramesh.has_data())
    {
        return kdop::raycast<T, 8>(
            ray, objA.m_tree, geoA.m_tetramesh.m_mesh, objA.m_X, objA.m_Y,
            objA.m_Z, geoA.m_tetramesh.m_surface_map, point, distance);
    }
    else
    {
        auto const bodyAtoWCS = CoordSysEigen<T>(tA, qA);

        if (!geoA.m_boxes.empty())
        {
            for (auto a = geoA.m_boxes.begin(); a != geoA.m_boxes.end(); ++a)
            {
                const CoordSysEigen<T> shapeAtobodyA = CoordSysEigen<T>(
                    (a->transform().T()), (a->transform().Q()));
                const CoordSysEigen<T> shapeAtoWCS
                    = prod(shapeAtobodyA, (bodyAtoWCS));

                geometry::OBBEigen<T> const obb = geometry::make_obb<T>(
                    shapeAtoWCS.T(), shapeAtoWCS.Q(), (a->half_extent()));

                T local_distance = std::numeric_limits<T>::max();
                EigenVector3<T> local_point = EigenVector3<T>(0, 0, 0);

                bool const did_hit = geometry::compute_raycast_obb(
                    (ray), obb, local_point, local_distance);

                distance = did_hit ? local_distance : distance;
                point = did_hit ? (local_point) : point;
            }
        }
        if (!geoA.m_spheres.empty())
        {
            for (auto a = geoA.m_spheres.begin(); a != geoA.m_spheres.end();
                 ++a)
            {
                auto const shapeAtobodyA = a->transform();
                auto const shapeAtoWCS = prod(shapeAtobodyA, bodyAtoWCS);

                geometry::Sphere<T> const sphere
                    = geometry::make_sphere<T>(shapeAtoWCS.T(), a->radius());

                T local_distance = std::numeric_limits<T>::max();
                EigenVector3<T> local_point = EigenVector3<T>(0, 0, 0);

                bool const did_hit = geometry::compute_raycast_sphere<T>(
                    (ray), sphere, local_point, local_distance);

                distance = did_hit ? local_distance : distance;
                point = did_hit ? (local_point) : point;
            }
        }
        if (!geoA.m_hulls.empty())
        {
            for (auto a = geoA.m_hulls.begin(); a != geoA.m_hulls.end(); ++a)
            {
                auto const shapeAtobodyA = a->transform();
                auto const shapeAtoWCS = prod(shapeAtobodyA, bodyAtoWCS);

                assert(false || !"not implemented yet");
                //          geometry::ConvexHull<V> const hull = geometry::make_convex_hull( shapeAtoWCS.T(), a->radius());

                T local_distance = std::numeric_limits<T>::max();
                EigenVector3<T> localPoint(0, 0, 0);

                bool const did_hit
                    = false; //geometry::compute_raycast_convex(ray, hull, local_point, local_distance);

                distance = did_hit ? local_distance : distance;
                point = did_hit ? localPoint : point;
            }
        }
    }

    return distance < std::numeric_limits<T>::max();
}

} //namespace narrow

// NARROW_RAYCAST_H
#endif
