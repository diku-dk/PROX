#ifndef GEOMETRY_COMPUTE_RAYCAST_TETRAHEDRON_H
#define GEOMETRY_COMPUTE_RAYCAST_TETRAHEDRON_H

#include <types/geometry_tetrahedron.h>
#include <types/geometry_triangle.h>
#include <types/geometry_ray.h>

#include <cmath>
#include <vector>

namespace geometry
{

/**
   * Specialized version of raycast of tetrahedron.
   *
   * Here it is assumed that not all vertices of the
   * tetrahedron are on the surface. Hence only the parts
   * of the tetrahedron made up by "surface" triangles can
   * be hit by the tetrahedron.
   *
   * @param surface_map     An array with a boolean flag for each
   *                        vertex of the tetrahedron. If a flag is
   *                        true then the corresponding vertex is
   *                        a surface vertex.
   */
template<typename T>
inline bool compute_raycast_tetrahedron(RayEigen<T> const & ray
                                        , TetrahedronEigen<T> const & tetrahedron
                                        , EigenVector3<T>& hit
                                        , T & length
                                        , std::vector<bool> const & surface_map
                                        )
{
    assert(surface_map.size() == 4u || !"compute_raycast_tetrahedron(): surface map must have four elements");


    length = std::numeric_limits<T>::max();
    hit    = EigenVector3<T>(0,0,0);

    bool status = false;

    std::vector<EigenVector3<T>>  n(4u, EigenVector3<T>(0,0,0) );
    std::vector<T>  w(4u, 0 );

    for(unsigned int m = 0u; m < 4u; ++m)
    {
        Triangle<T> const triangle = get_opposite_face( m, tetrahedron );

        const EigenVector3<T> n0    = ( triangle.p(1)- triangle.p(0)).cross( triangle.p(2) - triangle.p(0));
        const EigenVector3<T> n1    = ( triangle.p(2)- triangle.p(1)).cross( triangle.p(0) - triangle.p(1));
        const EigenVector3<T> n2    = ( triangle.p(0)- triangle.p(2)).cross( triangle.p(1) - triangle.p(2));

        T const w0    = dot( n0, triangle.p(0));
        T const w1    = dot( n1, triangle.p(1));
        T const w2    = dot( n2, triangle.p(2));

        T const test0 = dot(n0, n0);
        T const test1 = dot(n1, n1);
        T const test2 = dot(n2, n2);

        if( test0 >= test1 && test0 >= test2)
        {
            n[m] = n0;
            w[m] = w0;
        }
        else if( test1 >= test0 && test1 >= test2)
        {
            n[m] = n1;
            w[m] = w1;
        }
        else  if( test2 >= test0 && test2 >= test1)
        {
            n[m] = n2;
            w[m] = w2;
        }
        else
        {
            assert(false || !"compute_raycast_tetrahedron() : could not determine best plane wrt. numerical precision");
        }
    }

    for(unsigned int m = 0u; m < 4u; ++m)
    {
        bool const & is_surface = surface_map[m];

        if(!is_surface)
            continue;

        T const rn  = dot((ray.direction()), n[m] );
        T const on  = dot((ray.origin()),    n[m] );

        if(rn >= 0)  // Ray is parallel with plane or hitting form back-side
            continue;

        T const t = (w[m] - on) / rn;   // compute ray length for hitting the plane

        if( t <  0 )  // if we hit behind the ray origin we give up
            continue;

        const EigenVector3<T> q = (ray.origin() + t* ray.direction());   // we know we are hitting the ray in front of ray origin

        unsigned int const i = (m+1u) % 4;
        unsigned int const j = (m+2u) % 4;
        unsigned int const k = (m+3u) % 4;

        T const tst_i = dot( n[ i ], q) - w[ i ];
        T const tst_j = dot( n[ j ], q) - w[ j ];
        T const tst_k = dot( n[ k ], q) - w[ k ];

        if( tst_i > 0 )  // Hit poit in front of face i
            continue;

        if( tst_j > 0 ) // Hit poit in front of face j
            continue;

        if( tst_k > 0 ) // Hit poit in front of face k
            continue;

        if(t < length)
        {
            // We found a hit-point on the surface of the tetrahedron
            // which is better than any previous hit point
            length = t;
            hit    = (q);
            status = true;
        }

    }

    return status;
}


template<typename T>
inline bool compute_raycast_tetrahedron(
    RayEigen<T> const & ray
    , TetrahedronEigen<T> const & tetrahedron
    , EigenVector3<T>& hit
    , T & length
    )
{
    std::vector<bool> surface_map(4u,true);

    return compute_raycast_tetrahedron(ray, tetrahedron, hit, length, surface_map);
}

}// namespace geometry

// GEOMETRY_COMPUTE_RAYCAST_TETRAHEDRON_H
#endif
