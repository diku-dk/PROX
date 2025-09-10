#ifndef GEOMETRY_CLOSEST_POINTS_TETRAHEDRON_TETRAHEDRON_H
#define GEOMETRY_CLOSEST_POINTS_TETRAHEDRON_TETRAHEDRON_H

#include <types/geometry_tetrahedron.h>
#include <types/geometry_triangle.h>
#include <contacts/geometry_contacts_tetrahedron_tetrahedron.h>  // Needed for UnscaledPlane stuff

#include <eigenhelperall.h>

#include <vector>
#include <cassert>

namespace geometry
{
  template< typename T>
  inline void closest_points_tetrahedron_tetrahedron(
                                                     TetrahedronEigen<T> const & A
                                                     , TetrahedronEigen<T> const & B
                                                     , EigenVector3<T>& a
                                                     , EigenVector3<T>& b
                                                     , T& min_distance
                                                     )
  {

    unsigned int const edge[6][2] = {
      {0,  1}
      , {0,  2}
      , {0,  3}
      , {1,  2}
      , {1,  3}
      , {2,  3}
    };


    std::vector<Triangle<T> >                trianglesA(4u);
    std::vector<Triangle<T> >                trianglesB(4u);
    std::vector<details::UnscaledPlane<T> >  planesA(4u);
    std::vector<details::UnscaledPlane<T> >  planesB(4u);

    //--- Pre-computation ---- -------------------------------------------------
    for (unsigned int v =0u; v < 4u; ++v)
    {
      trianglesA[v] = get_opposite_face( v, A );
      trianglesB[v] = get_opposite_face( v, B );
      planesA[v]    = details::make_unscaled_plane_eigen(trianglesA[v]);
      planesB[v]    = details::make_unscaled_plane_eigen(trianglesB[v]);
    }

    min_distance = std::numeric_limits<T>::max();

    //--- Search for V-V cases -------------------------------------------------

    for (unsigned int i =0u; i < 4u; ++i)
    {
      const EigenVector3<T>& ai = A.p(i);

      for (unsigned int j =0u; j < 4u; ++j)
      {
        const EigenVector3<T>& bj = B.p(j);

          T const distance = (ai-bj).norm();

        a            = (distance < min_distance) ? ai : a;
        b            = (distance < min_distance) ? bj : b;
        min_distance =  distance < min_distance  ? distance : min_distance;
      }
    }

    //--- Search for V-E cases -------------------------------------------------
    for (unsigned int i =0u; i < 4u; ++i)
    {
      // Test vertex of A against edge from B
      for (unsigned int k =0u; k < 6u; ++k)
      {
        const EigenVector3<T>& p_a  = A.p(i);
        const EigenVector3<T>& p_b0 = B.p(edge[k][0]);
        const EigenVector3<T>& p_b1 = B.p(edge[k][1]);

        const EigenVector3<T> da    = p_a  - p_b0;
        const EigenVector3<T> db    = p_b1 - p_b0;

        T const t = (dot(da, db)/dot(db, db));

        if (t <= 0)
          continue;

        if (t >= 1)
          continue;

        const EigenVector3<T> p_b      = p_b0 + t * db;
        const EigenVector3<T> ortho    = p_a - p_b;
        T const distance = (ortho).norm();

        a            = (distance < min_distance) ? p_a      : a;
        b            = (distance < min_distance) ? p_b      : b;
        min_distance =  distance < min_distance  ? distance : min_distance;
      }

      // Test V of B against E of A
      for (unsigned int k =0u; k < 6u; ++k)
      {
        const EigenVector3<T> & p_b  = B.p(i);
        const EigenVector3<T>& p_a0 = A.p(edge[k][0]);
        const EigenVector3<T>& p_a1 = A.p(edge[k][1]);

        const EigenVector3<T> db    = p_b  - p_a0;
        const EigenVector3<T> da    = p_a1 - p_a0;

        T const t = (dot(da, db)/dot(da, da));

        if (t <= 0)
          continue;

        if (t >= 1)
          continue;

        const EigenVector3<T> p_a      = p_a0 + t * da;
        const EigenVector3<T> ortho    = p_a - p_b;
        T const distance = (ortho).norm();

        a            = (distance < min_distance) ? p_a      : a;
        b            = (distance < min_distance) ? p_b      : b;
        min_distance =  distance < min_distance  ? distance : min_distance;
      }

    }

    //--- Search for E-E cases -------------------------------------------------
    for (unsigned int k = 0u; k < 6u; ++k)
    {
      const EigenVector3<T>& ai = A.p(edge[k][0u]);
      const EigenVector3<T>& aj = A.p(edge[k][1u]);
      const  EigenVector3<T> da = (aj - ai).normalized();

      for (unsigned int m = 0u; m < 6u; ++m)
      {
        const EigenVector3<T>& bi = B.p(edge[m][0u]);
        const EigenVector3<T>& bj = B.p(edge[m][1u]);
        const EigenVector3<T> db = (bj - bi).normalized();

        const EigenVector3<T> r  =   bi - ai;
        T const k  =   dot(da, db);
        T const q1 =   dot(da,r);
        T const q2 = - dot(db,r);
        T const w  =   1 - k*k;

        // Test if edges too close to parallel
        if(fabs(w) < std::numeric_limits<T>::epsilon()*10)
          continue;

        // Compute edge-parameters corresponding to closest points
        T const t = (q1 + k*q2)/w;
        T const s = (q2 + k*q1)/w;

        // Test if closest points are interior on edges
        if (t<= 0)
          continue;

        if (s<= 0)
          continue;

        if (t>= (aj - ai).norm())
          continue;

        if (s>= (bj - bi).norm())
          continue;

        const EigenVector3<T> p_a      = ai + t*da;
        const EigenVector3<T> p_b      = bi + s*db;
        T const distance = (p_a - p_b).norm();

        a            = (distance < min_distance) ? p_a      : a;
        b            = (distance < min_distance) ? p_b      : b;
        min_distance =  distance < min_distance  ? distance : min_distance;
      }
    }

    //--- Search for V-F cases -------------------------------------------------
    for (unsigned int i = 0u; i < 4u; ++i)
    {
      // Vs of A against triangles of B
      for (unsigned int k = 0u; k < 4u; ++k)
      {
        const EigenVector3<T> p_a = A.p(i);

        bool inside_voronoi_planes = true;

        const EigenVector3<T>& n_b  = planesB[k].m_normal;

        for(unsigned int m = 0u; m < 3u;++m)
        {
          unsigned int const n = (m + 1) % 3;

          const EigenVector3<T>& p_b0 = trianglesB[k].p(m);
          const EigenVector3<T>& p_b1 = trianglesB[k].p(n);
          const  EigenVector3<T> e_b  = p_b1 - p_b0;

          const EigenVector3<T> vp_n  = (n_b).cross( e_b);
          T const tst  = ( p_a - p_b0).dot( vp_n);

          inside_voronoi_planes = tst < 0 ? false : inside_voronoi_planes;
        }

        if (!inside_voronoi_planes)
          continue;

        const EigenVector3<T>& q_b      = planesB[k].m_point;
        const  EigenVector3<T> p_b      =  p_a -  n_b * (n_b).dot( p_a - q_b)/dot(n_b, n_b);
        T const   distance = ( p_a - q_b).dot( n_b)/ ( n_b ).norm();

        if(distance < 0 )
          continue;

        a            = (distance < min_distance) ? p_a   : a;
        b            = (distance < min_distance) ? p_b   : b;
        min_distance =  distance < min_distance  ? distance : min_distance;
      }

      // Vs of B against triangles of A
      for (unsigned int k = 0u; k < 4u; ++k)
      {
        const EigenVector3<T> p_b = B.p(i);

        bool inside_voronoi_planes = true;

        const EigenVector3<T>& n_a  = planesA[k].m_normal;

        for(unsigned int m = 0u; m < 3u;++m)
        {
          unsigned int const n = (m + 1) % 3;

          const EigenVector3<T>& p_a0 = trianglesA[k].p(m);
          const EigenVector3<T>& p_a1 = trianglesA[k].p(n);
          const  EigenVector3<T> e_a  = p_a1 - p_a0;

          const EigenVector3<T> vp_n  = (n_a).cross( e_a);
          T const tst  = ( p_b - p_a0).dot( vp_n);

          inside_voronoi_planes = tst < 0 ? false : inside_voronoi_planes;
        }

        if (!inside_voronoi_planes)
          continue;

        const EigenVector3<T>& q_a      = planesA[k].m_point;
        const  EigenVector3<T> p_a      =  p_b -  n_a * (n_a).dot( p_b - q_a)/dot(n_a, n_a);
        T const   distance = ( p_b - q_a).dot( n_a)/ ( n_a ).norm();

        if(distance < 0 )
          continue;

        a            = (distance < min_distance) ? p_a   : a;
        b            = (distance < min_distance) ? p_b   : b;
        min_distance =  distance < min_distance  ? distance : min_distance;
      }
    }
  }

}// end namespace geometry

// GEOMETRY_CLOSEST_POINTS_TETRAHEDRON_TETRAHEDRON_H
#endif
