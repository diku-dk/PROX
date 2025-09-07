#ifndef GEOMETRY_CONTACTS_TETRAHEDRON_TETRAHEDRON_GROWTH_H
#define GEOMETRY_CONTACTS_TETRAHEDRON_TETRAHEDRON_GROWTH_H

#include <contacts/geometry_contacts_tetrahedron_tetrahedron.h>
#include <overlap/geometry_overlap_tetrahedron_tetrahedron.h>
#include <closest_points/geometry_closest_points_tetrahedron_tetrahedron.h>

#include <tiny.h>

#include <types/geometry_gauss_map_of_convex_polyhedra.h>


#include <vector>
#include <cmath>
#include <cassert>

namespace geometry
{
  namespace details
  {

    template< typename V>
    inline void compute_intersection_points(
                                     TetrahedronEigen<typename V::real_type> const & A
                                     , TetrahedronEigen<typename V::real_type> const & B
      , std::vector<EigenVector3<typename V::real_type>>  & intersections
                                     )
    {
        using T = typename V::real_type;
      intersections.clear();
      intersections.reserve(16);

      std::vector<Triangle<T> >                trianglesA(4u);
      std::vector<Triangle<T> >                trianglesB(4u);
      std::vector<details::UnscaledPlane<T> >  planesA(4u);
      std::vector<details::UnscaledPlane<T> >  planesB(4u);

      for (unsigned int v =0u; v < 4u; ++v)
      {
        trianglesA[v] = get_opposite_face( v, A );
        trianglesB[v] = get_opposite_face( v, B );
        planesA[v]    = details::make_unscaled_plane_eigen(trianglesA[v]);
        planesB[v]    = details::make_unscaled_plane_eigen(trianglesB[v]);
      }

      for (unsigned int v =0u; v < 4u; ++v)
      {
        if( details::inside_planes( B.p(v), planesA) )
          intersections.push_back(B.p(v));

        if( details::inside_planes( A.p(v), planesB) )
          intersections.push_back(A.p(v));
      }

      unsigned int edge_table[6][2] = {
        {0, 1}
        , {0, 2}
        , {0, 3}
        , {1, 2}
        , {1, 3}
        , {2, 3}
      };

      for (unsigned int e = 0u; e < 6u; ++e)
      {
        unsigned int const i = edge_table[e][0u];
        unsigned int const j = edge_table[e][1u];

        const EigenVector3<T>& Ai = A.p(i);
        const EigenVector3<T>& Aj = A.p(j);
        const EigenVector3<T>& Bi = B.p(i);
        const EigenVector3<T>& Bj = B.p(j);

        for (unsigned int p = 0u; p < 4u; ++p)
        {
          details::UnscaledPlane<T> const & planeA = planesA[p];
          details::UnscaledPlane<T> const & planeB = planesB[p];

          if (details::is_crossing_plane(Bi, Bj, planeA) )
          {
            const EigenVector3<T> qB = details::make_intersection(Bi, Bj, planeA);

            if( details::inside_planes(qB, planesA) )
              intersections.push_back(qB);
          }

          if (details::is_crossing_plane(Ai, Aj, planeB) )
          {
            const EigenVector3<T> qA = details::make_intersection( Ai, Aj, planeB);

            if( details::inside_planes(qA, planesB) )
              intersections.push_back(qA);
          }
        }
      }
    }

    template<typename V>
    inline void project_to_plane(V const & n, V const & p, std::vector<V> & intersections)
    {
      typedef typename std::vector<V>::iterator   iterator;

      iterator q   = intersections.begin();
      iterator end = intersections.end();

      for( ; q != end ; ++q )
      {
        (*q) = (*q) - tiny::inner_prod( n, ( (*q) - p ) ) * n;
      }
    }

    template<typename V>
    inline void estimate_overlap(
                                 V const & normal
                                 , std::vector<V> const & intersections
                                 , typename V::real_type & max_val
                                 , typename V::real_type & min_val
                                 )
    {
      using std::min;
      using std::max;

      typedef typename V::real_type                   T;
      typedef typename std::vector<V>::const_iterator iterator;

      iterator p   = intersections.begin();
      iterator end = intersections.end();

      min_val = std::numeric_limits<T>::max();
      max_val = std::numeric_limits<T>::lowest();

      for( ;  p!=end; ++p)
      {
        T const d = tiny::inner_prod( (*p), normal );

        min_val = min( min_val, d );
        max_val = max( max_val, d );
      }

    }

    template<typename T>
    inline void estimate_overlap(
        const EigenVector3<T>& normal
        , std::vector<EigenVector3<T>> const & intersections
        , T& max_val
        , T& min_val
        )
    {
        using std::min;
        using std::max;

        typedef typename std::vector<EigenVector3<T>>::const_iterator iterator;

        iterator p   = intersections.begin();
        iterator end = intersections.end();

        min_val = std::numeric_limits<T>::max();
        max_val = std::numeric_limits<T>::lowest();

        for( ;  p!=end; ++p)
        {
            T const d = dot( (*p), normal );

            min_val = min( min_val, d );
            max_val = max( max_val, d );
        }

    }

    template<typename V>
    inline void filter_unique(
                         std::vector<V> const & intersections
                       , std::vector<V>       & reduced
                       )
    {
      typedef typename V::real_type                     T;
      typedef typename std::vector<V>::const_iterator   const_iterator;
      typedef typename std::vector<V>::iterator         iterator;

      T const accuracy = std::numeric_limits<T>::epsilon()*10;

      reduced.clear();
      reduced.reserve( intersections.size() );

      const_iterator p   = intersections.begin();
      const_iterator end = intersections.end();

      for(; p!= end; ++p)
      {
        bool unique = true;

        iterator q    = reduced.begin();
        iterator endq = reduced.end();
        for(; q != endq; ++q)
        {
          T const similarity = tiny::norm_1( (*q) - (*p) );
          unique = similarity < accuracy ? false : unique;
        }

        if(unique)
          reduced.push_back( (*p) );
      }
    }

    template <typename V>
    inline void make_contacts(V const& normal,
                              typename V::real_type const& depth,
                              std::vector<V> const& positions,
                              ContactsCallback<V>& callback)
    {
        typedef typename std::vector<V>::const_iterator const_iterator;

        const_iterator p = positions.begin();
        const_iterator end = positions.end();

        for (; p != end; ++p) { callback(*p, normal, depth); }
    }

  }// end namespace details



  struct GROWTH {};


  template< typename V>
  inline bool contacts_tetrahedron_tetrahedron(
                                               TetrahedronEigen<typename V::real_type> const & A
                                               , TetrahedronEigen<typename V::real_type> const & B
                                               , ContactsCallback<V> & callback
                                               , std::vector<bool> const & surface_A
                                               , std::vector<bool> const & surface_B
                                               , GROWTH const & /*algorithm_tag*/
  )
  {
    using std::min;
    using std::max;

    typedef typename V::real_type    T;
    typedef typename V::value_traits VT;

    bool overlap = overlap_tetrahedron_tetrahedron(A, B);

    if(!overlap)
      return false;

    T delta = VT::numeric_cast(0.9);

    TetrahedronEigen<T> scaled_A;
    TetrahedronEigen<T> scaled_B;

    //--- Try to quickly find a scale, delta, that for sure results in separation
    while(overlap)
    {
      scaled_A = uniform_scale(delta , A);
      scaled_B = uniform_scale(delta , B);

      overlap = overlap_tetrahedron_tetrahedron(scaled_A, scaled_B);

      delta *= VT::numeric_cast(0.9);
    }

    //--- We now know that there exist a scale value, tau, such that
    //---
    //---    delta < tau < 1 = epsilon
    //---
    //--- Where delta results in separation and epislon gives penetration.
    //---
    //--- We want to find a tau-value such that the separation is in the
    //--- order of max_distance.
    //---
    //---
    T distance;
    EigenVector3<T> sA;
    EigenVector3<T> sB;


    closest_points_tetrahedron_tetrahedron<T>(scaled_A, scaled_B, sA, sB, distance );


    if( distance == std::numeric_limits<T>::max() )
      return false;

    T const max_distance = 2*min( distance, VT::numeric_cast(0.01)* (min(A.get_scale(), B.get_scale() )) );

    T tau = delta;

    const EigenVector3<T> dc   = (B.get_center() - A.get_center());
    T const dcdc = dot(dc, dc);

    unsigned int const max_iterations = 20u;
    unsigned int       iteration = 0u;

    while (distance > max_distance && iteration< max_iterations)  // Separation distance was too big.. we need to enlarge a little
    {
      const EigenVector3<T> ds             = sB - sA;
      T const gap_procentage =  dot(ds, dc) / dcdc;

      T const enlarge = 1 - gap_procentage;

      assert(enlarge < 1.0 );

      if(tau >= enlarge)  // We can not do better
        break;

      tau = enlarge;

      scaled_A = uniform_scale(tau , A);
      scaled_B = uniform_scale(tau , B);

      closest_points_tetrahedron_tetrahedron<T>(scaled_A, scaled_B, sA, sB, distance );

      ++iteration;
    }

    const EigenVector3<T> s = (sB-sA).normalized();

    static GaussMapOfConvexPolyhedra<T> G = make_gauss_map( TetrahedronEigen<T>() );

    //G.set_tolerance_in_degrees(50.0);

    EigenVector3<T> nA;
    EigenVector3<T> nB;

    std::vector<EigenVector3<T>> pointsA;
    std::vector<EigenVector3<T>> pointsB;

    update_gauss_map(G,(A));
    G.search_for_feature((s), pointsA, nA);

    update_gauss_map<T>(G,(B));
    G.search_for_feature(-(s), pointsB, nB);

    // Determine contact normal from highest dimensional feature
    const EigenVector3<T> normal = pointsA.size() >=  pointsB.size() ? (nA) : (-nB);

    std::vector<EigenVector3<T>> intersections;
    details::compute_intersection_points<V>(A, B, intersections);

    if(intersections.empty())
      return false;

    T max_val;
    T min_val;
    details::estimate_overlap(normal, intersections, max_val, min_val);

    T const depth = min_val - max_val;

    const EigenVector3<T> mid =  normal * (max_val + min_val)*0.5f;

    //details::project_to_plane(n, mid, intersections);


    std::vector<V> reduced;
    details::filter_unique(eigenvecofvecToTiny<V>(intersections), reduced);

    details::make_contacts(fromEigen(normal), depth, reduced, callback);

    return reduced.size() > 0u;
  }

}// end namespace geometry

// GEOMETRY_CONTACTS_TETRAHEDRON_TETRAHEDRON_GROWTH_H
#endif
