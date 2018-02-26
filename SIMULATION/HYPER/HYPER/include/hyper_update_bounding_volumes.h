#ifndef HYPER_UPDATE_BOUNDING_VOLUMES_H
#define HYPER_UPDATE_BOUNDING_VOLUMES_H

#include <hyper_body.h>

#include <geometry_make_intersection.h>

#include <cassert>
#include <vector>

namespace hyper
{

  struct fast {};
  struct smallest {};

  template<typename MT>
  inline void update_bounding_volumes( Body<MT> & body, fast const & /*tag*/ )
  {
    using std::fabs;

    typedef typename MT::real_type    T;
    typedef typename MT::vector3_type V;
    typedef typename MT::value_traits VT;

    //--- We assume BVH has been updated, then we get root volume and
    //--- from this we get radius and center information.

    size_t const K = 4u;  // Number of slabs
    size_t const N = 2*K; // Number of planes

    geometry::DOP<T,8>            const & dop = body.m_tree.m_root;
    geometry::DirectionTable<V,4> const   DT  = geometry::DirectionTableHelper<V,4>::make();

    //--- Generate array of all planes of kDOP ---------------------------------
    std::vector< geometry::Plane<V> > planes;
    planes.resize(N);

    for(size_t i = 0u; i < K; ++i)
    {
      geometry::Interval<T> const & slab = dop(i);

      T const l = slab.lower();
      T const u = slab.upper();
      V const n = DT(i);

      geometry::Plane<V> const P_lower = geometry::make_plane( -n, -l );
      geometry::Plane<V> const P_upper = geometry::make_plane(  n,  u );

      planes[ 2u*i + 0u ] = P_lower;
      planes[ 2u*i + 1u ] = P_upper;
    }

    //--- Generate intersection points of all planes ---------------------------
    std::vector< V > points;

    for(size_t i = 0u; i < N; ++i)
    {
      for(size_t j = i+1u; j < N; ++j)
      {
        for(size_t k = j+1u; k < N; ++k)
        {
          geometry::Plane<V> const & Pi = planes[i];
          geometry::Plane<V> const & Pj = planes[j];
          geometry::Plane<V> const & Pk = planes[k];

          bool const ij_is_parallel = fabs( inner_prod(Pi.n(), Pj.n()) ) > VT::numeric_cast(0.9999);
          bool const ik_is_parallel = fabs( inner_prod(Pi.n(), Pk.n()) ) > VT::numeric_cast(0.9999);
          bool const jk_is_parallel = fabs( inner_prod(Pj.n(), Pk.n()) ) > VT::numeric_cast(0.9999);

          if(ij_is_parallel)
            continue;

          if(ik_is_parallel)
            continue;

          if(jk_is_parallel)
            continue;

          V const p = geometry::make_intersection(Pi,Pj,Pk); // Get intersection point of the three planes

          points.push_back( p );

        }
      }
    }

    //--- Project points onto surface of kDOP ----------------------------------
    for(size_t i = 0u; i < N; ++i)
    {
      geometry::Plane<V> const & P = planes[i];

      for(size_t k = 0u; k < points.size(); ++k)
      {
        V const & q = points[k];
        T const   d = geometry::get_signed_distance( q, P);

        //--- If point outside plane then push it back
        points[k] = (d>VT::zero()) ? q - P.n()*d : q;
      }
    }

    //--- Now we can finally find radius, center, min and max corners of kDOP --
    V max_q = points[0];
    V min_q = points[0];
    V sum_q = points[0];
    for(size_t k = 1u; k < points.size(); ++k)
    {
      V const & q = points[k];

      max_q = max(q,max_q);
      min_q = min(q,min_q);
      sum_q = sum_q + q;
    }

    V const center = sum_q / VT::numeric_cast(  points.size()  );

    T radius = VT::zero();
    for(size_t k = 0u; k < points.size(); ++k)
    {
      V const & q = points[k];
      T const d = norm(q-center);
      radius = (d > radius) ? d : radius;
    }

    body.m_radius = radius;
    body.m_center = center;
    body.m_min_point = min_q;
    body.m_max_point = max_q;
  }

  template<typename MT>
  inline void update_bounding_volumes( Body<MT> & body, smallest const & /*tag*/ )
  {
    using std::fabs;

    typedef typename MT::real_type    T;
    typedef typename MT::vector3_type V;
    typedef typename MT::value_traits VT;

    V max_q = V::make( body.m_X(0u), body.m_Y(0u), body.m_Z(0u) );
    V min_q = V::make( body.m_X(0u), body.m_Y(0u), body.m_Z(0u) );
    V sum_q = V::make( body.m_X(0u), body.m_Y(0u), body.m_Z(0u) );
    for(size_t k = 1u; k < body.m_X.size(); ++k)
    {
      V const q = V::make( body.m_X(k), body.m_Y(k), body.m_Z(k) );

      max_q = max(q,max_q);
      min_q = min(q,min_q);
      sum_q = sum_q + q;
    }

    V const center = sum_q / VT::numeric_cast(  body.m_X.size()  );

    T radius = VT::zero();
    for(size_t k = 0u; k < body.m_X.size(); ++k)
    {
      V const q = V::make( body.m_X(k), body.m_Y(k), body.m_Z(k) );
      T const d = norm(q-center);
      radius = (d > radius) ? d : radius;
    }

    body.m_radius = radius;
    body.m_center = center;
    body.m_min_point = min_q;
    body.m_max_point = max_q;

  }

  template<typename MT>
  inline void update_bounding_volumes( Body<MT> & body )
  {
    //update_bounding_volumes( body, fast() );
    update_bounding_volumes( body, smallest() );
  }


} // namespace hyper

// HYPER_UPDATE_BOUNDING_VOLUMES_H
#endif 
