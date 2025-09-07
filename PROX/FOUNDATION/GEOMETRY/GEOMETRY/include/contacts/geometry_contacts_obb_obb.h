#ifndef GEOMETRY_CONTACTS_OBB_OBB_H
#define GEOMETRY_CONTACTS_OBB_OBB_H

#include <geometry_transform.h>
#include <geometry_inside.h>
#include <types/geometry_obb.h>
#include <contacts/geometry_contacts_callback.h>
#include <overlap/geometry_overlap_obb_obb.h>

namespace geometry
{

  namespace detail
  {

    class Edge2VertexLookupTable
    {
    protected:

      size_t m_value[12][2];

    public:

      Edge2VertexLookupTable()
      {
        m_value[0][0]  = 0u; m_value[0][1]  = 1u;
        m_value[1][0]  = 1u; m_value[1][1]  = 3u;
        m_value[2][0]  = 3u; m_value[2][1]  = 2u;
        m_value[3][0]  = 2u; m_value[3][1]  = 0u;

        m_value[4][0]  = 4u; m_value[4][1]  = 5u;
        m_value[5][0]  = 5u; m_value[5][1]  = 7u;
        m_value[6][0]  = 7u; m_value[6][1]  = 6u;
        m_value[7][0]  = 6u; m_value[7][1]  = 4u;

        m_value[8][0]  = 0u; m_value[8][1]  = 4u;
        m_value[9][0]  = 1u; m_value[9][1]  = 5u;
        m_value[10][0] = 3u; m_value[10][1] = 7u;
        m_value[11][0] = 2u; m_value[11][1] = 6u;
      }

      /**
       * An OBB has 12 edges the edge index indicates which one. Each edge
       * has two end points (0 or 1) endpoint_index indicates which one
       * we wish to retrive the vertex index for.
       */
      size_t const & operator()(size_t const & edge_idx, size_t const & endpoint_index) const
      {
        return m_value[edge_idx][endpoint_index];
      }

    };

    /**
     * Edge OBB Plane Intersection Test.
     * All points are computed/assumed to be in the local coordinate
     * frame of the OBB.
     *
     * @param a           Edge point 0 of the edge in the box coordinate frame.
     * @param b           Edge point 1 of the edge in the box coordinate frame.
     * @param face_idx    The face to be tested against (a number from 0-5)
     * @param box         The box to be tested against.
     * @param p           Upon return, this hold the intersection point (if any)
     *
     * @return            If intersection exists then the return value is true.
     */
    template<typename T>
    inline bool compute_edge_obb_face_intersection(
                                                   const EigenVector3<T>& a
                                                   , const EigenVector3<T>& b
                                                   , size_t const & face_idx
                                                   , OBBEigen<T> const & box
                                                   , EigenVector3<T>& p
                                                   )
    {
      using std::fabs;

      assert( face_idx<6u || !"compute_edge_obb_face_intersection(): logic error");


      // Decode face_idx argument to determine which face plane direction we clip against (x,y, or z)
      //
      //      -------+--------+----+----+----
      //          face_idx    |  i |  j |  k
      //       base10  base2  |    |    |
      //      -------+--------+----+----+----
      //         0   |  000   |  0 |  1 |  2
      //         1   |  001   |  0 |  1 |  2
      //         2   |  010   |  1 |  2 |  0
      //         3   |  011   |  1 |  2 |  0
      //         4   |  100   |  2 |  0 |  1
      //         5   |  101   |  2 |  0 |  1
      //      -------+--------+----+----+----
      //
      //
      //
      size_t const i = (face_idx >> 1);
      size_t const j = (i + 1) % 3;
      size_t const k = (i + 2) % 3;

      assert( i<3u     || !"compute_edge_obb_face_intersection(): logic error");
      assert( j<3u     || !"compute_edge_obb_face_intersection(): logic error");
      assert( k<3u     || !"compute_edge_obb_face_intersection(): logic error");
      assert( i!=j     || !"compute_edge_obb_face_intersection(): logic error");
      assert( i!=k     || !"compute_edge_obb_face_intersection(): logic error");
      assert( j!=k     || !"compute_edge_obb_face_intersection(): logic error");

      const EigenVector3<T>& e  = box.half_extent();  // just for readabiity
      T const   di = b(i) - a(i);
      T const   E  = (face_idx & 0x0001) ? -e(i) : e(i); // pick the positive or negative face plane orthgonal with i'th direction.

      if(
         ( (a(i) < E) && (b(i) > E) )
         ||
         ( (b(i) < E) && (a(i) > E) )
         )
      {
        //--- Intersection with plane found, now compute intersection point p
        p(i) = E;

        if( fabs(di) > 0)
        {
          T const dj = b(j) - a(j);
          T const dk = b(k) - a(k);

          p(j) = a(j) + (dj/di)*( E - a(i) );
          p(k) = a(k) + (dk/di)*( E - a(i) );
        }
        else
        {
          p(j) = a(j);
          p(k) = a(k);
        }

        //--- If p is inside rectangle of face plane then return true
        if( (-e(j) <= p(j)) && (p(j) <= e(j)) && (-e(k) <= p(k)) && (p(k) <= e(k))  )
        {
          return true;
        }

      }

      return false;
    }

  }// end namespace detail

  /**
   * OBB v OBB Intersection Test and Contact Point Generation.
   *
   * @param A         OBB A.
   * @param B         OOB B.
   * @param envelope  The collision enveleope.
   * @param callback  A pointer to a callback interface which is used to add contact point information
   *
   * @return          If intersection exists then the return value is true.
   */
  template <typename T>
  inline bool contacts_obb_obb(
      OBBEigen<T> const& obbEigenA, OBBEigen<T> const& obbEigenB,
      const T& envelope,
      ContactsCallback<typename tiny::MathTypes<T>::vector3_type>& callback)
  {
    // 2012-06-24 Kenny code review: What about collision envelope, it does nat appear to be used anywhere?
    using std::min;
    using std::max;

    std::vector<EigenVector3<T>> a(8u, EigenVector3<T>(0,0,0));
    a[0] = transform_from_obb( get_local_corner(0, obbEigenA), obbEigenA );
    a[1] = transform_from_obb( get_local_corner(1, obbEigenA), obbEigenA);
    a[2] = transform_from_obb( get_local_corner(2, obbEigenA), obbEigenA);
    a[3] = transform_from_obb( get_local_corner(3, obbEigenA), obbEigenA);
    a[4] = transform_from_obb( get_local_corner(4, obbEigenA), obbEigenA);
    a[5] = transform_from_obb( get_local_corner(5, obbEigenA), obbEigenA);
    a[6] = transform_from_obb( get_local_corner(6, obbEigenA), obbEigenA);
    a[7] = transform_from_obb( get_local_corner(7, obbEigenA), obbEigenA);

    std::vector<EigenVector3<T>> b(8, EigenVector3<T>(0,0,0));
    b[0] = transform_from_obb( get_local_corner(0, obbEigenB), obbEigenB);
    b[1] = transform_from_obb( get_local_corner(1, obbEigenB), obbEigenB);
    b[2] = transform_from_obb( get_local_corner(2, obbEigenB), obbEigenB);
    b[3] = transform_from_obb( get_local_corner(3, obbEigenB), obbEigenB);
    b[4] = transform_from_obb( get_local_corner(4, obbEigenB), obbEigenB);
    b[5] = transform_from_obb( get_local_corner(5, obbEigenB), obbEigenB);
    b[6] = transform_from_obb( get_local_corner(6, obbEigenB), obbEigenB);
    b[7] = transform_from_obb( get_local_corner(7, obbEigenB), obbEigenB);

    EigenVector3<T> n;
    bool const overlap = overlap_obb_obb(a,obbEigenA,b,obbEigenB,n);

    if(! overlap )
      return false;

    assert( is_number( n(0) )    || !"contacts_obb_obb(): nan");
    assert( is_finite( n(0) )    || !"contacts_obb_obb(): inf");
    assert( is_number( n(1) )    || !"contacts_obb_obb(): nan");
    assert( is_finite( n(1) )    || !"contacts_obb_obb(): inf");
    assert( is_number( n(2) )    || !"contacts_obb_obb(): nan");
    assert( is_finite( n(2) )    || !"contacts_obb_obb(): inf");
    assert( fabs(1 - dot(n,n)) < std::numeric_limits<T>::epsilon()*10 || !"contacts_obb_obb(): logic error");

    std::vector<bool> a_in_b( 8u, false);
    a_in_b[0] = inside_obb( a[0], obbEigenB);
    a_in_b[1] = inside_obb( a[1], obbEigenB);
    a_in_b[2] = inside_obb( a[2], obbEigenB);
    a_in_b[3] = inside_obb( a[3], obbEigenB);
    a_in_b[4] = inside_obb( a[4], obbEigenB);
    a_in_b[5] = inside_obb( a[5], obbEigenB);
    a_in_b[6] = inside_obb( a[6], obbEigenB);
    a_in_b[7] = inside_obb( a[7], obbEigenB);

    std::vector<bool> b_in_a( 8u, false);
    b_in_a[0] = inside_obb( b[0], obbEigenA);
    b_in_a[1] = inside_obb( b[1], obbEigenA);
    b_in_a[2] = inside_obb( b[2], obbEigenA);
    b_in_a[3] = inside_obb( b[3], obbEigenA);
    b_in_a[4] = inside_obb( b[4], obbEigenA);
    b_in_a[5] = inside_obb( b[5], obbEigenA);
    b_in_a[6] = inside_obb( b[6], obbEigenA);
    b_in_a[7] = inside_obb( b[7], obbEigenA);

    std::vector<EigenVector3<T>> b_local( 8u, EigenVector3<T>(0,0,0));
    b_local[0] = transform_to_obb( b[0], obbEigenA);
    b_local[1] = transform_to_obb( b[1], obbEigenA);
    b_local[2] = transform_to_obb( b[2], obbEigenA);
    b_local[3] = transform_to_obb( b[3], obbEigenA);
    b_local[4] = transform_to_obb( b[4], obbEigenA);
    b_local[5] = transform_to_obb( b[5], obbEigenA);
    b_local[6] = transform_to_obb( b[6], obbEigenA);
    b_local[7] = transform_to_obb( b[7], obbEigenA);

    std::vector<EigenVector3<T>> a_local( 8u, EigenVector3<T>(0,0,0) );
    a_local[0] = transform_to_obb( a[0], obbEigenB);
    a_local[1] = transform_to_obb( a[1], obbEigenB);
    a_local[2] = transform_to_obb( a[2], obbEigenB);
    a_local[3] = transform_to_obb( a[3], obbEigenB);
    a_local[4] = transform_to_obb( a[4], obbEigenB);
    a_local[5] = transform_to_obb( a[5], obbEigenB);
    a_local[6] = transform_to_obb( a[6], obbEigenB);
    a_local[7] = transform_to_obb( a[7], obbEigenB);

    std::vector<EigenVector3<T>> contacts;

    // If a vertex from OBB is inside the other OBB then generate a contact point
    unsigned int cnt_a_in_b = 0u;
    unsigned int cnt_b_in_a = 0u;

    for(size_t i = 0u; i < 8u; ++i)
    {
      if (a_in_b[i])
      {
        contacts.push_back( a[i] );
        ++cnt_a_in_b;
      }

      if (b_in_a[i])
      {
        contacts.push_back( b[i] );
        ++cnt_b_in_a;
      }
    }


    static detail::Edge2VertexLookupTable const LUT;

    // If an edge from one OBB intersects a face of the other OBB then generate a contact point.
    for(size_t edge_idx=0u; edge_idx<12u; ++edge_idx)
    {
      for(size_t face_idx=0u; face_idx<6u; ++face_idx)
      {
        EigenVector3<T> p;

        if( detail::compute_edge_obb_face_intersection( b_local[ LUT(edge_idx,0) ], b_local[ LUT(edge_idx,1) ], face_idx, obbEigenA, p) )
            contacts.push_back( (transform_from_obb<T>( p, obbEigenA)) );

        if( detail::compute_edge_obb_face_intersection( a_local[ LUT(edge_idx,0) ], a_local[ LUT(edge_idx,1) ], face_idx, obbEigenB, p) )
            contacts.push_back( (transform_from_obb<T>( p, obbEigenB)) );
      }
    }

    if(contacts.empty())
    {
      return false;
    }

    // Now we know we have found a bunch of contacts, we will now try to find
    // the best "contact plane" and project all contact points onto that
    // plane and compute penetration measures with respect to that plane too.
    T min_val = std::numeric_limits<T>::max();
    T max_val = std::numeric_limits<T>::lowest();

    {
      // 2015-11-12 Kenny code review: We assume we have points from both A
      //            and B making up the contact info.. if not depth will be
      //            off... instead it will measure the thinkness of the
      //            embedded box...

      for( auto p = contacts.begin(); p!= contacts.end(); ++p)
      {
        T const d = dot( (*p), n );

        min_val = min( min_val, d );
        max_val = max( max_val, d );
      }

    }

    T const depth = min_val - max_val;

    const EigenVector3<T> mid =  n * (max_val + min_val)*0.5f;

    for( auto p = contacts.begin(); p!= contacts.end(); ++p)
    {
      assert( is_number( (*p)(0) ) || !"contacts_obb_obb(): nan");
      assert( is_finite( (*p)(0) ) || !"contacts_obb_obb(): inf");
      assert( is_number( (*p)(1) ) || !"contacts_obb_obb(): nan");
      assert( is_finite( (*p)(1) ) || !"contacts_obb_obb(): inf");
      assert( is_number( (*p)(2) ) || !"contacts_obb_obb(): nan");
      assert( is_finite( (*p)(2) ) || !"contacts_obb_obb(): inf");

      (*p) = (*p) - ( n).dot( ( (*p) - mid ) ) * n;

      assert( is_number( (*p)(0) ) || !"contacts_obb_obb(): nan");
      assert( is_finite( (*p)(0) ) || !"contacts_obb_obb(): inf");
      assert( is_number( (*p)(1) ) || !"contacts_obb_obb(): nan");
      assert( is_finite( (*p)(1) ) || !"contacts_obb_obb(): inf");
      assert( is_number( (*p)(2) ) || !"contacts_obb_obb(): nan");
      assert( is_finite( (*p)(2) ) || !"contacts_obb_obb(): inf");
    }

    // Some contacts might have been projected to the same point in the
    // contact plane, so we filter away redundant information before
    // using the callback to report the computed contact point.
    typename std::vector<EigenVector3<T>>::iterator p    = contacts.begin();
    for(; p!= contacts.end(); ++p)
    {
      bool unique = true;

      typename std::vector<EigenVector3<T>>::iterator q = contacts.begin();
      for(; q != p; ++q)
      {
          if( ( (*q) - (*p) ).cwiseAbs().maxCoeff()  < std::numeric_limits<T>::epsilon()*10)
        {
          unique = false;
          break;
        }
      }

      if(unique)
          callback( fromEigen((*p)), fromEigen(n), depth);
    }

    return true;
  }

} //namespace geometry

// GEOMETRY_CONTACTS_OBB_OBB_H
#endif
