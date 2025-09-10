#ifndef GEOMETRY_OVERLAP_OBB_OBB
#define GEOMETRY_OVERLAP_OBB_OBB

#include <types/geometry_obb.h>

#include <eigenhelperall.h>

#include <cmath>
#include <vector>


///////////EIGEN
  namespace geometry
  {
  namespace detail
  {
  /**
     * This function generates all possible separating axes between the two specified OBBs.
     */
  template<typename T>
  inline void compute_obb_obb_sat_axes(
      OBBEigen<T> const & A
      , OBBEigen<T> const & B
      , std::vector<EigenVector3<T>> & axes
      )
  {
      using std::min;
      using std::max;
      using std::fabs;


      const EigenMatrix3<T> R_a(A.orientation());
      const EigenMatrix3<T> R_b(B.orientation());

      const EigenVector3<T> A0  = R_a.col(0);
      const EigenVector3<T> A1  = R_a.col(1);
      const EigenVector3<T> A2  = R_a.col(2);

      assert( fabs( 1 - dot(A0,A0) ) < std::numeric_limits<T>::epsilon()*10 || !"compute_obb_obb_sat_axes(): logic error");
      assert( fabs( 1 - dot(A1,A1) ) < std::numeric_limits<T>::epsilon()*10 || !"compute_obb_obb_sat_axes(): logic error");
      assert( fabs( 1 - dot(A2,A2) ) < std::numeric_limits<T>::epsilon()*10 || !"compute_obb_obb_sat_axes(): logic error");
      assert( fabs( dot(A0,A1) ) < std::numeric_limits<T>::epsilon()*10             || !"compute_obb_obb_sat_axes(): logic error");
      assert( fabs( dot(A0,A2) ) < std::numeric_limits<T>::epsilon()*10             || !"compute_obb_obb_sat_axes(): logic error");
      assert( fabs( dot(A1,A2) ) < std::numeric_limits<T>::epsilon()*10             || !"compute_obb_obb_sat_axes(): logic error");

      const EigenVector3<T> B0  = R_b.col(0);
      const EigenVector3<T> B1  = R_b.col(1);
      const EigenVector3<T> B2  = R_b.col(2);


      assert( fabs( 1 - dot(B0,B0) ) < std::numeric_limits<T>::epsilon()*10 || !"compute_obb_obb_sat_axes(): logic error");
      assert( fabs( 1 - dot(B1,B1) ) < std::numeric_limits<T>::epsilon()*10 || !"compute_obb_obb_sat_axes(): logic error");
      assert( fabs( 1 - dot(B2,B2) ) < std::numeric_limits<T>::epsilon()*10 || !"compute_obb_obb_sat_axes(): logic error");
      assert( fabs(dot(B0,B1) ) < std::numeric_limits<T>::epsilon()*10             || !"compute_obb_obb_sat_axes(): logic error");
      assert( fabs(dot(B0,B2) ) < std::numeric_limits<T>::epsilon()*10             || !"compute_obb_obb_sat_axes(): logic error");
      assert( fabs(dot(B1,B2) ) < std::numeric_limits<T>::epsilon()*10             || !"compute_obb_obb_sat_axes(): logic error");

      axes.resize(15u);

      axes[0] = A0;
      axes[1] = A1;
      axes[2] = A2;

      axes[3] = B0;
      axes[4] = B1;
      axes[5] = B2;

      axes[6]     = cross(A0, B0);
      axes[7]     = cross(A0, B1);
      axes[8]     = cross(A0, B2);

      axes[9]     = cross(A1, B0);
      axes[10]    = cross(A1, B1);
      axes[11]    = cross(A1, B2);

      axes[12]    = cross(A2, B0);
      axes[13]    = cross(A2, B1);
      axes[14]    = cross(A2, B2);


      T const l6  = norm( axes[6]  );
      T const l7  = norm( axes[7]  );
      T const l8  = norm( axes[8]  );
      T const l9  = norm( axes[9]  );
      T const l10 = norm( axes[10] );
      T const l11 = norm( axes[11] );
      T const l12 = norm( axes[12] );
      T const l13 = norm( axes[13] );
      T const l14 = norm( axes[14] );

      axes[6]     = (l6  > std::numeric_limits<T>::epsilon()*10 ) ? axes[6]  / l6  : A0;
      axes[7]     = (l7  > std::numeric_limits<T>::epsilon()*10 ) ? axes[7]  / l7  : A0;
      axes[8]     = (l8  > std::numeric_limits<T>::epsilon()*10 ) ? axes[8]  / l8  : A0;

      axes[9]     = (l9  > std::numeric_limits<T>::epsilon()*10 ) ? axes[9]  / l9  : A1;
      axes[10]    = (l10 > std::numeric_limits<T>::epsilon()*10 ) ? axes[10] / l10 : A1;
      axes[11]    = (l11 > std::numeric_limits<T>::epsilon()*10 ) ? axes[11] / l11 : A1;

      axes[12]    = (l12 > std::numeric_limits<T>::epsilon()*10 ) ? axes[12] / l12 : A2;
      axes[13]    = (l13 > std::numeric_limits<T>::epsilon()*10 ) ? axes[13] / l13 : A2;
      axes[14]    = (l14 > std::numeric_limits<T>::epsilon()*10 ) ? axes[14] / l14 : A2;

      assert( is_number( axes[0](0) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[0](0) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( is_number( axes[0](1) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[0](1) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( is_number( axes[0](2) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[0](2) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( fabs( 1 - dot(axes[0],axes[0]) ) < std::numeric_limits<T>::epsilon()*10 || !"compute_obb_obb_sat_axes(): logic error");

      assert( is_number( axes[1](0) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[1](0) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( is_number( axes[1](1) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[1](1) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( is_number( axes[1](2) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[1](2) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( fabs( 1 - dot(axes[1],axes[1]) ) < std::numeric_limits<T>::epsilon()*10 || !"compute_obb_obb_sat_axes(): logic error");

      assert( is_number( axes[2](0) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[2](0) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( is_number( axes[2](1) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[2](1) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( is_number( axes[2](2) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[2](2) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( fabs( 1 - dot(axes[2],axes[2]) ) < std::numeric_limits<T>::epsilon()*10 || !"compute_obb_obb_sat_axes(): logic error");

      assert( is_number( axes[3](0) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[3](0) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( is_number( axes[3](1) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[3](1) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( is_number( axes[3](2) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[3](2) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( fabs( 1 - dot(axes[3],axes[3]) ) < std::numeric_limits<T>::epsilon()*10 || !"ompute_obb_obb_sat_axes(): logic error");

      assert( is_number( axes[4](0) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[4](0) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( is_number( axes[4](1) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[4](1) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( is_number( axes[4](2) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[4](2) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( fabs( 1 - dot(axes[4],axes[4]) ) < std::numeric_limits<T>::epsilon()*10 || !"ompute_obb_obb_sat_axes(): logic error");

      assert( is_number( axes[5](0) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[5](0) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( is_number( axes[5](1) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[5](1) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( is_number( axes[5](2) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[5](2) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( fabs( 1 - dot(axes[5],axes[5]) ) < std::numeric_limits<T>::epsilon()*10 || !"ompute_obb_obb_sat_axes(): logic error");

      assert( is_number( axes[6](0) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[6](0) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( is_number( axes[6](1) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[6](1) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( is_number( axes[6](2) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[6](2) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( fabs( 1 - dot(axes[6],axes[6]) ) < std::numeric_limits<T>::epsilon()*10 || !"ompute_obb_obb_sat_axes(): logic error");

      assert( is_number( axes[7](0) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[7](0) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( is_number( axes[7](1) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[7](1) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( is_number( axes[7](2) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[7](2) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( fabs( 1 - dot(axes[7],axes[7]) ) < std::numeric_limits<T>::epsilon()*10 || !"ompute_obb_obb_sat_axes(): logic error");

      assert( is_number( axes[8](0) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[8](0) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( is_number( axes[8](1) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[8](1) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( is_number( axes[8](2) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[8](2) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( fabs( 1 - dot(axes[8],axes[8]) ) < std::numeric_limits<T>::epsilon()*10 || !"ompute_obb_obb_sat_axes(): logic error");

      assert( is_number( axes[9](0) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[9](0) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( is_number( axes[9](1) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[9](1) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( is_number( axes[9](2) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[9](2) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( fabs( 1 - dot(axes[9],axes[9]) ) < std::numeric_limits<T>::epsilon()*10 || !"ompute_obb_obb_sat_axes(): logic error");

      assert( is_number( axes[10](0) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[10](0) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( is_number( axes[10](1) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[10](1) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( is_number( axes[10](2) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[10](2) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( fabs( 1 - dot(axes[10],axes[10]) ) < std::numeric_limits<T>::epsilon()*10 || !"compute_obb_obb_sat_axes(): logic error");

      assert( is_number( axes[11](0) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[11](0) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( is_number( axes[11](1) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[11](1) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( is_number( axes[11](2) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[11](2) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( fabs( 1 - dot(axes[11],axes[11]) ) < std::numeric_limits<T>::epsilon()*10 || !"compute_obb_obb_sat_axes(): logic error");

      assert( is_number( axes[12](0) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[12](0) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( is_number( axes[12](1) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[12](1) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( is_number( axes[12](2) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[12](2) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( fabs( 1 - dot(axes[12],axes[12]) ) < std::numeric_limits<T>::epsilon()*10 || !"compute_obb_obb_sat_axes(): logic error");

      assert( is_number( axes[13](0) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[13](0) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( is_number( axes[13](1) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[13](1) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( is_number( axes[13](2) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[13](2) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( fabs( 1 - dot(axes[13],axes[13]) ) < std::numeric_limits<T>::epsilon()*10 || !"compute_obb_obb_sat_axes(): logic error");

      assert( is_number( axes[14](0) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[14](0) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( is_number( axes[14](1) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[14](1) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( is_number( axes[14](2) ) || !"compute_obb_obb_sat_axes(): nan");
      assert( is_finite( axes[14](2) ) || !"compute_obb_obb_sat_axes(): inf");
      assert( fabs( 1 - dot(axes[14],axes[14]) ) < std::numeric_limits<T>::epsilon()*10 || !"compute_obb_obb_sat_axes(): logic error");
  }

  }// end namespace detail

  /**
   * Extended OBB versus OBB Separating Axis(SAT) overlap Test. This version
   * provides geometric feedback information that can be helpful to determine
   * the "best" separation direction (the normal direction of the overlap).
   *
   * @param a   Corner points of OBB A.
   * @param A   OBB A.
   * @param b   Corner points of OBB B.
   * @param B   OOB B.
   * @param n   Upon return, this holds a potential normal for the overlap.
   *
   * @return    If overlap exists then the return value is true.
   */
  template<typename T>
  inline bool overlap_obb_obb(
      std::vector<EigenVector3<T>> const & a
      , OBBEigen<T> const & A
      , std::vector<EigenVector3<T>> const & b
      , OBBEigen<T> const & B
      , EigenVector3<T>& n
      )
  {
      using std::min;
      using std::max;
      using std::fabs;

      std::vector<EigenVector3<T>> axes;
      detail::compute_obb_obb_sat_axes( A,B, axes);

      std::vector<T> a_min( 15u, std::numeric_limits<T>::max() );
      std::vector<T> b_min( 15u, std::numeric_limits<T>::max() );
      std::vector<T> a_max( 15u, std::numeric_limits<T>::lowest()  );
      std::vector<T> b_max( 15u, std::numeric_limits<T>::lowest()  );

      T min_overlap = std::numeric_limits<T>::lowest();

      for(size_t i=0u;i < 15u; ++i)
      {

          for( typename std::vector<EigenVector3<T>>::const_iterator p_a = a.begin(); p_a != a.end(); ++p_a)
          {
              T const d = dot( (*p_a), axes[i] );
              a_min[i] = min( a_min[i], d);
              a_max[i] = max( a_max[i], d);
          }
          for( typename std::vector<EigenVector3<T>>::const_iterator p_b = b.begin(); p_b != b.end(); ++p_b)
          {
              T const d = dot( (*p_b), axes[i] );
              b_min[i] = min( b_min[i], d);
              b_max[i] = max( b_max[i], d);
          }

          assert( is_number( a_min[i] ) || !"overlap_obb_obb(); nan");
          assert( is_finite( a_min[i] ) || !"overlap_obb_obb(): inf");
          assert( is_number( a_max[i] ) || !"overlap_obb_obb(); nan");
          assert( is_finite( a_max[i] ) || !"overlap_obb_obb(): inf");

          assert( is_number( b_min[i] ) || !"overlap_obb_obb(); nan");
          assert( is_finite( b_min[i] ) || !"overlap_obb_obb(): inf");
          assert( is_number( b_max[i] ) || !"overlap_obb_obb(); nan");
          assert( is_finite( b_max[i] ) || !"overlap_obb_obb(): inf");

          if(a_max[i] < b_min[i])
              return false;

          if(b_max[i] < a_min[i])
              return false;

          if(a_min[i] <= b_min[i] &&  b_min[i] <= a_max[i])
          {
              T const overlap = b_min[i] - a_max[i];

              if(overlap > min_overlap)
              {
                  min_overlap = overlap;
                  n = axes[i];
              }
          }

          if(b_min[i] <= a_min[i] &&  a_min[i] <= b_max[i])
          {
              T const overlap = a_min[i] - b_max[i];

              if(overlap > min_overlap)
              {
                  min_overlap = overlap;
                  n = -axes[i];
              }
          }

          assert( is_number( min_overlap ) || !"overlap_obb_obb(); nan");
          assert( is_finite( min_overlap ) || !"overlap_obb_obb(): inf");
      }

      return (min_overlap <= 0);
  }

  /**
   * This version of the OBB verusus OBB overlap test is faster and does not
   * care about finding a sensible separation direction.
   */
  template<typename T>
  inline bool overlap_obb_obb(
      OBBEigen<T> const & A
      , OBBEigen<T> const & B
      )
  {
      using std::min;
      using std::max;
      using std::fabs;
      std::vector<EigenVector3<T>> axes;
      detail::compute_obb_obb_sat_axes( A,B, axes);

      std::vector<T> a_min( 15u, std::numeric_limits<T>::max() );
      std::vector<T> b_min( 15u, std::numeric_limits<T>::max() );
      std::vector<T> a_max( 15u, std::numeric_limits<T>::lowest()  );
      std::vector<T> b_max( 15u, std::numeric_limits<T>::lowest()  );

      const EigenMatrix3<T> R_a(A.orientation());
      const EigenMatrix3<T> R_b(B.orientation());

      const EigenVector3<T> ai_extent  = R_a.col(0) * A.half_extent()(0);
      const EigenVector3<T> aj_extent  = R_a.col(1) * A.half_extent()(1);
      const EigenVector3<T> ak_extent  = R_a.col(2) * A.half_extent()(2);
      const EigenVector3<T> bi_extent  = R_b.col(0) * B.half_extent()(0);
      const EigenVector3<T> bj_extent  = R_b.col(1) * B.half_extent()(1);
      const EigenVector3<T> bk_extent  = R_b.col(2) * B.half_extent()(2);

      for(size_t i=0u;i < 15u; ++i)
      {
          T const ae    = fabs(dot( ai_extent, axes[i] ))
                     + fabs(dot( aj_extent, axes[i] ))
                     + fabs(dot( ak_extent, axes[i] ))
              ;

          T const ac    = dot( A.center(), axes[i] );
          T const a_min = ac - ae;
          T const a_max = ac + ae;

          T const be    = fabs(dot( bi_extent, axes[i] ))
                     + fabs(dot( bj_extent, axes[i] ))
                     + fabs(dot( bk_extent, axes[i] ))
              ;

          T const bc    = dot( B.center(), axes[i] );
          T const b_min = bc - be;
          T const b_max = bc + be;

          assert( is_number( a_min ) || !"overlap_obb_obb(); nan");
          assert( is_finite( a_min ) || !"overlap_obb_obb(): inf");
          assert( is_number( a_max ) || !"overlap_obb_obb(); nan");
          assert( is_finite( a_max ) || !"overlap_obb_obb(): inf");

          assert( is_number( b_min ) || !"overlap_obb_obb(); nan");
          assert( is_finite( b_min ) || !"overlap_obb_obb(): inf");
          assert( is_number( b_max ) || !"overlap_obb_obb(); nan");
          assert( is_finite( b_max ) || !"overlap_obb_obb(): inf");

          if(a_max < b_min)
              return false;

          if(b_max < a_min)
              return false;
      }

      return true;
  }

  } //namespace geometry

  // GEOMETRY_OVERLAP_OBB_OBB
#endif
