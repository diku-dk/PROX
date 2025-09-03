#ifndef GEOMETRY_TRANSFORM_H
#define GEOMETRY_TRANSFORM_H

#include <types/geometry_obb.h>
#include <types/geometry_cylinder.h>

#include <tiny.h>

#include <cmath>

namespace geometry
{

  /**
   * TAG Dispatching type for selecting the proper kind of transformation.
   *
   * Use this TAG to specify that it is a "point" that is being transformed from one coordinate system to another.
   */
  struct TRANSFORM_POINT {};

  /**
   * TAG Dispatching type for selecting the proper kind of transformation.
   *
   * Use this TAG to specify that it is a "vector" that is being transformed from one coordinate system to another.
   */
  struct TRANSFORM_VECTOR {};

  template<typename MT>
  inline typename MT::vector3_type transform_from_obb( typename MT::vector3_type const & p, OBB<MT> const & box, TRANSFORM_POINT const & /*tag*/)
  {
    return rotate( box.orientation(), p ) + box.center();
  }

  template<typename MT>
  inline typename MT::vector3_type transform_from_obb( typename MT::vector3_type const & p, OBB<MT> const & box)
  {
    return transform_from_obb(p, box, TRANSFORM_POINT() );
  }

  template<typename MT>
  inline typename MT::vector3_type transform_to_obb( typename MT::vector3_type const & v, OBB<MT> const & box, TRANSFORM_VECTOR const & /*tag*/)
  {
    return rotate( conj( box.orientation() ), v );
  }

  template<typename MT>
  inline typename MT::vector3_type transform_to_obb( typename MT::vector3_type const & p, OBB<MT> const & box, TRANSFORM_POINT const & /*tag*/)
  {
    return rotate( conj( box.orientation() ), (p - box.center()) );
  }

  template<typename MT>
  inline typename MT::vector3_type transform_to_obb( typename MT::vector3_type const & p, OBB<MT> const & box)
  {
    return transform_to_obb( p, box, TRANSFORM_POINT() );
  }

  //EIGEN FUNCTIONS
  template<typename T>
  inline EigenVector3<T> transform_from_obb(const EigenVector3<T>& p, OBBEigen<T> const & box, TRANSFORM_POINT const & /*tag*/)
  {
      return rotate( box.orientation(), p ) + box.center();
  }

  template<typename T>
  inline EigenVector3<T> transform_from_obb( const EigenVector3<T>& p, OBBEigen<T> const & box)
  {
      return transform_from_obb(p, box, TRANSFORM_POINT() );
  }

  template<typename T>
  inline EigenVector3<T> transform_to_obb(const EigenVector3<T>& v, OBBEigen<T> const & box, TRANSFORM_VECTOR const & /*tag*/)
  {
      return rotate( ( box.orientation() ).conjugate(), v );
  }

  template<typename T>
  inline EigenVector3<T> transform_to_obb(const EigenVector3<T>& p, OBBEigen<T> const & box, TRANSFORM_POINT const & /*tag*/)
  {
      EigenQuaternion<T> val = ( box.orientation() ).conjugate();
      EigenVector3<T> val2 = (p - box.center());
      return rotate(val, val2 );
  }

  template<typename T>
  inline EigenVector3<T> transform_to_obb(const EigenVector3<T>& p, OBBEigen<T> const & box)
  {
      return transform_to_obb( p, box, TRANSFORM_POINT() );
  }

  template<typename T>
  inline EigenVector3<T> transform_to_cylinder( const EigenVector3<T> & p, CylinderEigen<T> const & cylinder, TRANSFORM_POINT const & /*tag*/)
  {
    using std::acos;
    using std::min;
    using std::max;



      // A cylinders default orientation is with axis aligned with V::k() axis
    const EigenVector3<T> axis = cross(EigenVector3<T>(0,0,1), cylinder.axis());

    const auto angle = acos(min<T>(1, max<T>(-1, dot({0,0,1}, cylinder.axis()))));

    const EigenQuaternion<T> orientation = Rotateu( angle, axis);

    return rotate( conj( orientation ), (p - cylinder.center()) );
  }

  template<typename T>
  inline EigenVector3<T> transform_to_cylinder(const EigenVector3<T>& p, CylinderEigen<T> const & cylinder )
  {
    return transform_to_cylinder( p, cylinder, TRANSFORM_POINT() );
  }

  template<typename T>
  inline EigenVector3<T> transform_from_cylinder( const EigenVector3<T>& p, CylinderEigen<T> const & cylinder, TRANSFORM_POINT const & /*tag*/)
  {
    using std::acos;
    using std::min;
    using std::max;




      // A cylinders default orientation is with axis aligned with V::k() axis
    const EigenVector3<T> axis        =  cross( EigenVector3<T>(0,0,1), cylinder.axis());

    T const angle = acos(min<T>(1, max<T>(-1, dot(EigenVector3<T>(0,0,1), cylinder.axis()))));

    const EigenQuaternion<T> orientation = Rotateu( angle, axis);

    return rotate( orientation , p ) + cylinder.center();
  }

  template<typename T>
  inline EigenVector3<T> transform_from_cylinder( const EigenVector3<T>& p, CylinderEigen<T> const & cylinder )
  {
    return transform_from_cylinder( p, cylinder, TRANSFORM_POINT());
  }

}// namespace geometry

  // GEOMETRY_TRANSFORM_H
#endif
