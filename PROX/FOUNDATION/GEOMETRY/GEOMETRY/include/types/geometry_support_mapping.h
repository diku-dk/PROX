#ifndef GEOMETRY_SUPPORT_MAPPING
#define GEOMETRY_SUPPORT_MAPPING

#include "tiny_math_types.h"

namespace geometry
{

  template<typename V>
  class SupportMapping
  {
  public:
      using T = typename V::real_type;
    /**
     * Get Support Point of Shape.
     *
     * A support point is simply the extreme point of a shape in a given direction.
     *
     * @param v    The support direction.
     *
     * @return     The support point.
     */
    virtual V get_support_point( V const & v ) const = 0;

      virtual EigenVector3<T> get_support_point(EigenVector3<T> dir) const = 0;

  };

} // namespace geometry

// GEOMETRY_SUPPORT_MAPPING
#endif
