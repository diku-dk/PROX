#ifndef GEOMETRY_SUPPORT_MAPPING
#define GEOMETRY_SUPPORT_MAPPING

#include "tiny_math_types.h"

namespace geometry
{

  template<typename T>
  class SupportMapping
  {
  public:
    /**
     * Get Support Point of Shape.
     *
     * A support point is simply the extreme point of a shape in a given direction.
     *
     * @param v    The support direction.
     *
     * @return     The support point.
     */
      virtual EigenVector3<T> get_support_point(EigenVector3<T> dir) const = 0;
  };

} // namespace geometry

// GEOMETRY_SUPPORT_MAPPING
#endif
