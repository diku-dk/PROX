#ifndef GEOMETRY_CONTACT_POINT_H
#define GEOMETRY_CONTACT_POINT_H

#include <tiny_math_types.h>

namespace geometry
{

template <typename Body, typename T> class ContactPoint
{
public:
    EigenVector3<T> position = {0, 0, 0};
    EigenVector3<T> normal = {0, 0, 0};
    T depth = 0;
    const Body* bodyI = nullptr;
    const Body* bodyJ = nullptr;

    auto operator<=>(const ContactPoint<Body, T>& other) const
    {
        // Match previous weird comparison order.
        return std::tie(bodyI, bodyJ, position.y(), position.x(), position.z(), normal.y(), normal.x(), normal.z(),
                        depth)
           <=> std::tie(other.bodyI, other.bodyJ, other.position.y(), other.position.x(), other.position.z(),
                        other.normal.y(), other.normal.x(), other.normal.z(), other.depth);
    }

    bool operator==(const ContactPoint<Body, T>& other) const = default;
};

}//namespace geometry

//GEOMETRY_CONTACT_POINT_H
#endif
