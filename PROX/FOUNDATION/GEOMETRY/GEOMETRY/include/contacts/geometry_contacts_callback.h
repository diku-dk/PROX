#ifndef GEOMETRY_CONTACTS_CALLBACK_H
#define GEOMETRY_CONTACTS_CALLBACK_H

#include <tiny_math_types.h>

namespace geometry
{

    template<typename V>
    class ContactsCallback
    {
    public:
        using T = typename V::real_type;
        /**
         * Callback interface for reporting newly found contact points.
         *
         * @param point     The contact point in WCS.
         * @param normal    The contact point normal in WCS.
         * @param distance  The penetration distance measure, negative if overlapping and positive if separation.
         */
        void operator()(const auto& point, const auto& normal, const auto& distance)
        {
            tempParenthesisOperatorImpl(toEigen(point), toEigen(normal), distance);
        }

        virtual void tempParenthesisOperatorImpl(const EigenVector3<T>& point, const EigenVector3<T>& normal,
                                                 const T& distance)
            = 0;
    };

}//namespace geometry

//GEOMETRY_CONTACTS_CALLBACK_H
#endif
