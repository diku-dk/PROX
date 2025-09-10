#ifndef GEOMETRY_CONTACTS_CALLBACK_H
#define GEOMETRY_CONTACTS_CALLBACK_H

#include <eigenhelperall.h>

namespace geometry
{

template <typename T>
    //requires(std::is_floating_point_v<T>)
class ContactsCallback
{
public:
        /**
         * Callback interface for reporting newly found contact points.
         *
         * @param point     The contact point in WCS.
         * @param normal    The contact point normal in WCS.
         * @param distance  The penetration distance measure, negative if overlapping and positive if separation.
         */
    void operator()(const EigenVector3<T>& point, const EigenVector3<T>& normal,
                    const T& distance)
    {
        tempParenthesisOperatorImpl((point), (normal), distance);
    }

    virtual void tempParenthesisOperatorImpl(const EigenVector3<T>& point,
                                             const EigenVector3<T>& normal,
                                             const T& distance)
        = 0;
};

}//namespace geometry

//GEOMETRY_CONTACTS_CALLBACK_H
#endif
