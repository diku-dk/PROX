#ifndef GEOMETRY_OBB_H
#define GEOMETRY_OBB_H

#include <types/geometry_box.h>

#include <cassert>

namespace geometry
{

/**
   * Oriented Bounding Box (OBBEigen).
   * This is a Box that are subjected to a rigid body transformation (rotation and translation).
   */
template<typename T>
class OBBEigen
{
protected:

    EigenVector3<T> m_center;
    EigenQuaternion<T> m_orientation;
    EigenVector3<T> m_half_extent;

public:

    const EigenVector3<T>& center()      const   { return this->m_center; }
    const EigenQuaternion<T>& orientation() const   { return this->m_orientation; }
    const EigenVector3<T>& half_extent() const   { return this->m_half_extent; }
    EigenVector3<T>& center()              { return this->m_center; }
    EigenQuaternion<T>& orientation()         { return this->m_orientation; }
    EigenVector3<T>& half_extent()         { return this->m_half_extent; }

public:

    OBBEigen()
        : m_center( EigenVector3<T>(0,0,0) )
        , m_orientation( EigenQuaternion<T>::Identity() )
        , m_half_extent( EigenVector3<T>(1,1,1) )
    {}

    OBBEigen(const EigenVector3<T>& center, const EigenQuaternion<T>& orientation, const EigenVector3<T>&  half_extent )
        : m_center( center )
        , m_orientation( orientation )
        , m_half_extent( half_extent )
    {}

    OBBEigen(OBBEigen const & box)
    {
        *this = box;
    }

    OBBEigen & operator=(OBBEigen const & box)
    {
        if( this != &box)
        {
            this->m_center      = box.m_center;
            this->m_orientation = box.m_orientation;
            this->m_half_extent = box.m_half_extent;
        }
        return *this;
    }

};

template<typename T>
inline OBBEigen<T> make_obb(
    const EigenVector3<T>& center
    , const EigenQuaternion<T>& orientation
    , const EigenVector3<T>&  half_extent
    )
{
    return OBBEigen<T>( center, orientation, half_extent);
}


template<typename T>
inline OBBEigen<T> make_obb(
    const EigenVector3<T>& center
    , const EigenQuaternion<T>& orientation
    , BoxEigen<T> const & box
    )
{
    return OBBEigen<T>( center, orientation, box.half_extent() );
}

/**
   * Returns the local-space coordiantes of a specifed corner point.
   *
   * To retrive the corner point in global (world) space use the transfrom_from_obb function.
   */
template<typename T>
inline EigenVector3<T> get_local_corner(size_t const & idx, OBBEigen<T> const & box)
{
    assert(idx < 8u || !"get_local_corner(): Idx was out of range");

    T const x = ((idx   ) & 0x0001u) ? box.half_extent()(0) : -box.half_extent()(0);
    T const y = ((idx>>1) & 0x0001u) ? box.half_extent()(1) : -box.half_extent()(1);
    T const z = ((idx>>2) & 0x0001u) ? box.half_extent()(2) : -box.half_extent()(2);

    return EigenVector3<T>(x,y,z);
}

}

// GEOMETRY_OBB_H
#endif
