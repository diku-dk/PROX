#ifndef NARROW_SPHERE_SPHERE_H
#define NARROW_SPHERE_SPHERE_H

#include "narrow_geometry.h"
#include <geometry.h>

#include <eigenhelperall.h>

namespace narrow
{

  namespace detail
  {

    /**
     * Dispatch pair-wise testing of shapes from two different objects.
     *
     * @tparam shape_containerA         A container type of shape types coming from the narrowphase object A.
     * @tparam shape_containerB         A container type of (possible another then A)  of shape types coming from the narrowphase object B.
     * @tparam M                        The math types used.
     * @tparam proximity_container      A container type of proximity info.
     * @tparam Params                   A Policy/Trait class holding parameter settings for the internal workings of the CONVEX algorithm.
     *
     * @param A        Shapes from object A.
     * @param B        Shapes from object B.
     * @param input    Current proximity information about the two narrowphase objects.
     * @param results  Upon return any new proximity information will have been added to this container.
     * @param tag      Tag dispatching used for transfering parameter settings to the internal settings.
     */

  template <typename T>
  inline void sphere_sphere(typename Geometry<T>::sphere_container const& A,
                            typename Geometry<T>::sphere_container const& B,
                            const EigenVector3<T>& tA,
                            const EigenQuaternion<T>& qA,
                            const EigenVector3<T>& tB,
                            const EigenQuaternion<T>& qB, T const& envelope,
                            typename geometry::ContactsCallback<T>& callback)
  {
      using std::min;

      typedef typename Geometry<T>::sphere_container::const_iterator
          sphere_iterator;

      assert( envelope > 0 || !"sphere_sphere(): collision envelope must be positive");

      if( A.empty() || B.empty())
        return;

      CoordSysEigen<T> bodyAtoWCS = CoordSysEigen<T>(tA, qA);
      CoordSysEigen<T> bodyBtoWCS = CoordSysEigen<T>(tB, qB);

      for( sphere_iterator a = A.begin(); a!=A.end(); ++a )
      {
        for( sphere_iterator b = B.begin(); b!=B.end(); ++b )
        {
            CoordSysEigen<T> shapeAtobodyA
                = CoordSysEigen<T>((a->transform().T()), (a->transform().Q()));
            CoordSysEigen<T> shapeBtobodyB
                = CoordSysEigen<T>((b->transform().T()), (b->transform().Q()));
            CoordSysEigen<T> shapeAtoWCS = prod(shapeAtobodyA, bodyAtoWCS);
            CoordSysEigen<T> shapeBtoWCS = prod(shapeBtobodyB, bodyBtoWCS);

            geometry::Sphere<T> const A
                = geometry::make_sphere((shapeAtoWCS.T()), a->radius());
            geometry::Sphere<T> const B
                = geometry::make_sphere((shapeBtoWCS.T()), b->radius());

            geometry::contacts_sphere_sphere<T>(
                A, B, envelope * min(a->scale(), b->scale()), callback);
        }
      }
  }

  } // namespace detail

} //namespace narrow

// NARROW_SPHERE_SPHERE_H
#endif
