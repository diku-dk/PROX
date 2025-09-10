#ifndef NARROW_BOX_BOX_H
#define NARROW_BOX_BOX_H

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
     * @tparam M                        The math types used.
     *
     * @param A        Shapes from object A.
     * @param B        Shapes from object B.
     * @param input    Current proximity information about the two narrowphase objects.
     * @param results  Upon return any new proximity information will have been added to this container.
     * @param tag      Tag dispatching used for transfering parameter settings to the internal settings.
     */

  template <typename T>
  inline void box_box(typename Geometry<T>::box_container const& A,
                      typename Geometry<T>::box_container const& B,
                      const EigenVector3<T>& tA, const EigenQuaternion<T>& qA,
                      const EigenVector3<T>& tB, const EigenQuaternion<T>& qB,
                      T const& envelope,
                      typename geometry::ContactsCallback<T>& callback)
  {
      using std::min;

      assert(envelope > 0 || !"box_box(): collision envelope must be positive");

      if( A.empty() || B.empty())
        return;

      CoordSysEigen<T> bodyAtoWCS = CoordSysEigen<T>((tA), (qA));
      CoordSysEigen<T> bodyBtoWCS = CoordSysEigen<T>((tB), (qB));

      for (auto a = A.begin(); a != A.end(); ++a)
      {
          for (auto b = B.begin(); b != B.end(); ++b)
          {
              CoordSysEigen<T> shapeAtobodyA = CoordSysEigen<T>(
                  (a->transform().T()), (a->transform().Q()));
              CoordSysEigen<T> shapeBtobodyB = CoordSysEigen<T>(
                  (b->transform().T()), (b->transform().Q()));
              CoordSysEigen<T> shapeAtoWCS = prod(shapeAtobodyA, bodyAtoWCS);
              CoordSysEigen<T> shapeBtoWCS = prod(shapeBtobodyB, bodyBtoWCS);

              geometry::OBBEigen<T> const A = geometry::make_obb<T>(
                  shapeAtoWCS.T(), shapeAtoWCS.Q(), (a->half_extent()));
              geometry::OBBEigen<T> const B = geometry::make_obb<T>(
                  shapeBtoWCS.T(), shapeBtoWCS.Q(), (b->half_extent()));

              geometry::contacts_obb_obb<T>(
                  A, B, envelope * std::min(a->scale(), b->scale()), callback);
        }
      }
    }

  } //namespace detail

} //namespace narrow

// NARROW_BOX_BOX_H
#endif
