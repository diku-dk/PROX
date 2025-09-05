#ifndef NARROW_BOX_BOX_H
#define NARROW_BOX_BOX_H

#include "narrow_geometry.h"
#include <geometry.h>

#include <tiny.h>

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

    template<typename M>
    inline void box_box(
                        typename Geometry<M>::box_container const & A
                        , typename Geometry<M>::box_container const & B
                        , typename M::vector3_type const & tA
                        , typename M::quaternion_type const & qA
                        , typename M::vector3_type const & tB
                        , typename M::quaternion_type const & qB
                        , typename M::real_type const & envelope
                        , typename geometry::ContactsCallback<typename M::vector3_type> & callback
                        )
    {
      using std::min;

      typedef typename M::value_traits    VT;
      typedef typename M::coordsys_type   C;
      using T = M::real_type;

      typedef typename Geometry<M>::box_container::const_iterator    box_iterator;

      assert( envelope > 0 || !"box_box(): collision envelope must be positive");

      if( A.empty() || B.empty())
        return;

      CoordSysEigen<T> bodyAtoWCS = CoordSysEigen<T>(toEigen(tA), toEigen(qA));
      CoordSysEigen<T> bodyBtoWCS = CoordSysEigen<T>(toEigen(tB), toEigen(qB));

      for( box_iterator a = A.begin(); a!=A.end(); ++a )
      {
        for( box_iterator b = B.begin(); b!=B.end(); ++b )
        {
            CoordSysEigen<T> shapeAtobodyA = CoordSysEigen<T>(toEigen(a->transform().T()), toEigen(a->transform().Q()));
            CoordSysEigen<T> shapeBtobodyB = CoordSysEigen<T>(toEigen(b->transform().T()), toEigen(b->transform().Q()));
          CoordSysEigen<T> shapeAtoWCS = prod(shapeAtobodyA, bodyAtoWCS);
          CoordSysEigen<T> shapeBtoWCS = prod(shapeBtobodyB, bodyBtoWCS);

          geometry::OBBEigen<T> const A = geometry::make_obb<T>( shapeAtoWCS.T(), shapeAtoWCS.Q(), toEigen(a->half_extent()));
          geometry::OBBEigen<T> const B = geometry::make_obb<T>( shapeBtoWCS.T(), shapeBtoWCS.Q(), toEigen(b->half_extent()));

          geometry::contacts_obb_obb<M>(
                                     A
                                     , B
                                     , envelope * min(a->scale(), b->scale())
                                     , callback
                                     );
        }
      }
    }

  } //namespace detail

} //namespace narrow

// NARROW_BOX_BOX_H
#endif
