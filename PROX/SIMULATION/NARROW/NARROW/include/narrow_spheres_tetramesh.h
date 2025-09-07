#ifndef NARROW_SPHERES_TETRAMESH_H
#define NARROW_SPHERES_TETRAMESH_H

#include <narrow_object.h>
#include <narrow_geometry.h>

#include <kdop_single_traversal.h>

namespace narrow
{

  namespace detail
  {

  template <typename M>
  inline void spheres_tetramesh(
      typename Geometry<typename M::real_type>::sphere_container const& A,
      const EigenVector3<typename M::real_type>& tA,
      const EigenQuaternion<typename M::real_type>& qA,
      Object<typename M::real_type> const& objB,
      Geometry<typename M::real_type> const& geoB,
      typename geometry::ContactsCallback<typename M::vector3_type>& callback,
      bool const& should_flip)
  {
      typedef typename M::vector3_type    V;
      typedef typename M::real_type       T;

      typedef typename Geometry<typename M::real_type>::sphere_container::
          const_iterator sphere_iterator;

      if( A.empty() )
        return;

      CoordSysEigen<T> const bodyAtoWCS = CoordSysEigen<T>(tA, qA);

      for( sphere_iterator a = A.begin(); a!=A.end(); ++a )
      {
          const CoordSysEigen<T> shapeAtobodyA
              = CoordSysEigen<T>((a->transform().T()), (a->transform().Q()));

          const CoordSysEigen<T> shapeAtoWCS = prod(shapeAtobodyA, bodyAtoWCS);

          geometry::Sphere<typename V::real_type> const sphere
              = geometry::make_sphere((shapeAtoWCS.T()), a->radius());

          kdop::single_traversal<8, T>(
              sphere, objB.m_tree, geoB.m_tetramesh.m_mesh, objB.m_X, objB.m_Y,
              objB.m_Z, geoB.m_tetramesh.m_surface_map, callback, should_flip);
      }
  }

  } // namespace detail

} //namespace narrow

// NARROW_SPHERES_TETRAMESH_H
#endif
