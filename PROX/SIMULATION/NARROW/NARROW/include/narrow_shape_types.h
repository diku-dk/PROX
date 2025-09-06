#ifndef NARROW_SHAPE_TYPES_H
#define NARROW_SHAPE_TYPES_H

#include <convex_shapes.h> // Needed for convex::ConvexHull data type
#include <mesh_array.h>    // needed for mesh_array data types
#include <kdop.h>

#include <tiny_vector_functions.h>  // needed for tiny::norm

namespace narrow
{
  namespace detail
  {

    template<typename M>
    class ShapeTypes
    {
    protected:

      typedef typename M::coordsys_type C;
      typedef typename M::real_type     T;
      typedef typename M::vector3_type  V;
      typedef typename M::value_traits  VT;

    protected:

      class ShapeBase
      {

      protected:

        C m_transform; ///< A coordinate transformation from the local shape frame to the local frame of the object.

      public:

        /**
         * Shape Transformation Accessor.
         * The shape transformation places the shape in the local object frame.
         *
         * @return     A reference to the current shape transformation
         */
        C const & transform(  ) const {  return m_transform; }
        C& transform() { return m_transform; }

        auto eigenTransform() const { return coordSysToEigen(m_transform); }

        /**
         * Get a support map function of the shape.
         *
         * @return    A pointer to a support map functor.
         */

      public:

        ShapeBase( )
        : m_transform( C::identity() )
        {}


      };

    public:

      class Sphere : public ShapeBase
      {
      protected:

        // 2015-02-01 Kenny code review: Why not use convex::Sphere (or geometry::Sphere) as a member? Like for ConvexHull?

        T m_radius;

      public:

        Sphere() {}
        virtual ~Sphere() {}

      public:

        T const & radius() const { return m_radius; }
        T       & radius()       { return m_radius; }
        T const & scale()  const { return m_radius; }

      };

      class Box : public ShapeBase
      {
      protected:

        // 2015-02-01 Kenny code review: Why not use convex::Box (or geometry::Box) as a member? Like for ConvexHull?

        V m_half_ext;

      public:

        Box() {}

        virtual ~Box() {}

      public:

        V const & half_extent() const { return m_half_ext;      }
        V       & half_extent()       { return m_half_ext;      }
        T         scale()       const { return min(m_half_ext); }

      };

      class ConvexHull : public ShapeBase
      {
      public:

          typedef convex::ConvexHull<typename M::real_type> data_type;

      protected:

        data_type m_data;    ///< Use ConvexHull data representation from the CONVEX library

      public:

        ConvexHull() {}

        virtual ~ConvexHull() {}

      public:

        data_type const & data()  const { return m_data;             }
        data_type       & data()        { return m_data;             }
        T                 scale() const { return m_data.get_scale(); }
      };

      class Tetramesh
      {
      public:

        mesh_array::T4Mesh      m_mesh;   ///< Volume mesh of object

        mesh_array::VertexAttribute<T,mesh_array::T4Mesh> m_X0;     ///< Undeformed (material) x-coordinate
        mesh_array::VertexAttribute<T,mesh_array::T4Mesh> m_Y0;     ///< Undeformed (material) y-coordinate
        mesh_array::VertexAttribute<T,mesh_array::T4Mesh> m_Z0;     ///< Undeformed (material) z-coordinate

        mesh_array::TetrahedronAttribute<mesh_array::TetrahedronSurfaceInfo,mesh_array::T4Mesh> m_surface_map;

        T                       m_mesh_radius;
        T                       m_mesh_scale;

      public:

        Tetramesh()
        : m_mesh()
        , m_X0()
        , m_Y0()
        , m_Z0()
        , m_surface_map()
        , m_mesh_radius( 0 )
        , m_mesh_scale(0 )
        {}

        virtual ~Tetramesh() {}

      public:

        T scale() const { return m_mesh_scale; }

        void set_tetramesh_shape( mesh_array::T4Mesh const & mesh
                                 , mesh_array::VertexAttribute<T,mesh_array::T4Mesh> const & X
                                 , mesh_array::VertexAttribute<T,mesh_array::T4Mesh> const & Y
                                 , mesh_array::VertexAttribute<T,mesh_array::T4Mesh> const & Z
                                 )
        {
          using tiny::norm;
          using std::max;
          using tiny::min;
          using tiny::max;

          kdop::mesh_reorder( mesh
                             , X
                             , Y
                             , Z
                             , m_mesh
                             , m_X0
                             , m_Y0
                             , m_Z0
                             );

          mesh_array::compute_surface_map( m_mesh
                                    , m_X0
                                    , m_Y0
                                    , m_Z0
                                    , m_surface_map
                                    );

          m_mesh_radius = 0;

          size_t const N = m_mesh.vertex_size();


          V min_coord = V(std::numeric_limits<T>::max());
          V max_coord = V(std::numeric_limits<T>::lowest());

          for(size_t n = 0u; n < N;++n)
          {
            mesh_array::Vertex const & v = m_mesh.vertex(n);

            V const r0 = V::make(  m_X0(v), m_Y0(v), m_Z0(v) );

            min_coord = min(min_coord,r0);
            max_coord = max(max_coord,r0);

            m_mesh_radius = max( m_mesh_radius, norm(r0) );
          }

          m_mesh_scale = min(  max_coord-min_coord );
        }


        bool has_data() const
        {
          return ((m_mesh.vertex_size() > 0u) && (m_mesh.tetrahedron_size() > 0u));
        }

        void clear()
        {
          m_X0.release();
          m_Y0.release();
          m_Z0.release();
          m_surface_map.release();

          m_mesh.clear();
          m_mesh_radius = 0;
          m_mesh_scale  = 0;
        }


      };

    };

  } // namespace detail

  } // namespace narrow

  //==============================================================================
  namespace narrow
  {
  namespace detail
  {

  template <typename T> class ShapeTypesEigen
  {
  protected:
  protected:
      class ShapeBaseEigen
      {

      protected:
          CoordSysEigen<T>
              m_transform; ///< A coordinate transformation from the local shape frame to the local frame of the object.

      public:
          /**
         * Shape Transformation Accessor.
         * The shape transformation places the shape in the local object frame.
         *
         * @return     A reference to the current shape transformation
         */
          const CoordSysEigen<T>& transform() const { return m_transform; }

          CoordSysEigen<T>& transform() { return m_transform; }

          /**
         * Get a support map function of the shape.
         *
         * @return    A pointer to a support map functor.
         */

      public:
          ShapeBaseEigen()
              : m_transform(CoordSysEigen<T>::identity())
          {
          }
      };

  public:
      class Sphere : public ShapeBaseEigen
      {
      protected:
          // 2015-02-01 Kenny code review: Why not use convex::Sphere (or geometry::Sphere) as a member? Like for ConvexHull?

          T m_radius;

      public:
          Sphere() {}

          virtual ~Sphere() {}

      public:
          T const& radius() const { return m_radius; }

          T& radius() { return m_radius; }

          T const& scale() const { return m_radius; }
      };

      class Box : public ShapeBaseEigen
      {
      protected:
          // 2015-02-01 Kenny code review: Why not use convex::Box (or geometry::Box) as a member? Like for ConvexHull?

          EigenVector3<T> m_half_ext;

      public:
          Box() {}

          virtual ~Box() {}

      public:
          const EigenVector3<T>& half_extent() const { return m_half_ext; }

          EigenVector3<T>& half_extent() { return m_half_ext; }

          T scale() const { return min(m_half_ext); }
      };

      class ConvexHull : public ShapeBaseEigen
      {
      public:
          typedef convex::ConvexHull<T> data_type;

      protected:
          data_type m_data; ///< Use ConvexHull data representation from the CONVEX library

      public:
          ConvexHull() {}

          virtual ~ConvexHull() {}

      public:
          data_type const& data() const { return m_data; }

          data_type& data() { return m_data; }

          T scale() const { return m_data.get_scale(); }
      };

      class Tetramesh
      {
      public:
          mesh_array::T4Mesh m_mesh; ///< Volume mesh of object

          mesh_array::VertexAttribute<T, mesh_array::T4Mesh> m_X0; ///< Undeformed (material) x-coordinate
          mesh_array::VertexAttribute<T, mesh_array::T4Mesh> m_Y0; ///< Undeformed (material) y-coordinate
          mesh_array::VertexAttribute<T, mesh_array::T4Mesh> m_Z0; ///< Undeformed (material) z-coordinate

          mesh_array::TetrahedronAttribute<mesh_array::TetrahedronSurfaceInfo, mesh_array::T4Mesh> m_surface_map;

          T m_mesh_radius;
          T m_mesh_scale;

      public:
          Tetramesh()
              : m_mesh()
              , m_X0()
              , m_Y0()
              , m_Z0()
              , m_surface_map()
              , m_mesh_radius(0)
              , m_mesh_scale(0)
          {
          }

          virtual ~Tetramesh() {}

      public:
          T scale() const { return m_mesh_scale; }

          void set_tetramesh_shape(mesh_array::T4Mesh const& mesh,
                                   mesh_array::VertexAttribute<T, mesh_array::T4Mesh> const& X,
                                   mesh_array::VertexAttribute<T, mesh_array::T4Mesh> const& Y,
                                   mesh_array::VertexAttribute<T, mesh_array::T4Mesh> const& Z)
          {
              using std::max;
              using tiny::max;
              using tiny::min;
              using tiny::norm;

              kdop::mesh_reorder(mesh, X, Y, Z, m_mesh, m_X0, m_Y0, m_Z0);

              mesh_array::compute_surface_map(m_mesh, m_X0, m_Y0, m_Z0, m_surface_map);

              m_mesh_radius = 0;

              size_t const N = m_mesh.vertex_size();

              EigenVector3<T> min_coord = EigenVector3<T>(std::numeric_limits<T>::max());
              EigenVector3<T> max_coord = EigenVector3<T>(std::numeric_limits<T>::lowest());

              for (size_t n = 0u; n < N; ++n)
              {
                  mesh_array::Vertex const& v = m_mesh.vertex(n);

                  const EigenVector3<T> r0 = EigenVector3<T>(m_X0(v), m_Y0(v), m_Z0(v));

                  min_coord = min(min_coord, r0);
                  max_coord = max(max_coord, r0);

                  m_mesh_radius = max(m_mesh_radius, norm(r0));
              }

              m_mesh_scale = min(max_coord - min_coord);
          }

          bool has_data() const { return ((m_mesh.vertex_size() > 0u) && (m_mesh.tetrahedron_size() > 0u)); }

          void clear()
          {
              m_X0.release();
              m_Y0.release();
              m_Z0.release();
              m_surface_map.release();

              m_mesh.clear();
              m_mesh_radius = 0;
              m_mesh_scale = 0;
          }
      };
  };

  } // namespace detail

  } // namespace narrow

// NARROW_SHAPE_TYPES_H
#endif
