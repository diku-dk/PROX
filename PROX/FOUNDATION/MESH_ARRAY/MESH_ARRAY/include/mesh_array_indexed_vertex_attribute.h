#ifndef MESH_ARRAY_INDEXED_VERTEX_ATTRIBUTE_H
#define MESH_ARRAY_INDEXED_VERTEX_ATTRIBUTE_H

#include <mesh_array_vertex_attribute.h>
#include <mesh_array_vertex.h>

#include <cstddef>

namespace mesh_array
{

  /**
   * Usefull for making a component of a vector3_type attribute behave like a "single" vertex attribute.
   */
  template<typename V, typename M>
  class IndexedVertexAttribute
  : public VertexAttribute<typename V::real_type, M>
  {
  protected:

    VertexAttribute<V,M> & m_data;
    size_t                 m_idx;

  private:

    IndexedVertexAttribute();  // Not allowed

  public:

    ~IndexedVertexAttribute();

    IndexedVertexAttribute(size_t const & idx, VertexAttribute<V,M> & data);

    IndexedVertexAttribute(IndexedVertexAttribute<V,M> const & o);

    IndexedVertexAttribute<V,M> const & operator=(IndexedVertexAttribute<V,M> const & o);

  public:

    M const * mesh() const;

    void bind(M const & m);

    void release();

    size_t size() const;

    typename V::real_type const & operator()(Vertex const & v) const;

    typename V::real_type & operator()(Vertex const & v);

    typename V::real_type const & operator[](Vertex const & v) const;

    typename V::real_type & operator[](Vertex const & v);

  };

} // end namespace mesh_array

  namespace mesh_array
  {

  /**
   * Usefull for making a component of a vector3_type attribute behave like a "single" vertex attribute.
   */
  template<typename T, typename M>
  class IndexedVertexAttributeEigen
      : public VertexAttribute<T, M>
  {
  public:

      VertexAttribute<T,M> & m_data;
      size_t                 m_idx;

  private:

      IndexedVertexAttributeEigen();  // Not allowed

  public:

      ~IndexedVertexAttributeEigen();

      IndexedVertexAttributeEigen(size_t const & idx, VertexAttribute<T,M> & data);

      IndexedVertexAttributeEigen(IndexedVertexAttributeEigen<T,M> const & o);

      IndexedVertexAttributeEigen<T,M> const & operator=(IndexedVertexAttributeEigen<T,M> const & o);

  public:

      M const * mesh() const;

      void bind(M const & m);

      void release();

      size_t size() const;

      const T& operator()(Vertex const & v) const;

      T& operator()(Vertex const & v);

      const T& operator[](Vertex const & v) const;

      T& operator[](Vertex const & v);

  };

  } // end namespace mesh_array

// MESH_ARRAY_INDEXED_VERTEX_ATTRIBUTE_H
#endif
