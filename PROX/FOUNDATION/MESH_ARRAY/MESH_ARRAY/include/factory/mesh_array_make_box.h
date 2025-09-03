#ifndef MESH_ARRAY_MAKE_BOX_H
#define MESH_ARRAY_MAKE_BOX_H

#include <mesh_array_t3mesh.h>
#include <mesh_array_vertex_attribute.h>

namespace mesh_array
{

  template<typename T>
  void make_box(
                const T& width
                , const T& height
                , const T& depth
                , T3Mesh & mesh
                , VertexAttribute<T,T3Mesh> & X
                , VertexAttribute<T,T3Mesh> & Y
                , VertexAttribute<T,T3Mesh> & Z
                );

  template<typename T>
  void make_box(
                const T& width
                , const T& height
                , const T& depth
                , T3Mesh & mesh
                , VertexAttribute<T,T3Mesh> & coords
                );

} //namespace mesh_array

// MESH_ARRAY_MAKE_BOX_H
#endif
