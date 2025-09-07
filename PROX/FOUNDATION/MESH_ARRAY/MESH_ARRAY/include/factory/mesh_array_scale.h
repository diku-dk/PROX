#ifndef MESH_ARRAY_SCALE_H
#define MESH_ARRAY_SCALE_H

#include <mesh_array_t3mesh.h>
#include <mesh_array_t4mesh.h>
#include <mesh_array_vertex_attribute.h>

namespace mesh_array
{

  template <typename T>
  inline void scale(const T& sx, const T& sy, const T& sz, T3Mesh const& mesh, VertexAttribute<T, T3Mesh>& X,
                    VertexAttribute<T, T3Mesh>& Y, VertexAttribute<T, T3Mesh>& Z)
  {

      for (size_t i = 0u; i < mesh.vertex_size(); ++i)
      {
          Vertex const v = mesh.vertex(i);

          const EigenVector3<T> p = EigenVector3<T>(sx * X(v), sy * Y(v), sz * Z(v));

          X(v) = p(0);
          Y(v) = p(1);
          Z(v) = p(2);
      }
  }

  template <typename T>
  inline void scale(const T& sx, const T& sy, const T& sz, T4Mesh const& mesh, VertexAttribute<T, T4Mesh>& X,
                    VertexAttribute<T, T4Mesh>& Y, VertexAttribute<T, T4Mesh>& Z)
  {

      for (size_t i = 0u; i < mesh.vertex_size(); ++i)
      {
          Vertex const v = mesh.vertex(i);

          const EigenVector3<T> p = EigenVector3<T>(sx * X(v), sy * Y(v), sz * Z(v));

          X(v) = p(0);
          Y(v) = p(1);
          Z(v) = p(2);
      }
  }

} //namespace mesh_array

// MESH_ARRAY_SCALE_H
#endif
