#ifndef MESH_ARRAY_TRANSFORM_H
#define MESH_ARRAY_TRANSFORM_H

#include <mesh_array_t3mesh.h>
#include <mesh_array_t4mesh.h>
#include <mesh_array_vertex_attribute.h>
#include <eigenhelperall.h>

namespace mesh_array
{

  template<typename T>
  inline void transform(
      EigenVector3<T> const & d
      , EigenQuaternion<T> const & R
      , T3Mesh const & mesh
      , VertexAttribute<T,T3Mesh> & X
      , VertexAttribute<T,T3Mesh> & Y
      , VertexAttribute<T,T3Mesh> & Z
      )
  {

      for(size_t i=0u; i < mesh.vertex_size(); ++i)
      {
          Vertex const v = mesh.vertex(i);

          EigenVector3<T> vec( X(v), Y(v), Z(v) );
          EigenVector3<T> p = ::rotate(R, vec) + d;

          X(v) = p(0);
          Y(v) = p(1);
          Z(v) = p(2);
      }
  }


  template<typename T>
  inline void transform(
      EigenVector3<T> const & d
      , EigenQuaternion<T> const & R
      , T4Mesh const & mesh
      , VertexAttribute<T,T4Mesh> & X
      , VertexAttribute<T,T4Mesh> & Y
      , VertexAttribute<T,T4Mesh> & Z
      )
  {

      for(size_t i=0u; i < mesh.vertex_size(); ++i)
      {
          Vertex const v = mesh.vertex(i);

          const EigenVector3<T> p = rotate(R, EigenVector3<T>( X(v), Y(v), Z(v) ) ) + d;

          X(v) = p(0);
          Y(v) = p(1);
          Z(v) = p(2);
      }
  }


} //namespace mesh_array

// MESH_ARRAY_TRANSFORM_H
#endif
