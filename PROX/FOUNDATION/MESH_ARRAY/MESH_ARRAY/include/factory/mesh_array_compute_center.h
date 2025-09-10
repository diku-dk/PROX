#ifndef MESH_ARRAY_COMPUTE_CENTER_H
#define MESH_ARRAY_COMPUTE_CENTER_H

#include <mesh_array_t3mesh.h>
#include <mesh_array_t4mesh.h>
#include <mesh_array_vertex_attribute.h>
#include <eigenhelperall.h>

namespace mesh_array
{

  template<typename T>
inline EigenVector3<T> compute_center(
                                           T3Mesh const & mesh
                                           , VertexAttribute<T,T3Mesh> const & X
                                           , VertexAttribute<T,T3Mesh> const & Y
                                           , VertexAttribute<T,T3Mesh> const & Z
                                           )
  {
    EigenVector3<T> center = {0,0,0};

    for(size_t i=0u; i < mesh.vertex_size(); ++i)
    {
      Vertex const v = mesh.vertex(i);
      center += EigenVector3<T>( X(v), Y(v), Z(v) );
    }

    center /= mesh.vertex_size();

    return center;
  }

  template<typename T>
  inline EigenVector3<T> compute_center(
                                                  T4Mesh const & mesh
                                                  , VertexAttribute<T,T4Mesh> const & X
                                                  , VertexAttribute<T,T4Mesh> const & Y
                                                  , VertexAttribute<T,T4Mesh> const & Z
                                                  )
  {

      EigenVector3<T> center = {0,0,0};

    for(size_t i=0u; i < mesh.vertex_size(); ++i)
    {
      Vertex const v = mesh.vertex(i);
      center += EigenVector3<T>( X(v), Y(v), Z(v) );
    }

    center /= mesh.vertex_size();

    return center;
  }




} //namespace mesh_array

// MESH_ARRAY_COMPUTE_CENTER_H
#endif
