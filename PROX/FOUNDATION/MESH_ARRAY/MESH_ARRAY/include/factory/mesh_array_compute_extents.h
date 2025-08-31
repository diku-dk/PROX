#ifndef MESH_ARRAY_COMPUTE_EXTENTS_H
#define MESH_ARRAY_COMPUTE_EXTENTS_H

#include <mesh_array_t3mesh.h>
#include <mesh_array_t4mesh.h>
#include <mesh_array_vertex_attribute.h>

#include <limits>

namespace mesh_array
{

  template<typename T>
  inline void compute_extents(
                       T3Mesh const & mesh
                       , VertexAttribute<T,T3Mesh> const & C
                       , T& low
                       , T& hi
                       )
  {

    T h = std::numeric_limits<T>::lowest();
    T l = std::numeric_limits<T>::max();

    for(size_t i=0u; i < mesh.vertex_size(); ++i)
    {
      Vertex const v = mesh.vertex(i);
      h = std::max(h, C(v));
      l = std::min(l, C(v));
    }
    low = l;
    hi  = h;
  }


  template<typename T>
  inline void compute_extents(
                              T4Mesh const & mesh
                              , VertexAttribute<T,T4Mesh> const & C
                              , T& low
                              , T& hi
                              )
  {
    T h = std::numeric_limits<T>::lowest();
    T l = std::numeric_limits<T>::max();

    for(size_t i=0u; i < mesh.vertex_size(); ++i)
    {
      Vertex const v = mesh.vertex(i);
      h = std::max(h, C(v));
      l = std::min(l, C(v));
    }
    low = l;
    hi  = h;
  }


} //namespace mesh_array

// MESH_ARRAY_COMPUTE_CENTER_H
#endif
