#ifndef MESH_ARRAY_PROFILE_SWEEP_H
#define MESH_ARRAY_PROFILE_SWEEP_H

#include <mesh_array_t3mesh.h>
#include <mesh_array_vertex_attribute.h>
#include <tiny_math_types.h>

namespace mesh_array
{

  /**
   * Sweep a profile around the y-axis to create a full closed mesh surface.
   */
  template<typename T>
  void profile_sweep(
                     std::vector<EigenVector3<T>> const & profile
                     , size_t const & slices
                     , T3Mesh & mesh
                     , VertexAttribute<T,T3Mesh> & X
                     , VertexAttribute<T,T3Mesh> & Y
                     , VertexAttribute<T,T3Mesh> & Z
                     );

} //namespace mesh_array

// MESH_ARRAY_PROFILE_SWEEP_H
#endif
