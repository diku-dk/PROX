#ifndef MESH_ARRAY_MAKE_CUBOID_H
#define MESH_ARRAY_MAKE_CUBOID_H

#include <mesh_array_t3mesh.h>
#include <mesh_array_vertex_attribute.h>
#include <tiny_math_types.h>

namespace mesh_array
{

  /**
   * The vertices of the eight corners of a cuboid are specified.
   * The first four vertices are the front face in CCW order. The last four gives
   * the back face vertices in same order as the front face. That is in CW order if
   * one is looking at the front-side of the back face.
   */
template<typename T>
void make_cuboid(EigenVector3<T> const & v0
                 , EigenVector3<T> const & v1
                 , EigenVector3<T> const & v2
                 , EigenVector3<T> const & v3
                 , EigenVector3<T> const & v4
                 , EigenVector3<T> const & v5
                 , EigenVector3<T> const & v6
                 , EigenVector3<T> const & v7
                 , T3Mesh & mesh
                 , VertexAttribute<T,T3Mesh> & X
                 , VertexAttribute<T,T3Mesh> & Y
                 , VertexAttribute<T,T3Mesh> & Z
                 );

} //namespace mesh_array

// MESH_ARRAY_MAKE_CUBOID_H
#endif
