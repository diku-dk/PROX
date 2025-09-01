#ifndef MESH_ARRAY_MAKE_TETRAHEDRON_H
#define MESH_ARRAY_MAKE_TETRAHEDRON_H

#include <mesh_array_t3mesh.h>
#include <mesh_array_vertex_attribute.h>
#include <tiny_math_types.h>

namespace mesh_array
{

template<typename T>
void make_tetrahedron(
    EigenVector3<T> const & p0
    , EigenVector3<T> const & p1
    , EigenVector3<T> const & p2
    , EigenVector3<T> const & p3
    , T3Mesh & mesh
    , VertexAttribute<T,T3Mesh> & X
    , VertexAttribute<T,T3Mesh> & Y
    , VertexAttribute<T,T3Mesh> & Z
    );

} //namespace mesh_array

// MESH_ARRAY_MAKE_TETRAHEDRON_H
#endif
