#ifndef MESH_ARRAY_MAKE_CONVEX_H
#define MESH_ARRAY_MAKE_CONVEX_H

#include <mesh_array_t3mesh.h>
#include <mesh_array_vertex_attribute.h>

namespace mesh_array
{
template<typename T>
void make_convex(
    std::vector<EigenVector3<T>> const & vertices
    , T3Mesh & mesh
    , VertexAttribute<T,T3Mesh> & X
    , VertexAttribute<T,T3Mesh> & Y
    , VertexAttribute<T,T3Mesh> & Z
    )
;

} //namespace mesh_array

// MESH_ARRAY_MAKE_CONVEX_H
#endif
