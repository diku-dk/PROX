#ifndef MESH_ARRAY_MAKE_CYLINDER_H
#define MESH_ARRAY_MAKE_CYLINDER_H

#include <mesh_array_t3mesh.h>
#include <mesh_array_vertex_attribute.h>

namespace mesh_array
{

template<typename T>
void make_cylinder(
    T const & radius
    , T const & height
    , size_t const & slices
    , T3Mesh & mesh
    , VertexAttribute<T,T3Mesh> & X
    , VertexAttribute<T,T3Mesh> & Y
    , VertexAttribute<T,T3Mesh> & Z
    );

} //namespace mesh_array

// MESH_ARRAY_MAKE_CYLINDER_H
#endif
