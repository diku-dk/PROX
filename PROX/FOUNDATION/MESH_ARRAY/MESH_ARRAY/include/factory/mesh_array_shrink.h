#ifndef MESH_ARRAY_SHRINK_H
#define MESH_ARRAY_SHRINK_H

#include <mesh_array_t3mesh.h>
#include <mesh_array_t4mesh.h>
#include <mesh_array_vertex_attribute.h>

#include <factory/mesh_array_translate.h>
#include <factory/mesh_array_scale.h>
#include <factory/mesh_array_compute_center.h>

#include <tiny_math_types.h>

namespace mesh_array
{

template <typename T>
inline void shrink(const T& factor, T3Mesh const& mesh,
                   VertexAttribute<T, T3Mesh>& X, VertexAttribute<T, T3Mesh>& Y,
                   VertexAttribute<T, T3Mesh>& Z)
{
    const EigenVector3<T> d = compute_center<T>(mesh, X, Y, Z);

    translate<T>(-d, mesh, X, Y, Z);
    scale<T>(factor, factor, factor, mesh, X, Y, Z);
    translate<T>(d, mesh, X, Y, Z);
}

template <typename T>
inline void shrink(const T& factor, T4Mesh const& mesh,
                   VertexAttribute<T, T4Mesh>& X, VertexAttribute<T, T4Mesh>& Y,
                   VertexAttribute<T, T4Mesh>& Z)
{
    const EigenVector3<T> d = (compute_center<T>(mesh, X, Y, Z));

    translate<T>(-d, mesh, X, Y, Z);
    scale<T>(factor, factor, factor, mesh, X, Y, Z);
    translate<T>(d, mesh, X, Y, Z);
}

} //namespace mesh_array

// MESH_ARRAY_SHRINK_H
#endif
