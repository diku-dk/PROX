#ifndef MESH_ARRAY_MAKE_ELLIPSOID_H
#define MESH_ARRAY_MAKE_ELLIPSOID_H

#include "mesh_array_make_sphere.h"
#include "tiny_math_types.h"
#include <mesh_array_t3mesh.h>
#include <mesh_array_vertex_attribute.h>
#include <factory/mesh_array_scale.h>

namespace mesh_array
{

template <typename Number>
void make_ellipsoid(Number a, Number b, Number c,
                    size_t slices, size_t segments,
                    T3Mesh& mesh,
                    VertexAttribute<Number, T3Mesh>& xAttr,
                    VertexAttribute<Number, T3Mesh>& yAttr,
                    VertexAttribute<Number, T3Mesh>& zAttr)
{
    make_sphere<Number>(1, slices, segments, mesh, xAttr, yAttr, zAttr);
    scale<Number>(a, b, c, mesh, xAttr, yAttr, zAttr);
}

} //namespace mesh_array

// MESH_ARRAY_MAKE_ELLIPSOID_H
#endif
