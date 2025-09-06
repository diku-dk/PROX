#ifndef NARROW_OBJECT_H
#define NARROW_OBJECT_H

#include <narrow_geometry.h>
#include <narrow_params.h>

#include <kdop_make_tree.h>
#include <kdop_tags.h>
#include <kdop_tree.h>

#include <mesh_array.h>

#include <limits>

namespace narrow
{

/**
   * Narrow Phase Object Interface.
   * This class defines a common interface for all narrow phase collision objects.
   *
   * @tparam M    The math types
   */
template <typename T>
requires std::is_floating_point_v<T>
class Object
{
public:
    mesh_array::VertexAttribute<T, mesh_array::T4Mesh> m_X;      ///< Deformed (spatial) x-coordinate
    mesh_array::VertexAttribute<T, mesh_array::T4Mesh> m_Y;      ///< Deformed (spatial) y-coordinate
    mesh_array::VertexAttribute<T, mesh_array::T4Mesh> m_Z;      ///< Deformed (spatial) z-coordinate

    kdop::Tree<T, 8> m_tree; ///< kDOP BVH tree used for collision detection of deformed (spatial) mesh

protected:
    size_t m_geometry_idx; ///< A geometry index.

public:
    void set_geometry_idx(size_t const& idx) { this->m_geometry_idx = idx; }

    size_t const& get_geometry_idx() const { return this->m_geometry_idx; }

public:
    Object()
        : m_X()
        , m_Y()
        , m_Z()
        , m_tree()
        , m_geometry_idx(0u)
    {
    }

    ~Object() {}

    Object(Object const& obj) { *this = obj; }

    Object& operator=(Object const& obj)
    {
        if (this != &obj)
        {
            this->m_X = obj.m_X;
            this->m_Y = obj.m_Y;
            this->m_Z = obj.m_Z;
            this->m_tree = obj.m_tree;
            this->m_geometry_idx = obj.m_geometry_idx;
        }
        return *this;
    }
};

template <typename T> void make_kdop_bvh(const Params<T>& params, Object<T>& object, const Geometry<T>& geometry)
{
    if (geometry.m_tetramesh.has_data())
    {
        size_t mem_bytes;

#ifdef HAS_DIKUCL
        if (params.use_open_cl() && params.use_gproximity())
        {
            // make sure there is only ever one chunk per object
            mem_bytes = std::numeric_limits<std::size_t>::max();
        }
        else
        {
#endif // HAS_DIKUCL

            // 8000 bytes if gProximity cannot/should not be used
            mem_bytes = params.get_chunk_bytes();

#ifdef HAS_DIKUCL
        }
#endif // HAS_DIKUCL

        object.m_tree = kdop::make_tree<8, T>(mem_bytes, geometry.m_tetramesh.m_mesh, geometry.m_tetramesh.m_X0,
                                              geometry.m_tetramesh.m_Y0, geometry.m_tetramesh.m_Z0, kdop::sequential());

        object.m_X.bind(geometry.m_tetramesh.m_mesh);
        object.m_Y.bind(geometry.m_tetramesh.m_mesh);
        object.m_Z.bind(geometry.m_tetramesh.m_mesh);
    }
}

} //namespace narrow

// NARROW_OBJECT_H
#endif
