#ifndef NARROW_SYSTEM_H
#define NARROW_SYSTEM_H

#include <narrow_geometry.h>
#include <narrow_object.h>
#include <narrow_params.h>

#include <cassert>
#include <vector>
#include <type_traits>

namespace narrow
{

template <typename T>
requires(std::is_floating_point_v<T>)
class System
{
protected:
    std::vector< Geometry<T> > m_geometries;
    Params<T> m_params;

public:
    typedef typename std::vector< Geometry<T> >::iterator geometry_iterator;

    geometry_iterator begin() { return m_geometries.begin(); }

    geometry_iterator end() { return m_geometries.end(); }

    const auto& params() const { return this->m_params; }

    auto& params() { return this->m_params; }

public:
    void clear() { m_geometries.clear(); }

    size_t size() const { return this->m_geometries.size(); }

    size_t create_geometry()
    {
        m_geometries.emplace_back();
        return m_geometries.size() - 1;
    }

    const auto& get_geometry(Object<T> const& obj) const { return m_geometries.at(obj.get_geometry_idx()); }

    auto& get_geometry(Object<T> const& obj) { return m_geometries.at(obj.get_geometry_idx()); }

    const auto& get_geometry(size_t const& idx) const { return m_geometries.at(idx); }

    auto& get_geometry(size_t const& idx) { return m_geometries.at(idx); }
};

} // namespace narrow

// NARROW_SYSTEM_H
#endif
