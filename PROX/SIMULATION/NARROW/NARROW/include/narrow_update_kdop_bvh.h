#ifndef NARROW_UPDATE_KDOP_BVH_H
#define NARROW_UPDATE_KDOP_BVH_H

#include <kdop_refit_tree.h>

#include <mesh_array.h>

#include <narrow_geometry.h>
#include <narrow_object.h>
#include <narrow_tags.h>

#include <eigenhelperall.h>

#include <util_profiling.h>

#include <cassert>
#include <vector>

namespace narrow
{

/**
   * This class encapsulates all information necessary for batch processing
   * of KDOP-BVH updates of objects.
   *
   * @tparam M a typebinder of math types.
   */
template <typename T>
requires std::is_floating_point_v<T>
class KDopBvhUpdateWorkItem
{
protected:
    Object<T>* m_object;
    Geometry<T> const* m_geometry;
    EigenVector3<T> m_p;
    EigenQuaternion<T> m_q;

private:
    bool is_valid() const
    {
        if (this->m_object == 0) return false;
        if (this->m_geometry == 0) return false;
        return true;
    }

public:
    auto& object() const
    {
        assert(this->is_valid()
               || "KDopBvhUpdateWorkItem::object(): null pointer");

        return *(this->m_object);
    }

    auto const& geometry() const
    {
        assert(this->is_valid()
               || "KDopBvhUpdateWorkItem::geometry() : null pointer");

        return *(this->m_geometry);
    }

    auto const& p() const
    {
        assert(this->is_valid() || "KDopBvhUpdateWorkItem::p(): null pointer");

        return this->m_p;
    }

    auto const& q() const
    {
        assert(this->is_valid() || "KDopBvhUpdateWorkItem::q(): null pointer");

        return this->m_q;
    }

public:
    KDopBvhUpdateWorkItem()
        : m_object(0)
        , m_geometry(0)
        , m_p()
        , m_q()
    {
    }

    ~KDopBvhUpdateWorkItem() {}

    KDopBvhUpdateWorkItem(Object<T>& object, Geometry<T> const& geometry,
                          EigenVector3<T> const& p, EigenQuaternion<T> const& q)
        : m_object(&object)
        , m_geometry(&geometry)
        , m_p(p)
        , m_q(q)
    {
    }

    KDopBvhUpdateWorkItem(KDopBvhUpdateWorkItem<T> const& item)
    {
        *this = item;
    }

    KDopBvhUpdateWorkItem<T>& operator=(KDopBvhUpdateWorkItem<T> const& item)
    {
        if (this != &item)
        {
            this->m_object = item.m_object;
            this->m_geometry = item.m_geometry;
            this->m_p = item.m_p;
            this->m_q = item.m_q;
        }
        return *this;
    }
};

template <typename T>
inline void update_kdop_bvh(std::vector<KDopBvhUpdateWorkItem<T>>& work_pool,
                            sequential const& /* tag */
)
{
    assert(!work_pool.empty()
           || "update_kdop_bvh : update_kdop_bvh_objects are empty");

    START_TIMER("refit_tree");

    for (auto& current : work_pool)
    {
        auto& object = current.object();
        auto const& geometry = current.geometry();
        auto const p = current.p();
        auto const q = current.q();
        size_t const N = geometry.m_tetramesh.m_mesh.vertex_size();

        if (N <= 0u) continue;

        for (size_t n = 0u; n < N; ++n)
        {
            mesh_array::Vertex const& v = geometry.m_tetramesh.m_mesh.vertex(n);

            EigenVector3<T> r0(geometry.m_tetramesh.m_X0(v),
                               geometry.m_tetramesh.m_Y0(v),
                               geometry.m_tetramesh.m_Z0(v));

            auto r = (rotate(q, r0) + p).eval();

            object.m_X(v) = r(0);
            object.m_Y(v) = r(1);
            object.m_Z(v) = r(2);
        }

        kdop::refit_tree<8, T>(object.m_tree, geometry.m_tetramesh.m_mesh,
                               object.m_X, object.m_Y, object.m_Z,
                               kdop::sequential());
    }

    STOP_TIMER("refit_tree");
}

} // namespace narrow

// NARROW_UPDATE_KDOP_BVH_H
#endif
