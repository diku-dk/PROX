#ifndef NARROW_TEST_PAIR_H
#define NARROW_TEST_PAIR_H

#include <narrow_system.h>
#include <narrow_object.h>

#include <geometry.h>

#include <tiny.h>
#include <cassert>

namespace narrow
{

/**
   * This class encapsulate information needed for the narrow phase collision
   * detection library to be able to deduce and dispatch the proper collision
   * handler for a pair of objects that needs to be tested for collision
   * detection.
   *
   * @tparam M   A typebinder of math types.
   */
template <typename T> class TestPair
{
public:
    using M = tiny::MathTypes<T>;

    typedef geometry::ContactsCallback<typename M::vector3_type> callback_type;

protected:
    Object<T> const* m_obj_a;
    Object<T> const* m_obj_b;

    EigenVector3<T> m_t_a;
    EigenVector3<T> m_t_b;

    EigenQuaternion<T> m_Q_a;
    EigenQuaternion<T> m_Q_b;

    callback_type* m_callback;

protected:
    bool validate() const
    {
        if (this->m_obj_a == 0) return false;
        if (this->m_obj_b == 0) return false;
        if (this->m_callback == 0) return false;
        return true;
    }

public:
    auto const& obj_a() const
    {
        assert(this->validate() || "TestPair::obj_a(): null pointer");

        return *(this->m_obj_a);
    }

    auto const& obj_b() const
    {
        assert(this->validate() || "TestPair::obj_b(): null pointer");

        return *(this->m_obj_b);
    }

    auto const& t_a() const
    {
        assert(this->validate() || "TestPair::t_a(): null pointer");

        return this->m_t_a;
    }

    auto const& t_b() const
    {
        assert(this->validate() || "TestPair::t_b(): null pointer");

        return this->m_t_b;
    }

    auto const& Q_a() const
    {
        assert(this->validate() || "TestPair::Q_a(): null pointer");

        return this->m_Q_a;
    }

    auto const& Q_b() const
    {
        assert(this->validate() || "TestPair::Q_b(): null pointer");

        return this->m_Q_b;
    }

    auto& callback()
    {
        assert(this->validate() || "TestPair::callback(): null pointer");

        return *(this->m_callback);
    }

public:
    TestPair()
        : m_obj_a(0)
        , m_obj_b(0)
        , m_t_a()
        , m_t_b()
        , m_Q_a()
        , m_Q_b()
        , m_callback(0)
    {
    }

    TestPair(Object<T> const& objA, Object<T> const& objB,
             const EigenVector3<T>& tA, const EigenQuaternion<T>& qA,
             const EigenVector3<T>& tB, const EigenQuaternion<T>& qB,
             geometry::ContactsCallback<typename M::vector3_type>& callback)
        : m_obj_a(&objA)
        , m_obj_b(&objB)
        , m_t_a(tA)
        , m_t_b(tB)
        , m_Q_a(qA)
        , m_Q_b(qB)
        , m_callback(&callback)
    {
    }
};

} //namespace narrow

// NARROW_TEST_PAIR_H
#endif
