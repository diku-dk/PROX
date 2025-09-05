#ifndef PROX_FORCE_CALLBACKS_H
#define PROX_FORCE_CALLBACKS_H

#include <prox_rigid_body.h>

#include <cassert>

#include <type_traits>

namespace prox
{
template <typename T>
requires(std::is_floating_point_v<T>)
class ForceCallback
{
protected:
    size_t m_idx;

public:
    size_t const& get_idx() const { return this->m_idx; }

    void set_idx(size_t const& idx) { this->m_idx = idx; }

public:
    ForceCallback()
        : m_idx(0u)
    {
    }

public:
    virtual void compute_force_and_torque(const RigidBody<T>& body, EigenVector3<T>& force,
                                          EigenVector3<T>& torque) const
        = 0;
};

template <typename T>
requires(std::is_floating_point_v<T>)
class Gravity : public ForceCallback<T>
{
protected:
    T m_acceleration;
    EigenVector3<T> m_up;

public:
    T const& acceleration() const { return this->m_acceleration; }

    T& acceleration() { return this->m_acceleration; }

    EigenVector3<T> const& up() const { return this->m_up; }

    EigenVector3<T>& up() { return this->m_up; }

public:
    Gravity()
        : m_acceleration(9.81)
        , m_up(0, 1, 0)
    {
    }

public:
    void compute_force_and_torque(const RigidBody<T>& body, EigenVector3<T>& force,
                                  EigenVector3<T>& torque) const
    {
        assert(this->m_acceleration >= 0
               || !"Gravity::compute_force_and_torque(): acceleration must be non-negative");

        force = -m_up * (body.get_mass() * this->m_acceleration);
        torque = 0;
    }
};

template <typename T>
requires(std::is_floating_point_v<T>)
class Damping : public ForceCallback<T>
{
protected:
    T m_linear;
    T m_angular;

public:
    T const& linear() const { return this->m_linear; }

    T& linear() { return this->m_linear; }

    T const& angular() const { return this->m_angular; }

    T& angular() { return this->m_angular; }

public:
    Damping()
        : m_linear(0.01)
        , m_angular(0.01)
    {
    }

public:
    void compute_force_and_torque(const RigidBody<T>& body, EigenVector3<T>& force,
                                  EigenVector3<T>& torque) const
    {
        assert(this->m_linear >= 0
               || !"Damping::compute_force_and_torque(): damping must be non-negative");
        assert(this->m_angular >= 0
               || !"Damping::compute_force_and_torque(): damping must be non-negative");

        force = -(body.get_velocity() * m_linear);
        torque = -(body.get_spin() * m_angular);
    }
};

template <typename T>
requires(std::is_floating_point_v<T>)
class Pin : public ForceCallback<T>
{
protected:
    T m_tau; ///< Characteristic damping time, use 1/6 of time-step size, usually gives nice results?
    EigenVector3<T> m_target; ///< Target position that must be hit
    EigenVector3<T>
        m_anchor; ///< Anchor point in local body space coordinates. This the point on the body that should be pinned to target.

public:
    T const& tau() const { return this->m_tau; }

    T& tau() { return this->m_tau; }

    EigenVector3<T> const& target() const { return this->m_target; }

    EigenVector3<T>& target() { return this->m_target; }

    EigenVector3<T> const& anchor() const { return this->m_anchor; }

    EigenVector3<T>& anchor() { return this->m_anchor; }

public:
    Pin()
        : m_tau(0.1)
        , m_target(0, 0, 0)
        , m_anchor(0, 0, 0)
    {
    }

public:
    void compute_force_and_torque(const RigidBody<T>& body, EigenVector3<T>& force,
                                  EigenVector3<T>& torque) const
    {
        assert(this->m_tau > 0 || !"MoveToPoint::compute_force_and_torque(): tau must be positive");

        EigenVector3<T> const r = rotate(body.get_orientation(), this->m_anchor);
        EigenVector3<T> const D = this->m_target - (r + body.get_position());
        T const l = norm(D);
        EigenVector3<T> const n = D.normalized();
        T const v = body.get_velocity().dot(n);
        T const m = body.get_mass();
        T const b = (2 * m) / this->m_tau;
        T const k = m / (this->m_tau * this->m_tau);

        force = (k * l - b * v) * n;

        torque = cross(r, force);
    }
};

} //namespace prox

// PROX_FORCE_CALLBACKS_H
#endif
