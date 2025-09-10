#ifndef PROX_PROPERTY_H
#define PROX_PROPERTY_H

#include "eigenhelperfunctions.h"
#include <cstdlib>  // for std::size_t
#include <type_traits>

namespace prox
{

template < typename T >
requires std::is_floating_point_v<T>
class Property
{
public:
protected:
    EigenVector3<T> m_mu;
    T m_e;

    //given in local frame of the master matrial body
    EigenVector3<T> m_s;
    size_t m_master_material_idx;

public:
    Property()
        : m_mu(0.5f, 0.5f, 0.5f)
        , m_e(0)
        , m_s(1, 0, 0)
        , m_master_material_idx(0u)
    {}

    const T& get_restitution_coefficient() const { return this->m_e; }

    void set_restitution_coefficient(T e) { this->m_e = e; }

    auto const& get_friction_coefficients() const { return this->m_mu; }

    void get_friction_coefficients(T& mux, T& muy, T& muz) const
    {
      mux = this->m_mu[0];
      muy = this->m_mu[1];
      muz = this->m_mu[2];
    }

    void set_friction_coefficients(const EigenVector3<T>& mu) { this->m_mu = mu; }

    void set_friction_coefficients(T mux, T muy, T muz)
    {
      this->m_mu[0] = mux;
      this->m_mu[1] = muy;
      this->m_mu[2] = muz;
    }

    size_t       const & get_master_material_idx() const {  return this->m_master_material_idx; }

    void set_master_material_idx(size_t const & i)       {  this->m_master_material_idx = i; }

    const EigenVector3<T>& get_s_vector() const { return this->m_s; }

    void get_s_vector(T& sx, T& sy, T& sz) const
    {
      sx = this->m_s[0];
      sy = this->m_s[1];
      sz = this->m_s[2];
    }

    void set_s_vector(const EigenVector3<T>& s_) { this->m_s = s_; }

    void set_s_vector(T sx, T sy, T sz)
    {
      this->m_s[0] = sx;
      this->m_s[1] = sy;
      this->m_s[2] = sz;
    }

    bool is_isotropic() const
    {
      if (m_mu(0) != m_mu(1))
        return false;
      return true;
    }
};
} // namespace prox

// PROX_PROPERTY_H
#endif
