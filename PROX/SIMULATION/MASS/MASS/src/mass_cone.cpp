#include <mass.h>

#include <eigenhelperall.h>

#include <cassert>

namespace mass
{

template <typename T> Properties<T> compute_cone(T const& density, T const& base_radius, T const& height)
{

    T const pi = std::numbers::pi_v<T>;
    T const k1 = T(1.0 / 10.0);
    T const k2 = T(3.0 / 20.0);
    T const k3 = T(3.0 / 10.0);

    assert(is_number(density) || !"density must be a number");
    assert(is_finite(density) || !"density must be a finite number");
    assert(density > 0 || !"density must be positive");
    assert(is_number(base_radius) || !"base_radius must be a  number");
    assert(is_finite(base_radius) || !"base_radius must be a finite number");
    assert(base_radius > 0 || !"base_radius must be positive");
    assert(is_number(height) || !"height must be a number");
    assert(is_finite(height) || !"height must be a finite number");
    assert(height > 0 || !"height must be positive");

    T const R = base_radius;
    T const h = height;
    T const volume = pi * R * R * h / 3;
    T const mass = density * volume;

    Properties<T> value;

    value.m_m = mass;
    value.m_Ixx = k1 * mass * h * h + k2 * mass * R * R;
    value.m_Iyy = k3 * mass * R * R;
    value.m_Izz = value.m_Ixx;
    value.m_y = h / 4;

    return value;
}

template Properties<float> compute_cone<float>(float const& density, float const& base_radius, float const& height);

template Properties<double> compute_cone<double>(double const& density, double const& base_radius,
                                                 double const& height);

} // namespace mass
