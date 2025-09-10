#include <mass.h>

#include <eigenhelperall.h>

#include <cassert>

namespace mass
{

template <typename T> Properties<T> compute_cylinder(T const& density, T const& radius, T const& half_height)
{

    assert(is_number(density) || !"density must be a number");
    assert(is_finite(density) || !"density must be a finite number");
    assert(density > 0 || !"density must be positive");
    assert(is_number(radius) || !"radius must be a  number");
    assert(is_finite(radius) || !"radius must be a finite number");
    assert(radius > 0 || !"radius must be positive");
    assert(is_number(half_height) || !"half_height must be a  number");
    assert(is_finite(half_height) || !"half_height must be a finite number");
    assert(half_height > 0 || !"half_height must be positive");

    T const& r = radius;
    T const& h = 2 * half_height;
    T const volume = std::numbers::pi_v<T> * r * r * h;
    T const mass = density * volume;

    Properties<T> value;

    value.m_m = mass;
    value.m_Ixx = (((r * r) / 4) + ((h * h) / 12)) * mass;
    value.m_Izz = value.m_Ixx;
    value.m_Iyy = mass * r * r * 0.5f;

    return value;
}

template Properties<float> compute_cylinder<float>(float const& density, float const& radius, float const& half_height);

template Properties<double> compute_cylinder<double>(double const& density, double const& radius,
                                                     double const& half_height);

} // namespace mass
