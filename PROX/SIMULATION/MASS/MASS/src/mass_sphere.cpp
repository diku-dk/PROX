#include <mass.h>

#include <eigenhelperall.h>

#include <cassert>

namespace mass
{

template <typename T> Properties<T> compute_sphere(T const& density, T const& radius)
{

    assert(is_number(density) || !"density must be a number");
    assert(is_finite(density) || !"density must be a finite number");
    assert(density > 0 || !"density must be positive");
    assert(is_number(radius) || !"radius must be a  number");
    assert(is_finite(radius) || !"radius must be a finite number");
    assert(radius > 0 || !"radius must be positive");

    T const volume = (4 * std::numbers::pi_v<T> * radius * radius * radius) / 3;
    T const mass = density * volume;
    T const five = T(5.0);

    Properties<T> value;

    value.m_m = mass;
    value.m_Ixx = (2 * mass * radius * radius) / five;
    value.m_Iyy = value.m_Ixx;
    value.m_Izz = value.m_Ixx;

    return value;
}

template Properties<float> compute_sphere<float>(float const& density, float const& radius);

template Properties<double> compute_sphere<double>(double const& density, double const& radius);

} // namespace mass
