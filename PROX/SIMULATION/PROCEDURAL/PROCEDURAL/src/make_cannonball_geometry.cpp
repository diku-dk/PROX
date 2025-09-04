#include <procedural_factory.h>

#include <tiny_math_types.h>

namespace procedural
{
template <typename T> GeometryHandleEigen<T> make_cannonball_geometry(content::API* engine, const T& radius)
{
    return create_geometry_handle_sphere<T>(engine, radius);
}

template GeometryHandleEigen<float> make_cannonball_geometry<float>(content::API* engine, const float& radius);
} // namespace procedural
