#include <factory/mesh_array_make_conical.h>
#include <factory/mesh_array_profile_sweep.h>

#include <tiny_math_types.h>

#include <vector>

namespace mesh_array
{

template <typename T>
void make_conical(T const& bottom_radius, T const& top_radius, T const& height, size_t const& slices, T3Mesh& mesh,
                  VertexAttribute<T, T3Mesh>& X, VertexAttribute<T, T3Mesh>& Y, VertexAttribute<T, T3Mesh>& Z)
{
    std::vector<EigenVector3<T>> profile;

    profile.resize(4u);

    profile[0] = EigenVector3<T>(0, 0, 0);
    profile[1] = EigenVector3<T>(bottom_radius, 0, 0);
    profile[2] = EigenVector3<T>(top_radius, height, 0);
    profile[3] = EigenVector3<T>(0, height, 0);

    profile_sweep<T>(profile, slices, mesh, X, Y, Z);
}

using MTf = tiny::MathTypes<float>;
using MTd = tiny::MathTypes<double>;

template void make_conical<float>(float const& bottom_radius, float const& top_radius, float const& height,
                                  size_t const& slices, T3Mesh& mesh, VertexAttribute<float, T3Mesh>& X,
                                  VertexAttribute<float, T3Mesh>& Y, VertexAttribute<float, T3Mesh>& Z);

template void make_conical<double>(double const& bottom_radius, double const& top_radius, double const& height,
                                   size_t const& slices, T3Mesh& mesh, VertexAttribute<double, T3Mesh>& X,
                                   VertexAttribute<double, T3Mesh>& Y, VertexAttribute<double, T3Mesh>& Z);

} //namespace mesh_array
