#include <factory/mesh_array_make_cone.h>
#include <factory/mesh_array_profile_sweep.h>

#include <eigenhelperall.h>

#include <vector>

namespace mesh_array
{

template <typename T>
void make_cone(T const& radius, T const& height, size_t const& slices, T3Mesh& mesh, VertexAttribute<T, T3Mesh>& X,
               VertexAttribute<T, T3Mesh>& Y, VertexAttribute<T, T3Mesh>& Z)
{

    std::vector<EigenVector3<T>> profile;

    profile.resize(3u);

    profile[0] = EigenVector3<T>(0, 0, 0);
    profile[1] = EigenVector3<T>(radius, 0, 0);
    profile[2] = EigenVector3<T>(0, height, 0);

    profile_sweep<T>(profile, slices, mesh, X, Y, Z);
}

template void make_cone<float>(float const& radius, float const& height, size_t const& slices, T3Mesh& mesh,
                               VertexAttribute<float, T3Mesh>& X, VertexAttribute<float, T3Mesh>& Y,
                               VertexAttribute<float, T3Mesh>& Z);

template void make_cone<double>(double const& radius, double const& height, size_t const& slices, T3Mesh& mesh,
                                VertexAttribute<double, T3Mesh>& X, VertexAttribute<double, T3Mesh>& Y,
                                VertexAttribute<double, T3Mesh>& Z);

} //namespace mesh_array
