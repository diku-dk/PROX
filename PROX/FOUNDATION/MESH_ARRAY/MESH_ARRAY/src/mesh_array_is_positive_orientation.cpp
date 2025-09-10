#include <mesh_array_is_positive_orientation.h>

#include <eigenhelperall.h>

#include <cassert>

namespace mesh_array
{

template <typename T>
bool is_positive_orientation(T4Mesh& mesh, VertexAttribute<T, T4Mesh>& X, VertexAttribute<T, T4Mesh>& Y,
                             VertexAttribute<T, T4Mesh>& Z)
{

    size_t const cntT = mesh.tetrahedron_size();

    for (size_t idx = 0u; idx < cntT; ++idx)
    {
        Tetrahedron const tetrahedron = mesh.tetrahedron(idx);

        Vertex const vi = mesh.vertex(tetrahedron.i());
        Vertex const vj = mesh.vertex(tetrahedron.j());
        Vertex const vk = mesh.vertex(tetrahedron.k());
        Vertex const vm = mesh.vertex(tetrahedron.m());

        const EigenVector3<T> pi = EigenVector3<T>(X(vi), Y(vi), Z(vi));
        EigenVector3<T> const pj = EigenVector3<T>(X(vj), Y(vj), Z(vj));
        EigenVector3<T> const pk = EigenVector3<T>(X(vk), Y(vk), Z(vk));
        EigenVector3<T> const pm = EigenVector3<T>(X(vm), Y(vm), Z(vm));

        EigenVector3<T> const e_ji = pj - pi;
        EigenVector3<T> const e_ki = pk - pi;
        EigenVector3<T> const e_mi = pm - pi;

        T const volume = (dot(e_mi, cross(e_ji, e_ki))) / T(6.0);

        if (volume <= 0) return false;
    }
    return true;
}

template bool is_positive_orientation<float>(T4Mesh& mesh, VertexAttribute<float, T4Mesh>& X,
                                             VertexAttribute<float, T4Mesh>& Y, VertexAttribute<float, T4Mesh>& Z);

template bool is_positive_orientation<double>(T4Mesh& mesh, VertexAttribute<double, T4Mesh>& X,
                                              VertexAttribute<double, T4Mesh>& Y, VertexAttribute<double, T4Mesh>& Z);

} // namespace mesh_array
