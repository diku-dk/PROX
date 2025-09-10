#include <factory/mesh_array_make_cuboid.h>

#include <eigenhelperall.h>

namespace mesh_array
{

template <typename T>
void make_cuboid(EigenVector3<T> const& v0, EigenVector3<T> const& v1, EigenVector3<T> const& v2,
                 EigenVector3<T> const& v3, EigenVector3<T> const& v4, EigenVector3<T> const& v5,
                 EigenVector3<T> const& v6, EigenVector3<T> const& v7, T3Mesh& mesh, VertexAttribute<T, T3Mesh>& X,
                 VertexAttribute<T, T3Mesh>& Y, VertexAttribute<T, T3Mesh>& Z)
{
    mesh.clear();
    X.release();
    Y.release();
    Z.release();

    mesh.set_capacity(8u, 12u);

    std::vector<Vertex> vertices;
    vertices.resize(8);

    vertices[0] = mesh.push_vertex();
    vertices[1] = mesh.push_vertex();
    vertices[2] = mesh.push_vertex();
    vertices[3] = mesh.push_vertex();
    vertices[4] = mesh.push_vertex();
    vertices[5] = mesh.push_vertex();
    vertices[6] = mesh.push_vertex();
    vertices[7] = mesh.push_vertex();

    unsigned int const quads[6][5] = {
        {0, 1, 2, 3},
        {4, 7, 6, 5},
        {4, 0, 3, 7},
        {1, 5, 6, 2},
        {7, 3, 2, 6},
        {5, 1, 0, 4},
    };

    for (size_t i = 0u; i < 6u; ++i)
    {
        Vertex const vi = mesh.vertex(quads[i][0]);
        Vertex const vj = mesh.vertex(quads[i][1]);
        Vertex const vk = mesh.vertex(quads[i][2]);
        Vertex const vm = mesh.vertex(quads[i][3]);

        mesh.push_triangle(vi, vj, vk);
        mesh.push_triangle(vi, vk, vm);
    }

    X.bind(mesh);
    Y.bind(mesh);
    Z.bind(mesh);

    X(vertices[0]) = v0(0);
    X(vertices[1]) = v1(0);
    X(vertices[2]) = v2(0);
    X(vertices[3]) = v3(0);
    X(vertices[4]) = v4(0);
    X(vertices[5]) = v5(0);
    X(vertices[6]) = v6(0);
    X(vertices[7]) = v7(0);

    Y(vertices[0]) = v0(1);
    Y(vertices[1]) = v1(1);
    Y(vertices[2]) = v2(1);
    Y(vertices[3]) = v3(1);
    Y(vertices[4]) = v4(1);
    Y(vertices[5]) = v5(1);
    Y(vertices[6]) = v6(1);
    Y(vertices[7]) = v7(1);

    Z(vertices[0]) = v0(2);
    Z(vertices[1]) = v1(2);
    Z(vertices[2]) = v2(2);
    Z(vertices[3]) = v3(2);
    Z(vertices[4]) = v4(2);
    Z(vertices[5]) = v5(2);
    Z(vertices[6]) = v6(2);
    Z(vertices[7]) = v7(2);
}

template void
make_cuboid<float>(EigenVector3<float> const& v0, EigenVector3<float> const& v1,
                   EigenVector3<float> const& v2, EigenVector3<float> const& v3,
                   EigenVector3<float> const& v4, EigenVector3<float> const& v5,
                   EigenVector3<float> const& v6, EigenVector3<float> const& v7,
                   T3Mesh& mesh, VertexAttribute<float, T3Mesh>& X,
                   VertexAttribute<float, T3Mesh>& Y,
                   VertexAttribute<float, T3Mesh>& Z);

template void make_cuboid<double>(
    EigenVector3<double> const& v0, EigenVector3<double> const& v1,
    EigenVector3<double> const& v2, EigenVector3<double> const& v3,
    EigenVector3<double> const& v4, EigenVector3<double> const& v5,
    EigenVector3<double> const& v6, EigenVector3<double> const& v7,
    T3Mesh& mesh, VertexAttribute<double, T3Mesh>& X,
    VertexAttribute<double, T3Mesh>& Y, VertexAttribute<double, T3Mesh>& Z);

} //namespace mesh_array
