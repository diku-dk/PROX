#include <factory/mesh_array_make_tetrahedron.h>

#include <tiny_math_types.h>

namespace mesh_array
{

  template<typename T>
  void make_tetrahedron(
                   EigenVector3<T> const & p0
                   , EigenVector3<T> const & p1
                   , EigenVector3<T> const & p2
                   , EigenVector3<T> const & p3
                   , T3Mesh & mesh
                   , VertexAttribute<T,T3Mesh> & X
                   , VertexAttribute<T,T3Mesh> & Y
                   , VertexAttribute<T,T3Mesh> & Z
                   )
  {
    mesh.clear();

    mesh.set_capacity( 4u, 4u );

    X.bind(mesh);
    Y.bind(mesh);
    Z.bind(mesh);

    Vertex const vi = mesh.push_vertex();
    Vertex const vj = mesh.push_vertex();
    Vertex const vk = mesh.push_vertex();
    Vertex const vm = mesh.push_vertex();

    X(vi) = p0(0);    Y(vi) = p0(1);    Z(vi) = p0(2);
    X(vj) = p1(0);    Y(vj) = p1(1);    Z(vj) = p1(2);
    X(vk) = p2(0);    Y(vk) = p2(1);    Z(vk) = p2(2);
    X(vm) = p3(0);    Y(vm) = p3(1);    Z(vm) = p3(2);

        mesh.push_triangle( vi, vk, vj );
        mesh.push_triangle( vi, vj, vm );
        mesh.push_triangle( vj, vk, vm );
        mesh.push_triangle( vk, vi, vm );
    }


    template void make_tetrahedron<float>(EigenVector3<float> const& p0, EigenVector3<float> const& p1,
                                        EigenVector3<float> const& p2, EigenVector3<float> const& p3, T3Mesh& mesh,
                                        VertexAttribute<float, T3Mesh>& X,
                                        VertexAttribute<float, T3Mesh>& Y,
                                        VertexAttribute<float, T3Mesh>& Z);
    template void make_tetrahedron<double>(EigenVector3<double> const& p0, EigenVector3<double> const& p1,
                                        EigenVector3<double> const& p2, EigenVector3<double> const& p3, T3Mesh& mesh,
                                        VertexAttribute<double, T3Mesh>& X,
                                        VertexAttribute<double, T3Mesh>& Y,
                                        VertexAttribute<double, T3Mesh>& Z);

} //namespace mesh_array
