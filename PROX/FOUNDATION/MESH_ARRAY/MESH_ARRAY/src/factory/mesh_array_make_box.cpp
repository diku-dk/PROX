#include <mesh_array_indexed_vertex_attribute.h>

#include <factory/mesh_array_make_box.h>
#include <factory/mesh_array_make_cuboid.h>

#include <tiny_math_types.h>

namespace mesh_array
{

template <typename T>
void make_box(const T& width, const T& height, const T& depth, T3Mesh& mesh, VertexAttribute<T, T3Mesh>& X,
              VertexAttribute<T, T3Mesh>& Y, VertexAttribute<T, T3Mesh>& Z)
{
    const EigenVector3<T> v0 = EigenVector3<T>(-width, -height, depth) * 0.5f;
    const EigenVector3<T> v1 = EigenVector3<T>(width, -height, depth) * 0.5f;
    const EigenVector3<T> v2 = EigenVector3<T>(width, height, depth) * 0.5f;
    const EigenVector3<T> v3 = EigenVector3<T>(-width, height, depth) * 0.5f;
    const EigenVector3<T> v4 = v0 - EigenVector3<T>(0, 0, depth);
    const EigenVector3<T> v5 = v1 - EigenVector3<T>(0, 0, depth);
    const EigenVector3<T> v6 = v2 - EigenVector3<T>(0, 0, depth);
    const EigenVector3<T> v7 = v3 - EigenVector3<T>(0, 0, depth);

    make_cuboid<T>((v0), (v1), (v2), (v3), (v4), (v5), (v6), (v7), mesh, X, Y, Z);
}

  /*template<typename T>
  void make_box(
                const T& width
                , const T& height
                , const T& depth
                , T3Mesh & mesh
                , VertexAttribute<T,T3Mesh> & coords
                )
  {
    IndexedVertexAttributeEigen<T,T3Mesh> X = IndexedVertexAttributeEigen<T,T3Mesh>(0, coords);
    IndexedVertexAttributeEigen<T,T3Mesh> Y = IndexedVertexAttributeEigen<T,T3Mesh>(1, coords);
    IndexedVertexAttributeEigen<T,T3Mesh> Z = IndexedVertexAttributeEigen<T,T3Mesh>(2, coords);

    make_box<T>(width,height,depth,mesh,X,Y,Z);
  }*/

template void make_box<float>(float const& width, float const& height, float const& depth, T3Mesh& mesh,
                              VertexAttribute<float, T3Mesh>& X, VertexAttribute<float, T3Mesh>& Y,
                              VertexAttribute<float, T3Mesh>& Z);

template void make_box<double>(double const& width, double const& height, double const& depth, T3Mesh& mesh,
                               VertexAttribute<double, T3Mesh>& X, VertexAttribute<double, T3Mesh>& Y,
                               VertexAttribute<double, T3Mesh>& Z);

/*  template
  void make_box<float>(
                     float const & width
                     , float const & height
                     , float const & depth
                     , T3Mesh & mesh
                     , VertexAttribute<float,T3Mesh> & coords
                     );

  template
  void make_box<double>(
                     double const & width
                     , double const & height
                     , double const & depth
                     , T3Mesh & mesh
                     , VertexAttribute<double,T3Mesh> & coords
                     );*/

} //namespace mesh_array
