#ifndef MESH_ARRAY_COMPUTE_VERTEX_NORMALS_H
#define MESH_ARRAY_COMPUTE_VERTEX_NORMALS_H

#include <mesh_array_t3mesh.h>
#include <mesh_array_vertex_attribute.h>

#include <cmath>
#include <vector>


namespace mesh_array
{

  template <typename T>
  inline void compute_vertex_normals(
                                     T3Mesh const & mesh
                                     , mesh_array::VertexAttribute<T, T3Mesh> const & X
                                     , mesh_array::VertexAttribute<T, T3Mesh> const & Y
                                     , mesh_array::VertexAttribute<T, T3Mesh> const & Z
                                     , std::vector<EigenVector3<T>> & normals
                                     )
  {
    using std::min;
    using std::max;
    using std::acos;


    unsigned int const cntV = mesh.vertex_size();
    unsigned int const cntT = mesh.triangle_size();

    normals.resize( cntV );

    for(unsigned int v = 0u; v < cntV; ++v)
    {
        EigenVector3<T> vecZero = {T(0),T(0),T(0)};
        normals[v] = vecZero;
    }

    for(unsigned int t = 0u; t < cntT; ++t)
    {
      mesh_array::Triangle const & tri = mesh.triangle( t );

      mesh_array::Vertex const & i = mesh.vertex( tri.i() );
      mesh_array::Vertex const & j = mesh.vertex( tri.j() );
      mesh_array::Vertex const & k = mesh.vertex( tri.k() );

      const EigenVector3<T> pi = EigenVector3<T>( X(i), Y(i), Z(i) );
      const EigenVector3<T> pj = EigenVector3<T>( X(j), Y(j), Z(j) );
      const EigenVector3<T> pk = EigenVector3<T>( X(k), Y(k), Z(k) );

      const EigenVector3<T> a_i      = (pj - pi).normalized();
      const EigenVector3<T> b_i      = (pk - pi).normalized();
      T const dot_i    = min<T>( 1, max<T>( -1, dot(a_i,b_i) ) );
      T const alpha_i  = acos( dot_i );

      const EigenVector3<T> a_j      = (pk - pj).normalized();
      const EigenVector3<T> b_j      = (pi - pj).normalized();
      T const dot_j    = min<T>( 1, max<T>( -1, dot(a_j,b_j) ) );
      T const alpha_j  = acos( dot_j );

      const EigenVector3<T> a_k      = (pi - pk).normalized();
      const EigenVector3<T> b_k      = (pj - pk).normalized();
      T const dot_k    = min<T>( 1, max<T>( -1, dot(a_k,b_k) ) );
      T const alpha_k  = acos( dot_k );

      const EigenVector3<T> a = pj - pi;
      const EigenVector3<T> b = pk - pi;
      const EigenVector3<T> n = ( (a.cross(b)).normalized() );

      normals[ tri.i() ] += n*alpha_i;
      normals[ tri.j() ] += n*alpha_j;
      normals[ tri.k() ] += n*alpha_k;
    }
    for(unsigned int v = 0u; v < cntV; ++v)
    {
        normals[v] = ( normals[v] ).normalized();
    }
  }

} // namespace mesh_array

// MESH_ARRAY_COMPUTE_VERTEX_NORMALS_H
#endif
