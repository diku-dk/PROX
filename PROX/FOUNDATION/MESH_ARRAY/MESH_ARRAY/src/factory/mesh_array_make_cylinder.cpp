#include <factory/mesh_array_make_cylinder.h>
#include <factory/mesh_array_profile_sweep.h>

#include <tiny_math_types.h>

#include <vector>

namespace mesh_array
{

  template<typename T>
  void make_cylinder(
                    T const & radius
                    , T const & height
                    , size_t const & slices
                    , T3Mesh & mesh
                    , VertexAttribute<T,T3Mesh> & X
                    , VertexAttribute<T,T3Mesh> & Y
                    , VertexAttribute<T,T3Mesh> & Z
                    )
  {

      std::vector<EigenVector3<T>> profile;

        profile.resize(4u);

        profile[ 0 ] = EigenVector3<T>( 0, -height/2, 0 );
        profile[ 1 ] = EigenVector3<T>(     radius, -height/2, 0 );
        profile[ 2 ] = EigenVector3<T>(     radius,  height/2, 0 );
        profile[ 3 ] = EigenVector3<T>( 0,  height/2, 0 );

        profile_sweep<T>( profile, slices, mesh, X, Y, Z  );
    }


    template void make_cylinder<float>(float const& radius, float const& height, size_t const& slices,
                                     T3Mesh& mesh, VertexAttribute<float, T3Mesh>& X,
                                     VertexAttribute<float, T3Mesh>& Y,
                                     VertexAttribute<float, T3Mesh>& Z);

    template void make_cylinder<double>(double const& radius, double const& height, size_t const& slices,
                                     T3Mesh& mesh, VertexAttribute<double, T3Mesh>& X,
                                     VertexAttribute<double, T3Mesh>& Y,
                                     VertexAttribute<double, T3Mesh>& Z);

} //namespace mesh_array
