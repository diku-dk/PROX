#include <factory/mesh_array_make_capsule.h>
#include <factory/mesh_array_profile_sweep.h>

#include <tiny_math_types.h>

#include <vector>
#include <numbers>
namespace mesh_array
{

template<typename T>
void make_capsule(
    const T& radius
    , const T& height
    , size_t const & slices
    , size_t const & segments
    , T3Mesh & mesh
    , VertexAttribute<T,T3Mesh> & X
    , VertexAttribute<T,T3Mesh> & Y
    , VertexAttribute<T,T3Mesh> & Z
    )
  {


        std::vector<EigenVector3<T>> profile;

        profile.resize(segments);

        T const dtheta = std::numbers::pi_v<T> / (segments-1);

        for(size_t i=0;i < segments; ++i )
        {
            T const theta = dtheta*i;

            EigenQuaternion<T> const R = Rotateu( theta ,EigenVector3<T>{0,0,1});
            T const dh = (i < (segments/2)) ? -height/2.0f : height/2.0f;

            profile[ i ] = rotate( R, EigenVector3<T>(  0, -radius, 0 ) ) + EigenVector3<T>(0,dh,0);
        }

        profile_sweep<T>( profile, slices, mesh, X, Y, Z );
    }

    template void make_capsule<float>(float const& radius, float const& height, size_t const& slices,
                                    size_t const& segments, T3Mesh& mesh, VertexAttribute<float, T3Mesh>& X,
                                    VertexAttribute<float, T3Mesh>& Y,
                                    VertexAttribute<float, T3Mesh>& Z);

    template void make_capsule<double>(double const& radius, double const& height, size_t const& slices,
                                    size_t const& segments, T3Mesh& mesh, VertexAttribute<double, T3Mesh>& X,
                                    VertexAttribute<double, T3Mesh>& Y,
                                    VertexAttribute<double, T3Mesh>& Z);

} //namespace mesh_array
