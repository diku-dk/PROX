#include <factory/mesh_array_make_sphere.h>
#include <factory/mesh_array_profile_sweep.h>

#include <tiny_math_types.h>
#include <numbers>

#include <cmath>
#include <cassert>

namespace mesh_array
{

template <typename T>
void make_sphere(T const& radius, size_t const& slices, size_t const& segments, T3Mesh& mesh,
                 VertexAttribute<T, T3Mesh>& X, VertexAttribute<T, T3Mesh>& Y, VertexAttribute<T, T3Mesh>& Z)
{
/*    typedef typename MT::real_type       T;
    typedef typename MT::value_traits    VT;
    typedef typename MT::vector3_type    V;
    typedef typename MT::quaternion_type Q;*/

    std::vector<EigenVector3<T>> profile;

    profile.resize(segments);

    T const dtheta = std::numbers::pi_v<T> / (segments - 1);

    for (size_t i = 0; i < segments; ++i)
    {
        T const theta = dtheta * i;

        const EigenQuaternion<T> R = Rotateu(theta, EigenVector3<T>(0, 0, 1));
        profile[i] = rotate(R, EigenVector3<T>(0, -radius, 0));
    }

    profile_sweep<T>(profile, slices, mesh, X, Y, Z);

//    typedef typename MT::real_type       T;
//    typedef typename MT::value_traits    VT;
//    typedef typename MT::vector3_type    V;
//    typedef typename MT::quaternion_type Q;
//
//		using std::cos;
//		using std::sin;
//
//		assert(slices>2u   || !"make_sphere(): must have at least 3 slices"  );
//		assert(segments>1u || !"make_sphere(): must have at least 2 segments");
//
//		size_t const no_quads = segments*slices;
//
//    mesh.clear();
//
//		mesh.set_capacity( no_quads*4, no_quads*2 );
//
//    X.bind(mesh);
//    Y.bind(mesh);
//    Z.bind(mesh);
//
//		T const delta_theta = 2*VT::pi()/slices;
//		T const delta_phi   =      VT::pi()/segments;
//
//		size_t vertex_offset = 0u;
//		for(size_t i=0u; i<segments; ++i)
//		{
//			for(size_t j=0u; j<slices; ++j)
//			{
//				T const theta = j*delta_theta;
//				T const phi   = i*delta_phi;
//
//        mesh.push_vertex();
//        mesh.push_vertex();
//        mesh.push_vertex();
//        mesh.push_vertex();
//
//				Vertex const vi = mesh.vertex(vertex_offset + 0);
//				Vertex const vj = mesh.vertex(vertex_offset + 1);
//				Vertex const vk = mesh.vertex(vertex_offset + 2);
//				Vertex const vm = mesh.vertex(vertex_offset + 3);
//
//        V const pi = radius * V::make(
//                                      cos(theta)*sin(phi+delta_phi)
//                                      ,sin(theta)*sin(phi+delta_phi)
//                                      , cos(phi+delta_phi)
//                                      );
//        V const pj =radius * V::make(
//                                     cos(theta+delta_theta)*sin(phi+delta_phi)
//                                     ,sin(theta+delta_theta)*sin(phi+delta_phi)
//                                     , cos(phi+delta_phi)
//                                     );
//        V const pk =radius * V::make(
//                                     cos(theta+delta_theta)*sin(phi)
//                                     ,sin(theta+delta_theta)*sin(phi)
//                                     , cos(phi)
//                                     );
//        V const pm =radius * V::make(
//                                     cos(theta)*sin(phi)
//                                     ,sin(theta)*sin(phi)
//                                     , cos(phi)
//                                     );
//
//        X(vi) = pi(0);        Y(vi) = pi(1);        Z(vi) = pi(2);
//        X(vj) = pj(0);        Y(vj) = pj(1);        Z(vj) = pj(2);
//        X(vk) = pk(0);        Y(vk) = pk(1);        Z(vk) = pk(2);
//        X(vm) = pm(0);        Y(vm) = pm(1);        Z(vm) = pm(2);
//
//        mesh.push_triangle( vi, vj, vk );
//        mesh.push_triangle( vi, vk, vm );
//
//				vertex_offset += 4u;
//			}
//		}
}

template void make_sphere<float>(float const& radius, size_t const& slices, size_t const& segments, T3Mesh& mesh,
                                 VertexAttribute<float, T3Mesh>& X, VertexAttribute<float, T3Mesh>& Y,
                                 VertexAttribute<float, T3Mesh>& Z);

template void make_sphere<double>(double const& radius, size_t const& slices, size_t const& segments, T3Mesh& mesh,
                                  VertexAttribute<double, T3Mesh>& X, VertexAttribute<double, T3Mesh>& Y,
                                  VertexAttribute<double, T3Mesh>& Z);

} //namespace mesh_array
