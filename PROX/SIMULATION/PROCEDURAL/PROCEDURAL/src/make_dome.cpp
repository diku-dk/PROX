#include <procedural_factory.h>

#include <eigenhelperall.h>

#include <cmath>
#include <cassert>

namespace procedural
{

template <typename T>
void make_dome(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
               const T& r_outer, const T& r_inner, size_t const& slices, size_t const& segments,
               MaterialInfo<T> mat_info)
{
    using std::cos;
    using std::sin;

    assert(r_inner > 0 || !"invalid inner radius value");
    assert(r_outer > r_inner || !"invalid outer radius value");
    assert(slices > 2u || !"at last three slices  must be used");
    assert(segments > 0u || !"at last one segment must be used");

    size_t const mid = get_material_id_eigen<T>(mat_info, "Stone");

        // The dome is created using spherical coordinates as defined here
        //
        // http://mathworld.wolfram.com/SphericalCoordinates.html
        //
        //   x = r cos( theta ) sin( phi )
        //   y = r sin( theta ) sin( phi )
        //   z = r cos( phi )
        //
        // Initially the up-direction of the dome corresponds to the positive z-axis.

    T const offset_theta = 0;
    T const delta_theta = 2 * std::numbers::pi_v<T> / T(slices);
    T const delta_phi = (std::numbers::pi_v<T> * 0.5f) / T(segments + 1u);
    T const stone_density = get_material_density_eigen<T>(mat_info, "Stone");
    T phi = delta_phi * T(1.5f);

    for (size_t v = 0u; v < segments; ++v, phi += delta_phi)
    {
        // When we create the stone we are looking down the negative x-axis direction and having the positive y-axis as our up-direction
        T const cos_phi_upper = cos(phi - delta_phi * 0.5f);
        T const sin_phi_upper = sin(phi - delta_phi * 0.5f);
        T const cos_phi_lower = cos(phi + delta_phi * 0.5f);
        T const sin_phi_lower = sin(phi + delta_phi * 0.5f);
        T const cos_theta_left = cos(-delta_theta / 2);
        T const sin_theta_left = sin(-delta_theta / 2);
        T const cos_theta_right = cos(delta_theta / 2);
        T const sin_theta_right = sin(delta_theta / 2);

        EigenVector3<T> vertices[8];

        // The front face in CCW order
        vertices[0]
            = r_outer * EigenVector3<T>(cos_theta_left * sin_phi_lower, sin_theta_left * sin_phi_lower, cos_phi_lower);
        vertices[1] = r_outer
                    * EigenVector3<T>(cos_theta_right * sin_phi_lower, sin_theta_right * sin_phi_lower, cos_phi_lower);
        vertices[2] = r_outer
                    * EigenVector3<T>(cos_theta_right * sin_phi_upper, sin_theta_right * sin_phi_upper, cos_phi_upper);
        vertices[3]
            = r_outer * EigenVector3<T>(cos_theta_left * sin_phi_upper, sin_theta_left * sin_phi_upper, cos_phi_upper);

        // The back face in CCW order
        vertices[4] = r_inner * vertices[0] / r_outer;
        vertices[5] = r_inner * vertices[1] / r_outer;
        vertices[6] = r_inner * vertices[2] / r_outer;
        vertices[7] = r_inner * vertices[3] / r_outer;

        //--- Compute geometric center of the stone
        const EigenVector3<T> geometric_center = (vertices[0] + vertices[1] + vertices[2] + vertices[3] + vertices[4]
                                                  + vertices[5] + vertices[6] + vertices[7])
                                               / 8;

        //--- Change geometry so geometric center is the model reference.

        for (size_t i = 0; i < 8; ++i) vertices[i] = vertices[i] - geometric_center;

        GeometryHandleEigen<T> stone_handle = create_geometry_handle_cuboid<T>(engine, vertices);

        T theta = offset_theta;
        for (size_t h = 0u; h < slices; ++h, theta += delta_theta)
        {
            const EigenVector3<T> T_b2m = (stone_handle.Tb2m());
            const EigenQuaternion<T> Q_b2m = (stone_handle.Qb2m());

            const EigenQuaternion<T> Q_m2l = Rotateu(theta, EigenVector3<T>(0, 0, 1));
            const EigenVector3<T> T_m2l = (rotate<T>((Q_m2l), geometric_center)); //V::make( x, y, z );

            const EigenVector3<T> T_l2w = (position);
            const EigenQuaternion<T> Q_l2w = (orientation);

            EigenVector3<T> T_b2w;
            EigenQuaternion<T> Q_b2w;

            compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w, Q_l2w, T_b2w, Q_b2w);

            create_rigid_body<T>(engine, T_b2w, Q_b2w, stone_handle, mid, stone_density);
        }
    }
}

template void make_dome<float>(content::API* engine, const EigenVector3<float>& position,
                               const EigenQuaternion<float>& orientation, const float& r_outer, const float& r_inner,
                               size_t const& slices, size_t const& segments, MaterialInfo<float> mat_info);

} //namespace procedural
