#include <procedural_factory.h>

#include <eigenhelperall.h>

#include <cmath> // needed for std::floor

namespace procedural
{

template <typename T>
void make_sphere_packing(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
                         const T& min_radius, const T& max_radius, const T& width, const T& height, const T& depth,
                         size_t const& number_of_spheres, MaterialInfo<T> mat_info)
{
    using std::ceil;
    using std::floor;

    T const stone_density = get_material_density_eigen<T>(mat_info, "Stone");
    size_t const mid = get_material_id_eigen<T>(mat_info, "Stone");

    std::vector< GeometryHandleEigen<T> > sphere_handle;
    sphere_handle.resize(10u);
    for (unsigned int c = 0; c < 10u; ++c)
    {
        T const radius = (c + 1u) * (max_radius - min_radius) / (10.0f) + min_radius;
        sphere_handle[c] = create_geometry_handle_sphere<T>(engine, radius);
    }

    unsigned int const I = std::floor(width / (max_radius * 2));
    unsigned int const K = std::floor(depth / (max_radius * 2));
    unsigned int const J = std::ceil(1.0f * number_of_spheres / (I * K));

    unsigned int sphere_count = 0u;

    for (unsigned int j = 0u; j < J; ++j)
        for (unsigned int k = 0u; k < K; ++k)
            for (unsigned int i = 0u; i < I; ++i)
            {
                if (sphere_count > number_of_spheres) return;

                T const x = -width * 0.5f + max_radius + i * max_radius * 2;
                T const z = -depth * 0.5f + max_radius + k * max_radius * 2;
                T const y = -height * 0.5f + max_radius + j * max_radius * 2;

                int choice = sphere_count % 10;

                const EigenVector3<T> T_b2m = sphere_handle[choice].Tb2m();
                const EigenQuaternion<T> Q_b2m = sphere_handle[choice].Qb2m();

                const EigenVector3<T> T_m2l = EigenVector3<T>(x, y, z);
                const EigenQuaternion<T> Q_m2l = Rotateu(-std::numbers::pi_v<T> * 0.5f, EigenVector3<T>(1, 0, 0));

                const EigenVector3<T> T_l2w = (position);
                const EigenQuaternion<T> Q_l2w = (orientation);

                EigenVector3<T> T_b2w;
                EigenQuaternion<T> Q_b2w;

                compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w, Q_l2w, T_b2w, Q_b2w);

                create_rigid_body<T>(engine, T_b2w, Q_b2w, sphere_handle[choice], mid, stone_density);
                ++sphere_count;
            }
}

template void make_sphere_packing<float>(content::API* engine, const EigenVector3<float>& position,
                                         const EigenQuaternion<float>& orientation, const float& min_radius,
                                         const float& max_radius, const float& width, const float& height,
                                         const float& depth, size_t const& number_of_spheres,
                                         MaterialInfo<float> mat_info);

} //namespace procedural
