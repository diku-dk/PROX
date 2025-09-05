#include <procedural_factory.h>

#include <tiny_math_types.h>

namespace procedural
{

template <typename T>
void make_sphere_layer(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
                       const T& sphere_radius, size_t const& spheres_width, size_t const& spheres_length,
                       MaterialInfo<T> mat_info)
{

    T const stone_density = get_material_density_eigen<T>(mat_info, "Stone");
    size_t const mid = get_material_id_eigen<T>(mat_info, "Stone");

    GeometryHandleEigen<T> sphere_handle = create_geometry_handle_sphere<T>(engine, sphere_radius);

    T x = 1;
    T y = sphere_radius;
    T z = 1;
    T const offset = sphere_radius * spheres_length * 0.5f;

    for (size_t i = 0; i < spheres_width; ++i)
    {
        x = 2 * i * sphere_radius - offset;
        for (size_t j = 0; j < spheres_length; ++j)
        {
            z = sphere_radius + 2 * j * sphere_radius;
            const EigenVector3<T> Tb = sphere_handle.Tb2m();
            const EigenQuaternion<T> Qb = sphere_handle.Qb2m();

            const EigenVector3<T> Tm = EigenVector3<T>(x, y, z);
            const EigenQuaternion<T> Qm = Rotateu(-std::numbers::pi_v<T> * 0.5f, EigenVector3<T>(1, 0, 0));

            const EigenVector3<T> Tw = rotate(Qm, Tb) + Tm;
            const EigenQuaternion<T> Qw = Qm * Qb;

            const EigenVector3<T> Tu = rotate((orientation), Tw) + (position);
            const EigenQuaternion<T> Qu = (orientation)*Qw;

            create_rigid_body<T>(engine, Tu, Qu, sphere_handle, mid, stone_density);
        }
    }
}

template void make_sphere_layer<float>(content::API* engine, const EigenVector3<float>& position,
                                       const EigenQuaternion<float>& orientation, const float& sphere_radius,
                                       size_t const& spheres_width, size_t const& spheres_length,
                                       MaterialInfo<float> mat_info);

} //namespace procedural
