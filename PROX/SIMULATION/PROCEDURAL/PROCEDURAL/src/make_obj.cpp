#include <mesh_array.h>

#include <procedural_factory.h>

#include <tiny_math_types.h>

namespace procedural
{

template <typename MT>
void make_obj(content::API* engine, std::string const& name, typename MT::real_type const& scale,
              typename MT::vector3_type const& position, typename MT::quaternion_type const& orientation,
              MaterialInfo<typename MT::real_type> mat_info, bool const fixed = false, bool const blind_copy = false,
              std::string const& material = "Stone",
              mesh_array::TetGenSettings tetset = mesh_array::tetgen_default_settings())
{
    typedef typename MT::real_type T;
    typedef typename MT::vector3_type V;
    typedef typename MT::quaternion_type Q;
    typedef typename MT::value_traits VT;

    T const stone_density = get_material_density_eigen<T>(mat_info, material);
    size_t const mid = get_material_id_eigen<T>(mat_info, material);

    GeometryHandleEigen<T> obj_handle = (create_geometry_handle_obj<T>(engine, name, scale, blind_copy, tetset));

    const EigenVector3<T> T_b2m = obj_handle.Tb2m();
    const EigenQuaternion<T> Q_b2m = obj_handle.Qb2m();

    const EigenVector3<T> T_m2l = EigenVector3<T>(0, 0.5f, 0);
    const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

    const EigenVector3<T> T_l2w = blind_copy ? EigenVector3<T>(0, 0, 0) : toEigen(position);
    const EigenQuaternion<T> Q_l2w = blind_copy ? EigenQuaternion<T>::Identity() : toEigen(orientation);

    EigenVector3<T> T_b2w;
    EigenQuaternion<T> Q_b2w;

    compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w, Q_l2w, T_b2w, Q_b2w);

    create_rigid_body<T>(engine, T_b2w, Q_b2w, obj_handle, mid, stone_density, fixed);
}

using MTf = tiny::MathTypes<float>;

template void make_obj<MTf>(content::API* engine, std::string const& name, MTf::real_type const& scale,
                            MTf::vector3_type const& position, MTf::quaternion_type const& orientation,
                            MaterialInfo<MTf::real_type> mat_info, bool const fixed, bool const blind_copy,
                            std::string const& material, mesh_array::TetGenSettings tetset);

} //namespace procedural
