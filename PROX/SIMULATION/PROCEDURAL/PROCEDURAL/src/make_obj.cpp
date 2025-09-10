#include <mesh_array.h>

#include <procedural_factory.h>

#include <eigenhelperall.h>

namespace procedural
{

template <typename T>
void make_obj(content::API* engine, std::string const& name, const T& scale, const EigenVector3<T>& position,
              const EigenQuaternion<T>& orientation, MaterialInfo<T> mat_info, bool const fixed = false,
              bool const blind_copy = false, std::string const& material = "Stone",
              mesh_array::TetGenSettings tetset = mesh_array::tetgen_default_settings())
{
    T const stone_density = get_material_density_eigen<T>(mat_info, material);
    size_t const mid = get_material_id_eigen<T>(mat_info, material);

    GeometryHandleEigen<T> obj_handle = (create_geometry_handle_obj<T>(engine, name, scale, blind_copy, tetset));

    const EigenVector3<T> T_b2m = obj_handle.Tb2m();
    const EigenQuaternion<T> Q_b2m = obj_handle.Qb2m();

    const EigenVector3<T> T_m2l = EigenVector3<T>(0, 0.5f, 0);
    const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

    const EigenVector3<T> T_l2w = blind_copy ? EigenVector3<T>(0, 0, 0) : (position);
    const EigenQuaternion<T> Q_l2w = blind_copy ? EigenQuaternion<T>::Identity() : (orientation);

    EigenVector3<T> T_b2w;
    EigenQuaternion<T> Q_b2w;

    compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w, Q_l2w, T_b2w, Q_b2w);

    create_rigid_body<T>(engine, T_b2w, Q_b2w, obj_handle, mid, stone_density, fixed);
}

template void make_obj<float>(content::API* engine, std::string const& name, const float& scale,
                              const EigenVector3<float>& position, const EigenQuaternion<float>& orientation,
                              MaterialInfo<float> mat_info, bool const fixed, bool const blind_copy,
                              std::string const& material, mesh_array::TetGenSettings tetset);

} //namespace procedural
