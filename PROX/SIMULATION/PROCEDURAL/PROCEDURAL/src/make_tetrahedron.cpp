#include <procedural_factory.h>

#include <tiny_math_types.h>

namespace procedural
{

template <typename T>
void make_tetrahedron(content::API* engine, const EigenVector3<T>& one, const EigenVector3<T>& two,
                      const EigenVector3<T>& three, const EigenVector3<T>& four, const EigenVector3<T>& position,
                      const EigenQuaternion<T>& orientation, MaterialInfo<T> mat_info, bool const fixed = false)
{;

    T const stone_density = get_material_density_eigen<T>(mat_info, "Stone");
    size_t const mid = get_material_id_eigen<T>(mat_info, "Stone");

    GeometryHandleEigen<T> tet_handle = create_geometry_handle_tetrahedron<T>(engine, (one), (two), (three), (four));

    const EigenVector3<T> T_b2m = (tet_handle.Tb2m());
    const EigenQuaternion<T> Q_b2m = (tet_handle.Qb2m());

    const EigenVector3<T> T_m2l = EigenVector3<T>(0, 0, 0);
    const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

    const EigenVector3<T> T_l2w = (position);
    const EigenQuaternion<T> Q_l2w = (orientation);

    EigenVector3<T> T_b2w;
    EigenQuaternion<T> Q_b2w;

    compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w, Q_l2w, T_b2w, Q_b2w);

    create_rigid_body<T>(engine, T_b2w, Q_b2w, tet_handle, mid, stone_density, fixed);
}

template void make_tetrahedron<float>(content::API* engine, const EigenVector3<float>& one,
                                      const EigenVector3<float>& two, const EigenVector3<float>& three,
                                      const EigenVector3<float>& four, const EigenVector3<float>& position,
                                      const EigenQuaternion<float>& orientation, MaterialInfo<float> mat_info,
                                      bool const fixed);

template <typename T>
void make_tetrahedron(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
                      MaterialInfo<T> mat_info, bool const fixed = false)
{
    EigenVector3<T> one = EigenVector3<T>(0.0, 0.0, 0.0);
    EigenVector3<T> two = EigenVector3<T>(1.0, 0.0, 0.0);
    EigenVector3<T> three = EigenVector3<T>(0.0, 0.0, 1.0);
    EigenVector3<T> four = EigenVector3<T>(0.0, 1.0, 0.0);

    make_tetrahedron<T>(engine, one, two, three, four, position, orientation, mat_info);
}

template void make_tetrahedron<float>(content::API* engine, const EigenVector3<float>& position,
                                      const EigenQuaternion<float>& orientation, MaterialInfo<float> mat_info,
                                      bool const fixed);

} //namespace procedural
