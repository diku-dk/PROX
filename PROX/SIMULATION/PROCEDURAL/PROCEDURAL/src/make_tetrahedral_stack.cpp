#include <procedural_factory.h>

#include <tiny_math_types.h>

namespace procedural
{

template <typename T>
void make_tetrahedral_stack(content::API* engine, const EigenVector3<T>& position,
                            const EigenQuaternion<T>& orientation, const T& cube_width, const T& cube_height,
                            const T& cube_depth, size_t const& sub_divisions, MaterialInfo<T> mat_info)
{

    T const stone_density = get_material_density_eigen<T>(mat_info, "Stone");
    size_t const mid = get_material_id_eigen<T>(mat_info, "Stone");

        //--- Find the size of each of the boxes
    T const box_width = cube_width / sub_divisions;
    T const box_height = cube_height / sub_divisions;
    T const box_depth = cube_depth / sub_divisions;

        //--- The dimensions of the box is in place now we will make the individual tetras that make up the box
    EigenVector3<T> one = EigenVector3<T>(-box_width / 2.0, -box_height / 2.0, -box_depth / 2.0);
    EigenVector3<T> two = EigenVector3<T>(box_width / 2.0, -box_height / 2.0, -box_depth / 2.0);
    EigenVector3<T> three = EigenVector3<T>(box_width / 2.0, box_height / 2.0, -box_depth / 2.0);
    EigenVector3<T> four = EigenVector3<T>(-box_width / 2.0, box_height / 2.0, -box_depth / 2.0);

    EigenVector3<T> five = EigenVector3<T>(-box_width / 2.0, -box_height / 2.0, box_depth / 2.0);
    EigenVector3<T> six = EigenVector3<T>(box_width / 2.0, -box_height / 2.0, box_depth / 2.0);
    EigenVector3<T> seven = EigenVector3<T>(box_width / 2.0, box_height / 2.0, box_depth / 2.0);
    EigenVector3<T> eight = EigenVector3<T>(-box_width / 2.0, box_height / 2.0, box_depth / 2.0);

    GeometryHandleEigen<T> tetra_handle1 = create_geometry_handle_tetrahedron<T>(engine, one, two, three, six);
    GeometryHandleEigen<T> tetra_handle2 = create_geometry_handle_tetrahedron<T>(engine, one, three, four, eight);
    GeometryHandleEigen<T> tetra_handle3 = create_geometry_handle_tetrahedron<T>(engine, five, eight, six, one);
    GeometryHandleEigen<T> tetra_handle4 = create_geometry_handle_tetrahedron<T>(engine, six, eight, seven, three);
    GeometryHandleEigen<T> tetra_handle5 = create_geometry_handle_tetrahedron<T>(engine, one, three, eight, six);

    T x = 1;
    T y = 1;
    T z = 1;

    for (size_t i = 0; i < sub_divisions; ++i)
    {
        x = i * box_width - (cube_width / 2) + (box_width / 2);
        for (size_t j = 0; j < sub_divisions; ++j)
        {
            y = (box_height / 2) + j * box_height;
            for (size_t k = 0; k < sub_divisions; ++k)
            {
                z = k * box_depth;

                const EigenVector3<T> Tb = (tetra_handle1.Tb2m());
                const EigenQuaternion<T> Qb = (tetra_handle1.Qb2m());
                const EigenVector3<T> Tm = EigenVector3<T>(x, y, z);
                const EigenQuaternion<T> Qm = EigenQuaternion<T>::Identity();
                const EigenVector3<T> Tu = rotate(Qm, Tb) + Tm;
                const EigenQuaternion<T> Qu = Qm * Qb;
                create_rigid_body<T>(engine, Tu, Qu, tetra_handle1, mid, stone_density);

                const EigenVector3<T> Tb2 = (tetra_handle2.Tb2m());
                const EigenQuaternion<T> Qb2 = (tetra_handle2.Qb2m());
                const EigenVector3<T> Tu2 = rotate(Qm, Tb2) + Tm;
                const EigenQuaternion<T> Qu2 = Qm * Qb2;
                create_rigid_body<T>(engine, Tu2, Qu2, tetra_handle2, mid, stone_density);

                const EigenVector3<T> Tb3 = (tetra_handle3.Tb2m());
                const EigenQuaternion<T> Qb3 = (tetra_handle3.Qb2m());
                const EigenVector3<T> Tu3 = rotate(Qm, Tb3) + Tm;
                const EigenQuaternion<T> Qu3 = Qm * Qb3;
                create_rigid_body<T>(engine, Tu3, Qu3, tetra_handle3, mid, stone_density);

                const EigenVector3<T> Tb4 = (tetra_handle4.Tb2m());
                const EigenQuaternion<T> Qb4 = (tetra_handle4.Qb2m());
                const EigenVector3<T> Tu4 = rotate(Qm, Tb4) + Tm;
                const EigenQuaternion<T> Qu4 = Qm * Qb4;
                create_rigid_body<T>(engine, Tu4, Qu4, tetra_handle4, mid, stone_density);

                const EigenVector3<T> Tb5 = (tetra_handle5.Tb2m());
                const EigenQuaternion<T> Qb5 = (tetra_handle5.Qb2m());
                const EigenVector3<T> Tu5 = rotate(Qm, Tb5) + Tm;
                const EigenQuaternion<T> Qu5 = Qm * Qb5;
                create_rigid_body<T>(engine, Tu5, Qu5, tetra_handle5, mid, stone_density);
            }
        }
    }
}

template void make_tetrahedral_stack<float>(content::API* engine, const EigenVector3<float>& position,
                                            const EigenQuaternion<float>& orientation, const float& cube_width,
                                            const float& cube_height, const float& cube_depth,
                                            size_t const& sub_divisions, MaterialInfo<float> mat_info);

} //namespace procedural
