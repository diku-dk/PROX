#include <procedural_factory.h>

#include <eigenhelperall.h>

namespace procedural
{

template <typename T>
void make_box_container(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
                        const T& width, const T& height, const T& depth, const T& wall_thickness,
                        MaterialInfo<T> mat_info)
{
    size_t const mid = get_material_id_eigen<T>(mat_info, "Ground");

    GeometryHandleEigen<T> bottom = create_geometry_handle_box<T>(engine, width, wall_thickness, depth);

    GeometryHandleEigen<T> left
        = create_geometry_handle_box<T>(engine, wall_thickness, height + wall_thickness, depth + 2 * wall_thickness);

    GeometryHandleEigen<T> right
        = create_geometry_handle_box<T>(engine, wall_thickness, height + wall_thickness, depth + 2 * wall_thickness);

    GeometryHandleEigen<T> front
        = create_geometry_handle_box<T>(engine, width, height + wall_thickness, wall_thickness);

    GeometryHandleEigen<T> back = create_geometry_handle_box<T>(engine, width, height + wall_thickness, wall_thickness);

    //--- Bottom ------------------------------------------
    {
        EigenVector3<T> const Pw
            = rotate(orientation, EigenVector3<T>(0, -height * 0.5f - wall_thickness * 0.5f, 0)) + position;

        EigenQuaternion<T> const Qw = orientation;

        create_rigid_body<T>(engine, Pw, Qw, bottom, mid, 1, true, "Visualizer/ground");
    }
    //--- Left ------------------------------------------
    {
        EigenVector3<T> const Pw
            = rotate(orientation, EigenVector3<T>(-width * 0.5f - wall_thickness * 0.5f, -wall_thickness * 0.5f, 0))
            + position;

        EigenQuaternion<T> const Qw = orientation;

        create_rigid_body<T>(engine, Pw, Qw, left, mid, 1, true, "Visualizer/ground");
    }
    //--- Right ------------------------------------------
    {
        EigenVector3<T> const Pw
            = rotate(orientation, EigenVector3<T>(width * 0.5f + wall_thickness * 0.5f, -wall_thickness * 0.5f, 0))
            + position;

        EigenQuaternion<T> const Qw = orientation;

        create_rigid_body<T>(engine, Pw, Qw, right, mid, 1, true, "Visualizer/ground");
    }
    //--- Front ------------------------------------------
    {
        EigenVector3<T> const Pw
            = rotate(orientation, EigenVector3<T>(0, -wall_thickness * 0.5f, depth * 0.5f + 0.5f * wall_thickness))
            + position;

        EigenQuaternion<T> const Qw = orientation;

        create_rigid_body<T>(engine, Pw, Qw, front, mid, 1, true, "Visualizer/ground");
    }
    //--- Back ------------------------------------------
    {
        EigenVector3<T> const Pw
            = rotate(orientation, EigenVector3<T>(0, -wall_thickness * 0.5f, -depth * 0.5f - 0.5f * wall_thickness))
            + position;

        EigenQuaternion<T> const Qw = orientation;

        create_rigid_body<T>(engine, Pw, Qw, back, mid, 1, true, "Visualizer/ground");
    }
}

template void make_box_container<float>(content::API* engine, const EigenVector3<float>& position,
                                        const EigenQuaternion<float>& orientation, const float& width,
                                        const float& height, const float& depth, const float& wall_thickness,
                                        MaterialInfo<float> mat_info);

} //namespace procedural
