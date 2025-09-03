#include <procedural_factory.h>

#include <tiny_math_types.h>

namespace procedural
{

template <typename MT>
void make_box_container(content::API* engine, typename MT::vector3_type const& position,
                        typename MT::quaternion_type const& orientation, typename MT::real_type const& width,
                        typename MT::real_type const& height, typename MT::real_type const& depth,
                        typename MT::real_type const& wall_thickness, MaterialInfo<typename MT::real_type> mat_info)
{
    typedef typename MT::vector3_type V;
    typedef typename MT::quaternion_type Q;
    typedef typename MT::value_traits VT;

    size_t const mid = get_material_id<MT>(mat_info, "Ground");

    GeometryHandle<MT> bottom = create_geometry_handle_box<MT>(engine, width, wall_thickness, depth);

    GeometryHandle<MT> left
        = create_geometry_handle_box<MT>(engine, wall_thickness, height + wall_thickness, depth + 2 * wall_thickness);

    GeometryHandle<MT> right
        = create_geometry_handle_box<MT>(engine, wall_thickness, height + wall_thickness, depth + 2 * wall_thickness);

    GeometryHandle<MT> front = create_geometry_handle_box<MT>(engine, width, height + wall_thickness, wall_thickness);

    GeometryHandle<MT> back = create_geometry_handle_box<MT>(engine, width, height + wall_thickness, wall_thickness);

    //--- Bottom ------------------------------------------
    {
        V const Pw = rotate(orientation, V::make(0, -height * 0.5f - wall_thickness * 0.5f, 0)) + position;

        Q const Qw = orientation;

        create_rigid_body<MT>(engine, Pw, Qw, bottom, mid, 1, true, "Visualizer/ground");
    }
    //--- Left ------------------------------------------
    {
        V const Pw
            = rotate(orientation, V::make(-width * 0.5f - wall_thickness * 0.5f, -wall_thickness * 0.5f, 0)) + position;

        Q const Qw = orientation;

        create_rigid_body<MT>(engine, Pw, Qw, left, mid, 1, true, "Visualizer/ground");
    }
    //--- Right ------------------------------------------
    {
        V const Pw
            = rotate(orientation, V::make(width * 0.5f + wall_thickness * 0.5f, -wall_thickness * 0.5f, 0)) + position;

        Q const Qw = orientation;

        create_rigid_body<MT>(engine, Pw, Qw, right, mid, 1, true, "Visualizer/ground");
    }
    //--- Front ------------------------------------------
    {
        V const Pw
            = rotate(orientation, V::make(0, -wall_thickness * 0.5f, depth * 0.5f + 0.5f * wall_thickness)) + position;

        Q const Qw = orientation;

        create_rigid_body<MT>(engine, Pw, Qw, front, mid, 1, true, "Visualizer/ground");
    }
    //--- Back ------------------------------------------
    {
        V const Pw
            = rotate(orientation, V::make(0, -wall_thickness * 0.5f, -depth * 0.5f - 0.5f * wall_thickness)) + position;

        Q const Qw = orientation;

        create_rigid_body<MT>(engine, Pw, Qw, back, mid, 1, true, "Visualizer/ground");
    }
}

using MTf = tiny::MathTypes<float>;

template void make_box_container<MTf>(content::API* engine, MTf::vector3_type const& position,
                                      MTf::quaternion_type const& orientation, MTf::real_type const& width,
                                      MTf::real_type const& height, MTf::real_type const& depth,
                                      MTf::real_type const& wall_thickness, MaterialInfo<MTf::real_type> mat_info);

} //namespace procedural
