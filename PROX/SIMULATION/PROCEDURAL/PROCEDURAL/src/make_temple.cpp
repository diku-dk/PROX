#include <procedural.h>
#include <procedural_factory.h>

#include <tiny_math_types.h>

namespace procedural
{

template <typename MT>
void make_temple(content::API* engine, typename MT::vector3_type const& position,
                 typename MT::quaternion_type const& orientation, typename MT::real_type const& temple_height,
                 typename MT::real_type const& pillar_width, size_t const& num_pillars_x, size_t const& num_pillars_z,
                 MaterialInfo<typename MT::real_type> mat_info)
{
    typedef typename MT::real_type T;
    typedef typename MT::vector3_type V;
    typedef typename MT::quaternion_type Q;
    typedef typename MT::value_traits VT;

    using std::atan;
    using std::floor;
    using std::tan;

        /// temple height ratios, must sum to one
    T const bottom_ratio = VT::numeric_cast(0.03f);//there are two of these
    T const pillar_ratio = VT::numeric_cast(0.61f);
    T const gable_ratio = VT::numeric_cast(0.3f);
    T const beam_ratio = VT::numeric_cast(0.03f);

    T const temple_width = (2 * num_pillars_x - 1) * pillar_width;
    T const temple_depth = (2 * num_pillars_z - 1) * pillar_width;

    T const beam_spacing_ratio = VT::numeric_cast(0.7f);
    T const spacing = (temple_depth - 4 * pillar_width) / floor(beam_spacing_ratio * num_pillars_z);
    T const pillar_height = temple_height * pillar_ratio;
    T const plane_height = temple_height * bottom_ratio;

    T const beam_height = temple_height * beam_ratio;
    T const beam_length_f = (temple_width - VT::numeric_cast(5.0f) * pillar_width) / (num_pillars_x - 3);
    T const beam_length_s = (temple_depth - VT::numeric_cast(5.0f) * pillar_width) / (num_pillars_z - 3);

    T const gable_height = temple_height * gable_ratio;
    T gable_num_brick = 1 * num_pillars_x - 2;
    T const gable_incline = atan(2 * gable_height / temple_width);
    T const gable_brick_w = (temple_width - 3 * pillar_width) / gable_num_brick;
    T const gable_brick_h = tan(gable_incline) * 0.5f * gable_brick_w;
    T gable_num_layers = (gable_height / gable_brick_h);

    size_t const mid = get_material_id<MT>(mat_info, "Stone");
    T const stone_density = get_material_density<MT>(mat_info, "Stone");

        /// y offsets
    T const plane_1_y = plane_height * 0.5f;
    T const plane_2_y = plane_1_y + plane_height;
    T const plane_3_y = plane_2_y + plane_height;
    T const outer_pillar_y = plane_2_y + plane_height * 0.5f;
    T const inner_pillar_y = outer_pillar_y + plane_height;
    T const beam_y = outer_pillar_y + pillar_height + 0.5f * beam_height;
    T const gable_y = beam_y + 0.5f * beam_height + 0.5f * gable_brick_h;

    size_t const pillar_segments = 5u;

    EigenVector3<T> vertices[8];

    V Tm, Tb, Tw, Tu;
    Q Qm, Qb, Qw, Qu;

    { /// bottom planes
        GeometryHandle<MT> plane_1 = create_geometry_handle_box<MT>(
            engine, pillar_width + temple_width, bottom_ratio * temple_height, pillar_width + temple_depth);

        V P = rotate(orientation, V::make(0, plane_1_y, 0)) + position;

        create_rigid_body<MT>(engine, P, orientation, plane_1, mid, stone_density);

        GeometryHandle<MT> plane_2
            = create_geometry_handle_box<MT>(engine, temple_width, bottom_ratio * temple_height, temple_depth);

        P = rotate(orientation, V::make(0, plane_2_y, 0)) + position;

        create_rigid_body<MT>(engine, P, orientation, plane_2, mid, stone_density);

        if ((num_pillars_z >= 4) && (num_pillars_x > 2))
        {
                /// foundation for inner pillars
            GeometryHandle<MT> plane_2 = create_geometry_handle_box<MT>(
                engine, temple_width - 4 * pillar_width, bottom_ratio * temple_height, temple_depth - 4 * pillar_width);

            P = rotate(orientation, V::make(0, plane_3_y, 0)) + position;

            create_rigid_body<MT>(engine, P, orientation, plane_2, mid, stone_density);
        }

    } /// bottom planes

    { /// outer pillars

            /// front and back pillar row
        for (size_t i = 0u; i < num_pillars_x; ++i)
        {
            T const xf = (pillar_width - temple_width) * 0.5f + 2 * i * pillar_width;
            T const zf = (-pillar_width + temple_depth) * 0.5f;
            V P = rotate(orientation, V::make(xf, outer_pillar_y, zf)) + position;
            make_greek_pillar<MT>(engine, P, orientation, pillar_width, pillar_height, pillar_width, pillar_segments,
                                  12u, mat_info);

            T const xb = (pillar_width - temple_width) * 0.5f + 2 * i * pillar_width;
            T const zb = (pillar_width - temple_depth) * 0.5f;
            P = rotate(orientation, V::make(xb, outer_pillar_y, zb)) + position;
            make_greek_pillar<MT>(engine, P, orientation, pillar_width, pillar_height, pillar_width, pillar_segments,
                                  12u, mat_info);
        }

            /// side row pillars
        for (size_t i = 1u; i < (num_pillars_z - 1); ++i)
        {

            T const xe = (pillar_width - temple_width) * 0.5f;
            T const ze = (-pillar_width + temple_depth) * 0.5f - 2 * i * pillar_width;
            V P = rotate(orientation, V::make(xe, outer_pillar_y, ze)) + position;
            make_greek_pillar<MT>(engine, P, orientation, pillar_width, pillar_height, pillar_width, pillar_segments,
                                  12u, mat_info);

            T const xw = (-pillar_width + temple_width) * 0.5f;
            T const zw = (-pillar_width + temple_depth) * 0.5f - 2 * i * pillar_width;
            P = rotate(orientation, V::make(xw, outer_pillar_y, zw)) + position;
            make_greek_pillar<MT>(engine, P, orientation, pillar_width, pillar_height, pillar_width, pillar_segments,
                                  12u, mat_info);
        }
    }/// outer pillars

    { /// inner pillars
        if ((num_pillars_z >= 4) && (num_pillars_x > 2))
        {
                /// end bricks
            vertices[0] = EigenVector3<T>(0, 0, 0);
            vertices[1] = EigenVector3<T>(VT::numeric_cast(1.1f) * pillar_width, 0, 0);
            vertices[2] = EigenVector3<T>(VT::numeric_cast(1.1f) * pillar_width, gable_brick_h, 0);
            vertices[3] = EigenVector3<T>(pillar_width, gable_brick_h, 0);
            vertices[4] = vertices[0] - EigenVector3<T>(0, 0, pillar_width);
            vertices[5] = vertices[1] - EigenVector3<T>(0, 0, pillar_width);
            vertices[6] = vertices[2] - EigenVector3<T>(0, 0, pillar_width);
            vertices[7] = vertices[3] - EigenVector3<T>(0, 0, pillar_width);

            GeometryHandleEigen<T> gable_end_brick = create_geometry_handle_cuboid<T>(engine, vertices);

                /// gable bricks
            GeometryHandle<MT> gable_brick
                = create_geometry_handle_box<MT>(engine, gable_brick_w, gable_brick_h, pillar_width);

                /// inner pillars
            for (size_t i = 0u; i < floor(beam_spacing_ratio * num_pillars_z); ++i)
            {
                T z = 0.5f * temple_depth - 2 * pillar_width - (0.5f + i) * spacing;

                for (size_t j = 0u; j < num_pillars_x - 2; ++j)
                {
                    T x = (pillar_width - temple_width) * 0.5f + (2 + 2 * j) * pillar_width;
                    V P = rotate(orientation, V::make(x, inner_pillar_y, z)) + position;

                    make_greek_pillar<MT>(engine, P, orientation, pillar_width, pillar_height, pillar_width,
                                          pillar_segments, 12u, mat_info);
                }

                T const xstart = -0.5f * temple_width + 0.5f * pillar_width;
                T y = gable_y;
                T bricks = gable_num_brick;

                for (size_t i = 0; i < gable_num_layers; ++i)
                {
            //layer no.

                    T x = xstart + 0.5f * i * gable_brick_w;

                    for (size_t j = 0; j < bricks; ++j)
                    {
              //brick placement

                        x += gable_brick_w;

                        Tm = V::make(x, y, z);

                        create_rigid_body<MT>(engine, Tm, orientation, gable_brick, mid, stone_density);
                    }

                    y += gable_brick_h;
                    --bricks;
                }
                for (size_t i = 0; i < gable_num_layers - 1; ++i)
                {
                        //place end bricks
                    Tm = V::make(-0.5f * temple_width + (VT::numeric_cast(0.4f) + i) * pillar_width,
                                 gable_y + (-0.5f + i) * gable_brick_h, z + 0.5f * pillar_width);
                    Qm = Q::identity();

                    Tb = fromEigen(gable_end_brick.Tb2m());
                    Qb = fromEigen(gable_end_brick.Qb2m());

                    Tw = rotate(Qm, Tb) + Tm;
                    Qw = Qm * Qb;

                    Tu = rotate(orientation, Tw) + position;
                    Qu = orientation * Qw;

                    create_rigid_body<MT>(engine, Tu, Qu, gable_end_brick, mid, stone_density);

                    Tm = V::make(0.5f * temple_width - (VT::numeric_cast(0.4f) + i) * pillar_width,
                                 gable_y + (-0.5f + i) * gable_brick_h, z - 0.5f * pillar_width);
                    Qm = Q::Ru(VT::pi(), V::j());

                    Tb = fromEigen(gable_end_brick.Tb2m());
                    Qb = fromEigen(gable_end_brick.Qb2m());

                    Tw = rotate(Qm, Tb) + Tm;
                    Qw = Qm * Qb;

                    Tu = rotate(orientation, Tw) + position;
                    Qu = orientation * Qw;

                    create_rigid_body<MT>(engine, Tu, Qu, gable_end_brick, mid, stone_density);
                }
            }
        }

    }/// inner pillars

    { /// roof 'beams'
        GeometryHandle<MT> beam_f = create_geometry_handle_box<MT>(engine, beam_length_f, beam_height, pillar_width);
        for (size_t i = 0; i < num_pillars_x - 3; ++i)
        {

            T x = -0.5f * temple_width + VT::numeric_cast(3.5f) * pillar_width + i * beam_length_f;
            Tm = rotate(orientation, V::make(x, beam_y, 0.5f * (temple_depth - pillar_width))) + position;

            create_rigid_body<MT>(engine, Tm, orientation, beam_f, mid, stone_density);

            Tm = rotate(orientation, V::make(x, beam_y, 0.5f * (-temple_depth + pillar_width))) + position;

            create_rigid_body<MT>(engine, Tm, orientation, beam_f, mid, stone_density);
        }

        GeometryHandle<MT> beam_s = create_geometry_handle_box<MT>(engine, pillar_width, beam_height, beam_length_s);
        for (size_t i = 0; i < num_pillars_z - 3; ++i)
        {

            T z = 0.5f * temple_depth - VT::numeric_cast(3.5f) * pillar_width - i * beam_length_s;

            Tm = rotate(orientation, V::make(0.5f * (-temple_width + pillar_width), beam_y, z)) + position;

            create_rigid_body<MT>(engine, Tm, orientation, beam_s, mid, stone_density);

            Tm = rotate(orientation, V::make(0.5f * (temple_width - pillar_width), beam_y, z)) + position;

            create_rigid_body<MT>(engine, Tm, orientation, beam_s, mid, stone_density);
        }

        vertices[0] = EigenVector3<T>(0, 0, 0);
        vertices[1] = EigenVector3<T>(VT::numeric_cast(2.5f) * pillar_width, 0, 0);
        vertices[2] = EigenVector3<T>(VT::numeric_cast(2.5f) * pillar_width, 0, pillar_width);
        vertices[3] = EigenVector3<T>(pillar_width, 0, pillar_width);
        vertices[4] = vertices[0] + EigenVector3<T>(0, beam_height, 0);
        vertices[5] = vertices[1] + EigenVector3<T>(0, beam_height, 0);
        vertices[6] = vertices[2] + EigenVector3<T>(0, beam_height, 0);
        vertices[7] = vertices[3] + EigenVector3<T>(0, beam_height, 0);

        GeometryHandleEigen<T> beam_c = create_geometry_handle_cuboid<T>(engine, vertices);

        T const corner_y = beam_y - 0.5f * beam_height;

        Tm = V::make(-0.5f * temple_width, corner_y, -0.5f * temple_depth);
        Qm = Q::identity();//Q::Ru( VT::pi(), V::i() );

        Tb = fromEigen(beam_c.Tb2m());
        Qb = fromEigen(beam_c.Qb2m());

        Tw = rotate(Qm, Tb) + Tm;
        Qw = Qm * Qb;

        Tu = rotate(orientation, Tw) + position;
        Qu = orientation * Qw;

        create_rigid_body<MT>(engine, Tu, Qu, beam_c, mid, stone_density);

        Tm = V::make(0.5f * temple_width, corner_y, -0.5f * temple_depth);
        Qm = Q::Ru(-VT::pi_half(), V::j()) * Qm;

        Tb = fromEigen(beam_c.Tb2m());
        Qb = fromEigen(beam_c.Qb2m());

        Tw = rotate(Qm, Tb) + Tm;
        Qw = Qm * Qb;

        Tu = rotate(orientation, Tw) + position;
        Qu = orientation * Qw;

        create_rigid_body<MT>(engine, Tu, Qu, beam_c, mid, stone_density);

        Tm = V::make(0.5f * temple_width, corner_y, 0.5f * temple_depth);
        Qm = Q::Ru(-VT::pi_half(), V::j()) * Qm;

        Tb = fromEigen(beam_c.Tb2m());
        Qb = fromEigen(beam_c.Qb2m());

        Tw = rotate(Qm, Tb) + Tm;
        Qw = Qm * Qb;

        Tu = rotate(orientation, Tw) + position;
        Qu = orientation * Qw;

        create_rigid_body<MT>(engine, Tu, Qu, beam_c, mid, stone_density);

        Tm = V::make(-0.5f * temple_width, corner_y, 0.5f * temple_depth);
        Qm = Q::Ru(-VT::pi_half(), V::j()) * Qm;

        Tb = fromEigen(beam_c.Tb2m());
        Qb = fromEigen(beam_c.Qb2m());

        Tw = rotate(Qm, Tb) + Tm;
        Qw = Qm * Qb;

        Tu = rotate(orientation, Tw) + position;
        Qu = orientation * Qw;

        create_rigid_body<MT>(engine, Tu, Qu, beam_c, mid, stone_density);

        Tm = V::make(-0.5f * temple_width, corner_y + beam_height, -0.5f * temple_depth);
        Qm = Q::Ru(VT::pi(), V::i()) * Qm;

        Tb = fromEigen(beam_c.Tb2m());
        Qb = fromEigen(beam_c.Qb2m());

        Tw = rotate(Qm, Tb) + Tm;
        Qw = Qm * Qb;

        Tu = rotate(orientation, Tw) + position;
        Qu = orientation * Qw;

        create_rigid_body<MT>(engine, Tu, Qu, beam_c, mid, stone_density);

        Tm = V::make(0.5f * temple_width, corner_y + beam_height, -0.5f * temple_depth);
        Qm = Q::Ru(-VT::pi_half(), V::j()) * Qm;

        Tb = fromEigen(beam_c.Tb2m());
        Qb = fromEigen(beam_c.Qb2m());

        Tw = rotate(Qm, Tb) + Tm;
        Qw = Qm * Qb;

        Tu = rotate(orientation, Tw) + position;
        Qu = orientation * Qw;

        create_rigid_body<MT>(engine, Tu, Qu, beam_c, mid, stone_density);

        Tm = V::make(0.5f * temple_width, corner_y + beam_height, 0.5f * temple_depth);
        Qm = Q::Ru(-VT::pi_half(), V::j()) * Qm;

        Tb = fromEigen(beam_c.Tb2m());
        Qb = fromEigen(beam_c.Qb2m());

        Tw = rotate(Qm, Tb) + Tm;
        Qw = Qm * Qb;

        Tu = rotate(orientation, Tw) + position;
        Qu = orientation * Qw;

        create_rigid_body<MT>(engine, Tu, Qu, beam_c, mid, stone_density);

        Tm = V::make(-0.5f * temple_width, corner_y + beam_height, 0.5f * temple_depth);
        Qm = Q::Ru(-VT::pi_half(), V::j()) * Qm;

        Tb = fromEigen(beam_c.Tb2m());
        Qb = fromEigen(beam_c.Qb2m());

        Tw = rotate(Qm, Tb) + Tm;
        Qw = Qm * Qb;

        Tu = rotate(orientation, Tw) + position;
        Qu = orientation * Qw;

        create_rigid_body<MT>(engine, Tu, Qu, beam_c, mid, stone_density);
    } /// roof 'beams'

    { /// gables
            /// end bricks
        vertices[0] = EigenVector3<T>(0, 0, 0);
        vertices[1] = EigenVector3<T>((1.1f) * pillar_width, 0, 0);
        vertices[2] = EigenVector3<T>((1.1f) * pillar_width, gable_brick_h, 0);
        vertices[3] = EigenVector3<T>(pillar_width, gable_brick_h, 0);

        vertices[4] = vertices[0] - EigenVector3<T>(0, 0, pillar_width * (0.7f));
        vertices[5] = vertices[1] - EigenVector3<T>(0, 0, pillar_width * (0.7f));
        vertices[6] = vertices[2] - EigenVector3<T>(0, 0, pillar_width * (0.7f));
        vertices[7] = vertices[3] - EigenVector3<T>(0, 0, pillar_width * (0.7f));

        GeometryHandleEigen<T> gable_end_brick = create_geometry_handle_cuboid<T>(engine, vertices);

            /// gable bricks
        GeometryHandle<MT> gable_brick
            = create_geometry_handle_box<MT>(engine, gable_brick_w, gable_brick_h, pillar_width * 0.7f);

        for (size_t i = 0; i < gable_num_layers - 1; ++i)
        {

            for (size_t j = 0; j < 2; ++j)
            {
                T z = 0.5f * temple_depth - VT::numeric_cast(0.15f) * pillar_width - j * (-pillar_width + temple_depth);
                    //place end bricks
                Tm = V::make(-0.5f * temple_width + (VT::numeric_cast(0.4f) + i) * pillar_width,
                             gable_y + (-0.5f + i) * gable_brick_h, z);
                Qm = Q::identity();

                Tb = fromEigen(gable_end_brick.Tb2m());
                Qb = fromEigen(gable_end_brick.Qb2m());

                Tw = rotate(Qm, Tb) + Tm;
                Qw = Qm * Qb;

                Tu = rotate(orientation, Tw) + position;
                Qu = orientation * Qw;

                create_rigid_body<MT>(engine, Tu, Qu, gable_end_brick, mid, stone_density);

                Tm = V::make(0.5f * temple_width - (VT::numeric_cast(0.4f) + i) * pillar_width,
                             gable_y + (-0.5f + i) * gable_brick_h, z - VT::numeric_cast(0.7f) * pillar_width);
                Qm = Q::Ru(VT::pi(), V::j());

                Tb = fromEigen(gable_end_brick.Tb2m());
                Qb = fromEigen(gable_end_brick.Qb2m());

                Tw = rotate(Qm, Tb) + Tm;
                Qw = Qm * Qb;

                Tu = rotate(orientation, Tw) + position;
                Qu = orientation * Qw;

                create_rigid_body<MT>(engine, Tu, Qu, gable_end_brick, mid, stone_density);
            }
        }
        for (size_t k = 0; k < 2; ++k)
        {
            T const xstart = -0.5f * temple_width + 0.5f * pillar_width;
            T y = gable_y;
            T z = 0.5f * (temple_depth - pillar_width) - k * (-pillar_width + temple_depth);
            T x = 0;
            T bricks = gable_num_brick;

            for (size_t i = 0; i < gable_num_layers; ++i)
            {//layer no.

                x = xstart + 0.5f * i * gable_brick_w;
                for (size_t j = 0; j < bricks; ++j)
                {//brick placement

                    x += gable_brick_w;

                    Tm = V::make(x, y, z);

                    create_rigid_body<MT>(engine, Tm, orientation, gable_brick, mid, stone_density);
                }

                y += gable_brick_h;
                --bricks;
            }
        }
    } /// gables

    { /// tiles
        V Tm, Tb, Tw, Tu;
        Q Qm, Qb, Qw, Qu;

            /// roof tiles
        vertices[0] = EigenVector3<T>(0, 0, 0);
        vertices[1] = EigenVector3<T>(VT::numeric_cast(0.4f) * pillar_width, 0, 0);
        vertices[2] = EigenVector3<T>(VT::numeric_cast(1.4f) * pillar_width, gable_brick_h, 0);
        vertices[3] = EigenVector3<T>(pillar_width, gable_brick_h, 0);

        vertices[4]
            = vertices[0]
            - EigenVector3<T>(0, 0, 0.5f * (temple_depth - (floor(beam_spacing_ratio * num_pillars_z) - 1) * spacing));
        vertices[5]
            = vertices[1]
            - EigenVector3<T>(0, 0, 0.5f * (temple_depth - (floor(beam_spacing_ratio * num_pillars_z) - 1) * spacing));
        vertices[6]
            = vertices[2]
            - EigenVector3<T>(0, 0, 0.5f * (temple_depth - (floor(beam_spacing_ratio * num_pillars_z) - 1) * spacing));
        vertices[7]
            = vertices[3]
            - EigenVector3<T>(0, 0, 0.5f * (temple_depth - (floor(beam_spacing_ratio * num_pillars_z) - 1) * spacing));

        GeometryHandleEigen<T> roof_end_tile = create_geometry_handle_cuboid<T>(engine, vertices);

        vertices[0] = EigenVector3<T>(0, 0, 0);
        vertices[1] = EigenVector3<T>(VT::numeric_cast(0.4f) * pillar_width, 0, 0);
        vertices[2] = EigenVector3<T>(VT::numeric_cast(1.4f) * pillar_width, gable_brick_h, 0);
        vertices[3] = EigenVector3<T>(pillar_width, gable_brick_h, 0);
        vertices[4] = vertices[0] - EigenVector3<T>(0, 0, spacing);
        vertices[5] = vertices[1] - EigenVector3<T>(0, 0, spacing);
        vertices[6] = vertices[2] - EigenVector3<T>(0, 0, spacing);
        vertices[7] = vertices[3] - EigenVector3<T>(0, 0, spacing);

        GeometryHandleEigen<T> roof_tile = create_geometry_handle_cuboid<T>(engine, vertices);

            /// regular roof tiles
        for (size_t i = 0u; i < floor(beam_spacing_ratio * num_pillars_z) - 1; ++i)
        {
            T z = 0.5f * temple_depth - 2 * pillar_width - (0.5f + i) * spacing;

            for (size_t i = 0; i < gable_num_layers - 1; ++i)
            {
                T x = -0.5f * temple_width + i * pillar_width;
                T y = beam_y + 0.5f * beam_height + i * gable_brick_h;

                Tm = V::make(x, y, z);
                Qm = Q::identity();

                Tb = fromEigen(roof_tile.Tb2m());
                Qb = fromEigen(roof_tile.Qb2m());

                Tw = rotate(Qm, Tb) + Tm;
                Qw = Qm * Qb;

                Tu = rotate(orientation, Tw) + position;
                Qu = orientation * Qw;

                create_rigid_body<MT>(engine, Tu, Qu, roof_tile, mid, stone_density);

                Tm = V::make(-x, y, -z);
                Qm = Q::Ru(VT::pi(), V::j());

                Tw = rotate(Qm, Tb) + Tm;
                Qw = Qm * Qb;

                Tu = rotate(orientation, Tw) + position;
                Qu = orientation * Qw;

                create_rigid_body<MT>(engine, Tu, Qu, roof_tile, mid, stone_density);
            }
        }
            /// end roof tiles
        for (size_t i = 0; i < gable_num_layers - 1; ++i)
        {
            T x = -0.5f * temple_width + i * pillar_width;
            T y = beam_y + 0.5f * beam_height + i * gable_brick_h;
            T z = 0.5f * temple_depth;

            Tm = V::make(x, y, -z + 0.5f * (temple_depth - (floor(beam_spacing_ratio * num_pillars_z) - 1) * spacing));
            Qm = Q::identity();

            Tb = fromEigen(roof_end_tile.Tb2m());
            Qb = fromEigen(roof_end_tile.Qb2m());

            Tw = rotate(Qm, Tb) + Tm;
            Qw = Qm * Qb;

            Tu = rotate(orientation, Tw) + position;
            Qu = orientation * Qw;

            create_rigid_body<MT>(engine, Tu, Qu, roof_end_tile, mid, stone_density);

            Tm = V::make(x, y, z);

            Tw = rotate(Qm, Tb) + Tm;
            Qw = Qm * Qb;

            Tu = rotate(orientation, Tw) + position;
            Qu = orientation * Qw;

            create_rigid_body<MT>(engine, Tu, Qu, roof_end_tile, mid, stone_density);
            Tm = V::make(-x, y, -z);
            Qm = Q::Ru(VT::pi(), V::j());

            Tw = rotate(Qm, Tb) + Tm;
            Qw = Qm * Qb;

            Tu = rotate(orientation, Tw) + position;
            Qu = orientation * Qw;

            create_rigid_body<MT>(engine, Tu, Qu, roof_end_tile, mid, stone_density);

            Tm = V::make(-x, y, z - 0.5f * (temple_depth - (floor(beam_spacing_ratio * num_pillars_z) - 1) * spacing));

            Tw = rotate(Qm, Tb) + Tm;
            Qw = Qm * Qb;

            Tu = rotate(orientation, Tw) + position;
            Qu = orientation * Qw;

            create_rigid_body<MT>(engine, Tu, Qu, roof_end_tile, mid, stone_density);
        }
    } /// roof tiles

    { /// top triangles

            /// top triangle brick
        vertices[0] = EigenVector3<T>(0, 0, 0);
        vertices[1] = EigenVector3<T>(pillar_width, 0, 0);
        vertices[2] = EigenVector3<T>(VT::numeric_cast(0.51f) * pillar_width, 0.5f * gable_brick_h, 0);
        vertices[3] = EigenVector3<T>(VT::numeric_cast(0.49f) * pillar_width, 0.5f * gable_brick_h, 0);

        vertices[4] = vertices[0] - EigenVector3<T>(0, 0, pillar_width);
        vertices[5] = vertices[1] - EigenVector3<T>(0, 0, pillar_width);
        vertices[6] = vertices[2] - EigenVector3<T>(0, 0, pillar_width);
        vertices[7] = vertices[3] - EigenVector3<T>(0, 0, pillar_width);

        GeometryHandleEigen<T> gable_top_brick = create_geometry_handle_cuboid<T>(engine, vertices);

        for (size_t i = 0; i < 2 * num_pillars_z - 1; ++i)
        {

            Tm = V::make(-0.5f * pillar_width, temple_height - 0.5f * gable_brick_h,
                         0.5f * temple_depth - i * pillar_width);
            Qm = Q::identity();

            Tb = fromEigen(gable_top_brick.Tb2m());
            Qb = fromEigen(gable_top_brick.Qb2m());

            Tw = rotate(Qm, Tb) + Tm;
            Qw = Qm * Qb;

            Tu = rotate(orientation, Tw) + position;
            Qu = orientation * Qw;

            create_rigid_body<MT>(engine, Tu, Qu, gable_top_brick, mid, stone_density);
        }
    } /// top triangles
}

using MTf = tiny::MathTypes<float>;

template void make_temple<MTf>(content::API* engine, MTf::vector3_type const& position,
                               MTf::quaternion_type const& orientation, MTf::real_type const& temple_height,
                               MTf::real_type const& pillar_width, size_t const& num_pillars_x,
                               size_t const& num_pillars_z, MaterialInfo<MTf::real_type> mat_info);

} //namespace procedural
