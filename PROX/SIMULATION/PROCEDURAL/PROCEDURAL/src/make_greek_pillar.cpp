#include <procedural_factory.h>

#include <eigenhelperall.h>

namespace procedural
{

template <typename T>
void make_greek_pillar(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
                       const T& pillar_width, const T& pillar_height, const T& pillar_depth,
                       size_t const& pillar_segments, size_t const& pillar_slices, MaterialInfo<T> mat_info,
                       mesh_array::TetGenSettings tetset)
{

    T const stone_density = get_material_density_eigen<T>(mat_info, "Stone");
    size_t const mid = get_material_id_eigen<T>(mat_info, "Stone");

    std::vector<EigenVector3<T>> segment_vertices(pillar_slices * 2u + 2u);

        /// pillar element height ratios, must sum to one
    T const b_box_ratio = (0.06f);
    T const b_conical_ratio = (0.06f);
    T const c_pillar_ratio = (0.80f);
    T const t_conical_ratio = (0.05f);
    T const t_box_ratio = (0.03f);

    // element heigths
    T const b_box_height = pillar_height * b_box_ratio;
    T const b_conical_height = pillar_height * b_conical_ratio;
    T const c_pillar_height = pillar_height * c_pillar_ratio;
    T const t_conical_height = pillar_height * t_conical_ratio;
    T const t_box_height = pillar_height * t_box_ratio;

    T y = b_box_height * 0.5f;
    {
        GeometryHandleEigen<T> bottom_box
            = create_geometry_handle_box<T>(engine, pillar_width, b_box_height, pillar_depth);

        const EigenVector3<T> local_translation = EigenVector3<T>(0, y, 0);

        EigenVector3<T> body_to_world_translation;
        EigenQuaternion<T> body_to_world_orientation;

        compute_body_to_world_transform<T>(bottom_box.Tb2m()    // body to model
                                           ,
                                           bottom_box.Qb2m()    // body to model
                                           ,
                                           local_translation   // model to local
                                           ,
                                           EigenQuaternion<T>::Identity()       // model to local
                                           ,
                                           (position)            // local to wolrd
                                           ,
                                           (orientation)         // local to world
                                           ,
                                           body_to_world_translation, body_to_world_orientation);

        create_rigid_body<T>(engine, body_to_world_translation, body_to_world_orientation, bottom_box, mid,
                             stone_density);
    }

    y += 0.5f * b_box_height;
    {

        GeometryHandleEigen<T> pillar_segment = create_geometry_handle_pillar_segment<T>(
            engine, (0.5f * pillar_width), ((0.5f * pillar_width) * .8f), b_conical_height, pillar_slices, tetset);

        const EigenVector3<T> local_translation = EigenVector3<T>(0, y, 0);

        EigenVector3<T> body_to_world_translation;
        EigenQuaternion<T> body_to_world_orientation;

        compute_body_to_world_transform<T>((pillar_segment.Tb2m())       // body 2 model
                                           ,
                                           (pillar_segment.Qb2m())       // body 2 model
                                           ,
                                           local_translation           // model 2 local
                                           ,
                                           EigenQuaternion<T>::Identity()               // model 2 local
                                           ,
                                           (position)                    // local 2 world
                                           ,
                                           (orientation)                 // local 2 world
                                           ,
                                           body_to_world_translation, body_to_world_orientation);

        create_rigid_body<T>(engine, body_to_world_translation, body_to_world_orientation, pillar_segment, mid,
                             stone_density);
    }

    y += b_conical_height;

    {

        T const pillar_segment_height = c_pillar_height / pillar_segments;
        T const s_bottom_radius = (0.8f) * (pillar_width * 0.5f);
        T top_radius = (0.6f) * (pillar_width * 0.5f);
        T bottom_radius = s_bottom_radius;
        T const alpha = (top_radius - bottom_radius) / c_pillar_height;

        for (size_t i = 0; i < pillar_segments; ++i)
        {
            top_radius = s_bottom_radius + alpha * pillar_segment_height * (1 + i);

            GeometryHandleEigen<T> c_pillar_segment = create_geometry_handle_pillar_segment<T>(
                engine, bottom_radius, top_radius, pillar_segment_height, pillar_slices, tetset);

            const EigenVector3<T> local_translation = EigenVector3<T>(0, y, 0);

            EigenVector3<T> body_to_world_translation;
            EigenQuaternion<T> body_to_world_orientation;

            compute_body_to_world_transform<T>((c_pillar_segment.Tb2m())  // body to model
                                               ,
                                               (c_pillar_segment.Qb2m())  // body to model
                                               ,
                                               local_translation        // model to local
                                               ,
                                               EigenQuaternion<T>::Identity()            // model to local
                                               ,
                                               (position)                 // local to world
                                               ,
                                               (orientation)              // local to world
                                               ,
                                               body_to_world_translation, body_to_world_orientation);

            create_rigid_body<T>(engine, body_to_world_translation, body_to_world_orientation, c_pillar_segment, mid,
                                 stone_density);

            bottom_radius = top_radius;
            y += pillar_segment_height;
        }
    }

    y += t_conical_height;
    {

        GeometryHandleEigen<T> pillar_segment = create_geometry_handle_pillar_segment<T>(
            engine, ((0.5f * pillar_width)), (0.6f * (0.5f * pillar_width)), t_conical_height, pillar_slices, tetset);

        const EigenVector3<T> local_translation = EigenVector3<T>(0, y, 0);

        EigenVector3<T> body_to_world_translation;
        EigenQuaternion<T> body_to_world_orientation;

        compute_body_to_world_transform<T>((pillar_segment.Tb2m())      // body to model
                                           ,
                                           (pillar_segment.Qb2m())      // body to model
                                           ,
                                           local_translation          // model to local
                                           ,
                                           Rotateu(std::numbers::pi_v<T>, EigenVector3<T>(1, 0, 0))  // model to local
                                           ,
                                           (position)                   // local to world
                                           ,
                                           (orientation)                // local to world
                                           ,
                                           body_to_world_translation, body_to_world_orientation);

        create_rigid_body<T>(engine, body_to_world_translation, body_to_world_orientation, pillar_segment, mid,
                             stone_density);
    }

    y += 0.5f * t_box_height;
    {
        GeometryHandleEigen<T> top_box
            = create_geometry_handle_box<T>(engine, pillar_width, t_box_height, pillar_depth);

        const EigenVector3<T> local_translation = EigenVector3<T>(0, y, 0);

        EigenVector3<T> body_to_world_translation;
        EigenQuaternion<T> body_to_world_orientation;

        compute_body_to_world_transform<T>(top_box.Tb2m()              // body to model
                                           ,
                                           top_box.Qb2m()              // body to model
                                           ,
                                           local_translation           // model to local
                                           ,
                                           EigenQuaternion<T>::Identity()               // model to local
                                           ,
                                           (position)                    // local to world
                                           ,
                                           (orientation)                 // local to world
                                           ,
                                           body_to_world_translation, body_to_world_orientation);

        create_rigid_body<T>(engine, body_to_world_translation, body_to_world_orientation, top_box, mid, stone_density);
    }
}

template void make_greek_pillar<float>(content::API* engine, const EigenVector3<float>& position,
                                       const EigenQuaternion<float>& orientation, const float& pillar_width,
                                       const float& pillar_height, const float& pillar_depth,
                                       size_t const& pillar_segments, size_t const& pillar_slices,
                                       MaterialInfo<float> mat_info, mesh_array::TetGenSettings tetset);

} //namespace procedural
