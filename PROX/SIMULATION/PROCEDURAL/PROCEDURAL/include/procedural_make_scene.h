#ifndef PROCEDURAL_MAKE_SCENE_H
#define PROCEDURAL_MAKE_SCENE_H

#include <content.h>
#include <procedural.h>
#include <procedural_factory.h>

#include <mesh_array.h>

#include <util_config_file.h>

#include <cstdlib>

namespace procedural
{

  /**
   * This function can be used to quickly generate known test scenes.
   *
   * @param scene        The named scene to create.
   * @param obj_path     The folder location to find obj-files.
   * @param engine       A pointer to the engine that upon return will have added the named scene to it.
   * @param params       Procedural parameter settings for controlling "some" of the name scenes. Parameters
   *                     could be size of scene or number of object, or object types etc..
   */
template <typename T>
inline void make_scene(std::string const& scene, std::string const& obj_path, content::API* engine,
                       util::ConfigFile const& params)
{

    mesh_array::TetGenSettings tetset = mesh_array::tetgen_default_settings();
    tetset.m_quality_ratio      = util::to_value<double>(params.get_value("tetgen_quality_ratio", "2.0"));
    tetset.m_maximum_volume     = util::to_value<double>(params.get_value("tetgen_maximum_volume", "0.0"));
    tetset.m_quiet_output       = util::to_value<bool>(params.get_value("tetgen_quiet_output", "true"));
    tetset.m_suppress_splitting = util::to_value<bool>(params.get_value("tetgen_suppress_splitting", "true"));

    procedural::MaterialInfo<T> mat_info = procedural::create_material_info<T>(engine);

    if (scene.compare("arch") == 0)
    {
        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), mat_info, 10.0f,
                                   1.0f, 10.0f);

        procedural::make_arch<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), 2      // r outer
                                 ,
                                 1.5    // r innter
                                 ,
                                 2      // pillar height
                                 ,
                                 2.0    // stone depth
                                 ,
                                 7      // arch slices
                                 ,
                                 3      // pillar segments
                                 ,
                                 mat_info);
    }
    if (scene.compare("pillar") == 0)
    {
        auto const ground_width = util::to_value<float>(params.get_value("procedural_param_1", "10.0"));
        auto const pillar_height = util::to_value<float>(params.get_value("procedural_param_2", "3.5"));
        auto const pillar_width = util::to_value<float>(params.get_value("procedural_param_3", "0.5"));

        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), mat_info,
                                   ground_width,
                                   ground_width / 10.0  // ground height
                                   ,
                                   ground_width);

        if (tetset.m_maximum_volume > 0.0) tetset.m_maximum_volume /= 1000.0;

        procedural::make_greek_pillar<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), pillar_width,
                                         pillar_height, pillar_width, 3, 6, mat_info, tetset);
    }
    if (scene.compare("twist") == 0)
    {
        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), mat_info, 10.0f,
                                   1.0f, 10.0f);

        auto const layers = util::to_value<unsigned int>(params.get_value("procedural_param_1", "5"));
        auto const degree = util::to_value<float>(params.get_value("procedural_param_2", "25.0"));

        procedural::make_twisted_stack<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), 1, layers,
                                          degree, mat_info);
    }
    if (scene.compare("stack") == 0)
    {
        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), mat_info, 10.0f,
                                   1.0f, 10.0f);

        auto const layers = util::to_value<unsigned int>(params.get_value("procedural_param_1", "5"));

        procedural::make_stack<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), 1, layers,
                                  mat_info);
    }
    if (scene.compare("wall") == 0)
    {
        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), mat_info, 10.0f,
                                   1.0f, 10.0f);

        T const width = util::to_value<T>(params.get_value("procedural_param_1", "5.0"));
        T const height = util::to_value<T>(params.get_value("procedural_param_2", "4.0"));
        T const depth = util::to_value<T>(params.get_value("procedural_param_3", "1.0"));
        auto const layers = util::to_value<unsigned int>(params.get_value("procedural_param_4", "5"));
        auto const span = util::to_value<unsigned int>(params.get_value("procedural_param_5", "5"));

        procedural::make_wall<T>(engine, EigenVector3<T>(-0.5f * width, 0, 0), EigenQuaternion<T>::Identity(), width,
                                 height, depth, layers, span, mat_info);
    }
    if (scene.compare("tower") == 0)
    {
        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), mat_info, 10.0f,
                                   1.0f, 10.0f);

        T const outer_radius = util::to_value<T>(params.get_value("procedural_param_1", "3.0"));
        T const inner_radius = util::to_value<T>(params.get_value("procedural_param_2", "2.5"));
        T const height = util::to_value<T>(params.get_value("procedural_param_3", "3.0"));
        auto const slices = util::to_value<unsigned int>(params.get_value("procedural_param_4", "12"));
        auto const segments = util::to_value<unsigned int>(params.get_value("procedural_param_5", "4"));

        procedural::make_tower<T>(engine, EigenVector3<T>(0, 0, 0), Rotatex(-std::numbers::pi_v<T>), outer_radius,
                                  inner_radius, height, slices, segments, mat_info, true);
    }
    if (scene.compare("cuboid_tower") == 0)
    {
        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), mat_info, 10.0f,
                                   1.0f, 10.0f);

        T const outer_radius = util::to_value<T>(params.get_value("procedural_param_1", "3.0"));
        T const inner_radius = util::to_value<T>(params.get_value("procedural_param_2", "2.5"));
        T const height = util::to_value<T>(params.get_value("procedural_param_3", "3.0"));
        auto const slices = util::to_value<unsigned int>(params.get_value("procedural_param_4", "12"));
        auto const segments = util::to_value<unsigned int>(params.get_value("procedural_param_5", "4"));

        procedural::make_tower<T>(engine, EigenVector3<T>(0, 0, 0), Rotatex(-std::numbers::pi_v<T>), outer_radius,
                                  inner_radius, height, slices, segments, mat_info, false);
    }
    if (scene.compare("dome") == 0)
    {
        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), mat_info, 10.0f,
                                   1.0f, 10.0f);

        T const outer_radius = util::to_value<T>(params.get_value("procedural_param_1", "3.0"));
        T const inner_radius = util::to_value<T>(params.get_value("procedural_param_2", "2.5"));
        auto const slices = util::to_value<unsigned int>(params.get_value("procedural_param_3", "10"));
        auto const segments = util::to_value<unsigned int>(params.get_value("procedural_param_4", "8"));

        procedural::make_dome<T>(engine, EigenVector3<T>(0, 0, 0), Rotatex(-std::numbers::pi_v<T> * 0.5f), outer_radius,
                                 inner_radius, slices, segments, mat_info);
    }
    if (scene.compare("pantheon") == 0)
    {
        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), mat_info, 10.0f,
                                   1.0f, 10.0f);

        T const outer_radius = util::to_value<T>(params.get_value("procedural_param_1", "3.0"));
        T const inner_radius = util::to_value<T>(params.get_value("procedural_param_2", "2.5"));
        T const height = util::to_value<T>(params.get_value("procedural_param_3", "2.0"));
        auto const slices = util::to_value<unsigned int>(params.get_value("procedural_param_4", "11"));
        auto const segments = util::to_value<unsigned int>(params.get_value("procedural_param_5", "11"));

        procedural::make_pantheon<T>(engine, EigenVector3<T>(0, 0, 0), Rotatex(-std::numbers::pi_v<T>), outer_radius,
                                     inner_radius, height, slices, segments, mat_info);
    }
    if (scene.compare("temple") == 0)
    {
        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), mat_info, 10.0f,
                                   1.0f, 10.0f);

        T const temple_height = util::to_value<T>(params.get_value("procedural_param_1", "3.0"));
        T const pillar_width = util::to_value<T>(params.get_value("procedural_param_2", "0.5"));
        auto const cnt_pillars_x = util::to_value<unsigned int>(params.get_value("procedural_param_3", "6"));
        auto const cnt_pillars_y = util::to_value<unsigned int>(params.get_value("procedural_param_4", "8"));

        procedural::make_temple<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), temple_height,
                                   pillar_width, cnt_pillars_x, cnt_pillars_y, mat_info);
    }
    if(scene.compare("colosseum") == 0)
    {
      T const outer_radius                 = util::to_value<T>( params.get_value( "procedural_param_1", "13.0") );
      T const inner_radius                 = util::to_value<T>( params.get_value( "procedural_param_2", "10.0") );
      auto const number_of_arches = util::to_value<unsigned int>(params.get_value("procedural_param_3", "14"));
      auto const layers_of_archs = util::to_value<unsigned int>(params.get_value("procedural_param_4", "3"));

      procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), mat_info,
                                 outer_radius * 3, 1.0f, outer_radius * 3);

      procedural::make_colosseum<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(),
                                    outer_radius // outer radius
                                    ,
                                    inner_radius // inner radius
                                    ,
                                    number_of_arches // number of arches
                                    ,
                                    layers_of_archs // layers of archs
                                    ,
                                    mat_info);
    }
    if (scene.compare("scripted_motions") == 0)
    {
        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), mat_info, 10.0f,
                                   1.0f, 10.0f);

        T const width = util::to_value<T>(params.get_value("procedural_param_1", "5.0"));
        T const layers = util::to_value<T>(params.get_value("procedural_param_2", "5.0"));
        T const span = util::to_value<T>(params.get_value("procedural_param_3", "5.0"));

        procedural::make_wall<T>(engine, EigenVector3<T>(-0.5f * width, 0, 0), EigenQuaternion<T>::Identity(),
                                 width      // width
                                 ,
                                 layers * 2.0f      // height
                                 ,
                                 2.0f       // depth
                                 ,
                                 layers    // layers
                                 ,
                                 span      // span
                                 ,
                                 mat_info);

        size_t const oscillation_motion_idx = engine->create_oscilation_scripted_motion();

        engine->set_scripted_oscilation_paramters(
            oscillation_motion_idx, width * 2.0f // amplitude
            ,
            (std::numbers::pi_v<T> * 2.0) / 2.5f // frequency   omega = 2 pi / T; where T is the period
            ,
            0.0 // phase
            ,
            1.0 // dir x
            ,
            1.0 // dir y
            ,
            0.0 // dir z
            ,
            0.0 // ref x
            ,
            0.0 // ref y
            ,
            0.0 // ref z
        );

        size_t const keyframe_motion_idx = engine->create_key_frame_scripted_motion();

        engine->set_scripted_key_position(keyframe_motion_idx, 0.0f   // time
                                          ,
                                          -width // x
                                          ,
                                          0.0f   // y
                                          ,
                                          0.0f   // z
        );

        engine->set_scripted_key_position(keyframe_motion_idx, 1.0f   // time
                                          ,
                                          width  // x
                                          ,
                                          0.0f   // y
                                          ,
                                          0.0f   // z
        );

        engine->set_scripted_key_position(keyframe_motion_idx, 2.0f        // time
                                          ,
                                          width       // x
                                          ,
                                          2.0f * width  // y
                                          ,
                                          0.0f        // z
        );

        engine->set_scripted_key_position(keyframe_motion_idx, 3.0f        // time
                                          ,
                                          -width      // x
                                          ,
                                          2.0f * width  // y
                                          ,
                                          0.0f        // z
        );

        engine->set_scripted_key_position(keyframe_motion_idx, 4.0f    // time
                                          ,
                                          -width  // x
                                          ,
                                          0.0f    // y
                                          ,
                                          0.0f    // z
        );

        std::vector<size_t> rids;

        rids.resize(engine->get_number_of_rigid_bodies());

        engine->get_rigid_body_indices(&rids[0]);

        engine->connect_scripted_motion(rids[rids.size() - 1u], keyframe_motion_idx);
        engine->connect_scripted_motion(rids[rids.size() - 2u], oscillation_motion_idx);
    }
    if(scene.compare("funnel_dims") == 0)
    {
        auto const funnel_size = util::to_value<float>(params.get_value("procedural_param_1", "1.0"));
        auto const funnel_height = util::to_value<float>(params.get_value("procedural_param_2", "10.0"));
        auto const object_size = util::to_value<float>(params.get_value("procedural_param_3", "1.0"));
        auto const ground_size = util::to_value<float>(params.get_value("procedural_param_4", "10.0"));
        auto const object_spacing = util::to_value<float>(params.get_value("procedural_param_5", "0.1"));
        auto const number_of_objects_in_x = util::to_value<float>(params.get_value("procedural_param_6", "20.0"));
        auto const number_of_objects_in_y = util::to_value<float>(params.get_value("procedural_param_7", "20.0"));
        auto const number_of_objects_in_z = util::to_value<float>(params.get_value("procedural_param_8", "20.0"));

        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), mat_info,
                                   ground_size, 1.0f, ground_size);

        // A funnel is 1.0x0.34x1.0
        // It is also upside-down.
        procedural::make_obj<T>(engine, obj_path + "funnel.obj", funnel_size, EigenVector3<T>(0, funnel_height, 0),
                                Rotatez(std::numbers::pi_v<T>), mat_info, true, false, "Stone", tetset);

        std::vector<std::string> obj_names;
        obj_names.push_back(obj_path + "dims.obj");

        procedural::make_obj_packing<T>(engine, EigenVector3<T>(0.0f, funnel_size + funnel_height, 0.0f),
                                        EigenQuaternion<T>::Identity(), obj_names, number_of_objects_in_x,
                                        number_of_objects_in_y, number_of_objects_in_z, object_size, object_spacing,
                                        mat_info, tetset);
    }
    if (scene.compare("earthquake") == 0)
    {
        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), mat_info, 10.0f,
                                   1.0f, 10.0f);

        T const temple_height = util::to_value<T>(params.get_value("procedural_param_1", "3.0"));
        T const pillar_width = util::to_value<T>(params.get_value("procedural_param_2", "0.5"));
        auto const cnt_pillars_x = util::to_value<unsigned int>(params.get_value("procedural_param_3", "6"));
        auto const cnt_pillars_y = util::to_value<unsigned int>(params.get_value("procedural_param_4", "8"));

        procedural::make_temple<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), temple_height,
                                   pillar_width, cnt_pillars_x, cnt_pillars_y, mat_info);

        size_t const oscillation_motion_idx = engine->create_oscilation_scripted_motion();

        T const amplitude = util::to_value<T>(params.get_value("procedural_param_5", "0.01"));
        T const frequency = util::to_value<T>(params.get_value("procedural_param_6", "3.0"));
        T const phase = util::to_value<T>(params.get_value("procedural_param_7", "0.0"));
        T const dir_x = util::to_value<T>(params.get_value("procedural_param_8", "1.0"));
        T const dir_y = util::to_value<T>(params.get_value("procedural_param_9", "0.0"));
        T const dir_z = util::to_value<T>(params.get_value("procedural_param_10", "0.0"));
        T const ref_x = util::to_value<T>(params.get_value("procedural_param_11", "0.0"));
        T const ref_y = util::to_value<T>(params.get_value("procedural_param_12", "0.0"));
        T const ref_z = util::to_value<T>(params.get_value("procedural_param_13", "0.0"));

        engine->set_scripted_oscilation_paramters(oscillation_motion_idx, amplitude, frequency, phase, dir_x, dir_y,
                                                  dir_z, ref_x, ref_y, ref_z);

        std::vector<size_t> rids;
        rids.resize(engine->get_number_of_rigid_bodies());
        engine->get_rigid_body_indices(&rids[0]);
        engine->connect_scripted_motion(rids[0u], oscillation_motion_idx);
    }
    if(scene.compare("shoot") == 0)
    {
      T     const radius     = util::to_value<T>(     params.get_value( "procedural_param_1", "1.0" ) );
      auto const pos_x = util::to_value<float>(params.get_value("procedural_param_2", "0"));
      auto const pos_y = util::to_value<float>(params.get_value("procedural_param_3", "0"));
      auto const pos_z = util::to_value<float>(params.get_value("procedural_param_4", "0"));

      auto const vel_x = util::to_value<float>(params.get_value("procedural_param_5", "0"));
      auto const vel_y = util::to_value<float>(params.get_value("procedural_param_6", "0"));
      auto const vel_z = util::to_value<float>(params.get_value("procedural_param_7", "10"));

      auto const colosseum_radius = util::to_value<float>(params.get_value("procedural_param_8", "10"));
      auto const colosseum_width = util::to_value<float>(params.get_value("procedural_param_9", "2"));
      auto const colosseum_arches = util::to_value<float>(params.get_value("procedural_param_10", "12"));
      auto const colosseum_layers = util::to_value<float>(params.get_value("procedural_param_11", "3"));

      procedural::make_cannonball<T>(engine, radius, EigenVector3<T>(pos_x, pos_y, pos_z),
                                     EigenQuaternion<T>::Identity(), EigenVector3<T>(vel_x, vel_y, vel_z), mat_info);

      procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), mat_info,
                                 2.0f * (colosseum_radius + colosseum_width) + 2.0f, 1.0f,
                                 2.0f * (colosseum_radius + colosseum_width) + 2.0f);

      procedural::make_colosseum<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(),
                                    colosseum_radius + colosseum_width, colosseum_radius, colosseum_arches,
                                    colosseum_layers, mat_info);
    }
    if (scene.compare("glass_glasses") == 0)
    {
        auto const large_glass_size = util::to_value<float>(params.get_value("procedural_param_1", "10.0"));
        auto const small_glass_size = util::to_value<float>(params.get_value("procedural_param_2", "0.9"));

        auto const number_of_glasses_in_x = util::to_value<unsigned int>(params.get_value("procedural_param_3", "5.0"));
        auto const number_of_glasses_in_y = util::to_value<unsigned int>(params.get_value("procedural_param_4", "5.0"));
        auto const number_of_glasses_in_z = util::to_value<unsigned int>(params.get_value("procedural_param_5", "5.0"));

        auto const glass_spacing = util::to_value<float>(params.get_value("procedural_param_6", "0.1"));
        auto const ground_width = util::to_value<float>(params.get_value("procedural_param_7", "20.0"));

        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), mat_info,
                                   ground_width, 1.0f, ground_width);

        // A glass is 0.89x1.0x0.89
        // Draw glass second because it may be see through and will blend only
        // with the rest of the scene if drawn last.
        procedural::make_obj<T>(engine, obj_path + "glass.obj", large_glass_size,
                                EigenVector3<T>(0, 0.5 * large_glass_size, 0), EigenQuaternion<T>::Identity(), mat_info,
                                true, false, "Stone", tetset);

        std::vector<std::string> obj_names;
        obj_names.push_back(obj_path + "glass.obj");

        procedural::make_obj_packing<T>(engine, EigenVector3<T>(0.0f, large_glass_size * 1.5f, 0.0f),
                                        EigenQuaternion<T>::Identity(), obj_names, number_of_glasses_in_x,
                                        number_of_glasses_in_y, number_of_glasses_in_z, small_glass_size, glass_spacing,
                                        mat_info, tetset);
    }
    if (scene.compare("glass_dims") == 0)
    {
        auto const large_glass_size = util::to_value<float>(params.get_value("procedural_param_1", "10.0"));
        auto const small_dims_size = util::to_value<float>(params.get_value("procedural_param_2", "0.9"));

        auto const number_of_dims_in_x = util::to_value<unsigned int>(params.get_value("procedural_param_3", "5.0"));
        auto const number_of_dims_in_y = util::to_value<unsigned int>(params.get_value("procedural_param_4", "5.0"));
        auto const number_of_dims_in_z = util::to_value<unsigned int>(params.get_value("procedural_param_5", "5.0"));

        auto const dims_spacing = util::to_value<float>(params.get_value("procedural_param_6", "0.1"));
        auto const ground_width = util::to_value<float>(params.get_value("procedural_param_7", "20.0"));

        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), mat_info,
                                   ground_width, 1.0f, ground_width);

        // A glass is 0.89x1.0x0.89
        // Draw glass second because it may be see through and will blend only
        // with the rest of the scene if drawn last.
        procedural::make_obj<T>(engine, obj_path + "glass.obj", large_glass_size,
                                EigenVector3<T>(0, 0.5 * large_glass_size, 0), EigenQuaternion<T>::Identity(), mat_info,
                                true, false, "Stone", tetset);

        std::vector<std::string> obj_names;
        obj_names.push_back(obj_path + "dims.obj");

        procedural::make_obj_packing<T>(engine, EigenVector3<T>(0.0f, large_glass_size * 1.5f, 0.0f),
                                        EigenQuaternion<T>::Identity(), obj_names, number_of_dims_in_x,
                                        number_of_dims_in_y, number_of_dims_in_z, small_dims_size, dims_spacing,
                                        mat_info, tetset);
    }
    if (scene.compare("glass_spheres") == 0)
    {
        auto const glass_scale = util::to_value<float>(params.get_value("procedural_param_1", "1.0"));
        auto const ground_width = util::to_value<float>(params.get_value("procedural_param_2", "10.0"));
        auto const sphere_radius = util::to_value<float>(params.get_value("procedural_param_3", "0.25"));
        int const cube_width = util::to_value<int>(params.get_value("procedural_param_4", "4"));

        procedural::make_sphere_cube<T>(engine, EigenVector3<T>(0, glass_scale + sphere_radius * cube_width, 0),
                                        EigenQuaternion<T>::Identity(), sphere_radius, cube_width, mat_info);

        procedural::make_obj<T>(engine, obj_path + "glass.obj", glass_scale, EigenVector3<T>(0, 0.5 * glass_scale, 0),
                                EigenQuaternion<T>::Identity(), mat_info, true, false, "Cannonball", tetset);

        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), mat_info,
                                   ground_width, 1.0f, ground_width);
    }
    if (scene.compare("wall_pins") == 0)
    {
        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), mat_info, 10.0f,
                                   1.0f, 10.0f);

        T const width = util::to_value<T>(params.get_value("procedural_param_1", "5.0"));
        T const layers = util::to_value<T>(params.get_value("procedural_param_2", "5.0"));
        T const span = util::to_value<T>(params.get_value("procedural_param_3", "5.0"));

        procedural::make_wall<T>(engine, EigenVector3<T>(-0.5f * width, 0, 0), EigenQuaternion<T>::Identity(),
                                 width // width
                                 ,
                                 layers * 2.0f // height
                                 ,
                                 2.0f // depth
                                 ,
                                 layers // layers
                                 ,
                                 span // span
                                 ,
                                 mat_info);

        std::vector<size_t> pin_indices;

        pin_indices.push_back(engine->create_pin_force());
        pin_indices.push_back(engine->create_pin_force());

        engine->set_pin_target(pin_indices[0], 10, 10, 0);
        engine->set_pin_target(pin_indices[1], -10, 10, 0);

        std::vector<size_t> rids;
        rids.resize(engine->get_number_of_rigid_bodies());
        engine->get_rigid_body_indices(&rids[0]);

        for (size_t i = 0u; i < rids.size(); ++i)
        {
            size_t const rid = rids[i];

            if (engine->get_rigid_body_fixed(rid)) continue;

            size_t const pin_idx = pin_indices[i % 2];

            engine->connect_force(rid, pin_idx);
        }
    }
    if (scene.compare("slide") == 0)
    {
        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), mat_info, 10.0f,
                                   1.0f, 10.0f);

        T const degree = util::to_value<T>(params.get_value("procedural_param_1", "25.0"));

        procedural::make_slide<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), degree, mat_info);
    }
    if (scene.compare("bunny_boxes") == 0)
    {
        auto const obstacle_scale = util::to_value<float>(params.get_value("procedural_param_1", "1.0"));
        auto const obstacle_height = util::to_value<float>(params.get_value("procedural_param_2", "10.0"));
        auto const object_size = util::to_value<float>(params.get_value("procedural_param_3", "1.0"));
        auto const ground_size = util::to_value<float>(params.get_value("procedural_param_4", "10.0"));
        auto const object_spacing = util::to_value<float>(params.get_value("procedural_param_5", "0.1"));
        auto const number_of_objects_in_x = util::to_value<float>(params.get_value("procedural_param_6", "20.0"));
        auto const number_of_objects_in_y = util::to_value<float>(params.get_value("procedural_param_7", "20.0"));
        auto const number_of_objects_in_z = util::to_value<float>(params.get_value("procedural_param_8", "20.0"));

        // The bunny is 0.31 x 1.0 x 0.14
        std::string const obstacle_obj_filename = params.get_value("procedural_param_9", obj_path + "bunny.obj");
        std::string const objects_obj_filename = params.get_value("procedural_param_10", obj_path + "box.obj");

        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), mat_info,
                                   ground_size, 1.0f, ground_size);

        procedural::make_obj<T>(engine, obstacle_obj_filename, obstacle_scale, EigenVector3<T>(0, obstacle_height, 0),
                                Rotatez<T>(0), mat_info, false, false, "Cannonball", tetset);

        std::vector<std::string> obj_names;
        obj_names.push_back(objects_obj_filename);

        procedural::make_obj_packing<T>(
            engine,
            EigenVector3<T>(
                0.0f, number_of_objects_in_y * (object_size + object_spacing) + obstacle_scale + obstacle_height, 0.0f),
            EigenQuaternion<T>::Identity(), obj_names, number_of_objects_in_x, number_of_objects_in_y,
            number_of_objects_in_z, object_size, object_spacing, mat_info, tetset);
    }
    if (scene.compare("sphere_cube") == 0)
    {
        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), mat_info, 10.0f,
                                   1.0f, 10.0f);

        T const sphere_radius = util::to_value<T>(params.get_value("procedural_param_1", "0.5"));
        auto const number_of_spheres = util::to_value<unsigned int>(params.get_value("procedural_param_2", "4"));

        procedural::make_sphere_cube<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), sphere_radius,
                                        number_of_spheres, mat_info);
    }
    if (scene.compare("pile") == 0)
    {
      unsigned int const & number_of_objects_in_x =  util::to_value<unsigned int>( params.get_value( "procedural_param_1", "5") );
      unsigned int const & number_of_objects_in_y =  util::to_value<unsigned int>( params.get_value( "procedural_param_2", "5") );
      unsigned int const & number_of_objects_in_z =  util::to_value<unsigned int>( params.get_value( "procedural_param_3", "5") );
      T            const & object_size            =  util::to_value<T>(            params.get_value( "procedural_param_4", "1.0") );
      T            const & spacing                =  util::to_value<T>(            params.get_value( "procedural_param_5", "0.2") );

      std::vector<std::string> obj_names;

      obj_names.push_back(obj_path + "bunny.obj");
      obj_names.push_back(obj_path + "deer.obj");
      obj_names.push_back(obj_path + "cow.obj");

      procedural::make_box_container<T>(engine, EigenVector3<T>(0.0f, -5.0f, 0.0f), EigenQuaternion<T>::Identity(),
                                        25.0f // width x height x depth of space inside box
                                        ,
                                        1.0f, 25.0f, 1.0f // wall thickness
                                        ,
                                        mat_info);

      procedural::make_obj_packing<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), obj_names,
                                      number_of_objects_in_x, number_of_objects_in_y, number_of_objects_in_z,
                                      object_size, spacing // total number of objects
                                      ,
                                      mat_info, tetset);
    }
    if (scene.compare("heavy_light") == 0)
    {
        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), mat_info, 10.0f,
                                   1.0f, 10.0f);

        auto const sphere_radius = util::to_value<float>(params.get_value("procedural_param_1", "1.0"));
        auto const scale = util::to_value<float>(params.get_value("procedural_param_2", "1000.0"));

        procedural::make_heavy_sphere_light_sphere<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(),
                                                      sphere_radius, scale, mat_info);
    }
    if (scene.compare("dropping") == 0)
    {
        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), mat_info, 10.0f,
                                   1.0f, 10.0f);

        auto const sphere_radius = util::to_value<float>(params.get_value("procedural_param_1", "0.15"));
        auto const width = util::to_value<float>(params.get_value("procedural_param_2", "10.0"));
        auto const height = util::to_value<float>(params.get_value("procedural_param_3", "10.0"));
        auto const depth = util::to_value<float>(params.get_value("procedural_param_4", "10.0"));

        procedural::make_dropping_spheres<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(),
                                             sphere_radius, width, height, depth, mat_info);
    }
    if (scene.compare("funnel") == 0)
    {

        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), mat_info, 10.0f,
                                   1.0f, 10.0f);

        procedural::make_obj<T>(engine, obj_path + "funnel.obj", 4.25, EigenVector3<T>(0, 5, 0),
                                Rotatez(std::numbers::pi_v<T>), mat_info, true, false, "Stone", tetset);

        auto const sphere_radius = util::to_value<float>(params.get_value("procedural_param_1", "0.2"));
        auto const spheres_width = util::to_value<unsigned int>(params.get_value("procedural_param_2", "5"));
        auto const spheres_length = util::to_value<unsigned int>(params.get_value("procedural_param_3", "5"));

        procedural::make_sphere_layer<T>(engine, EigenVector3<T>(0, 8, 0), EigenQuaternion<T>::Identity(),
                                         sphere_radius, spheres_width, spheres_length, mat_info);
    }
    if (scene.compare("propella_glass") == 0)
    {
        auto const large_glass_size = util::to_value<float>(params.get_value("procedural_param_1", "10.0"));
        auto const small_props_size = util::to_value<float>(params.get_value("procedural_param_2", "0.9"));
        auto const number_of_props_in_x = util::to_value<unsigned int>(params.get_value("procedural_param_3", "5.0"));
        auto const number_of_props_in_y = util::to_value<unsigned int>(params.get_value("procedural_param_4", "5.0"));
        auto const number_of_props_in_z = util::to_value<unsigned int>(params.get_value("procedural_param_5", "5.0"));
        auto const props_spacing = util::to_value<float>(params.get_value("procedural_param_6", "0.1"));
        auto const ground_width = util::to_value<float>(params.get_value("procedural_param_7", "20.0"));

        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), mat_info,
                                   ground_width, 1.0f, ground_width);

        // A glass is 0.89x1.0x0.89
        // Draw glass second because it may be see through and will blend only
        // with the rest of the scene if drawn last.
        procedural::make_obj<T>(engine, obj_path + "glass.obj", large_glass_size,
                                EigenVector3<T>(0, 0.5 * large_glass_size, 0), EigenQuaternion<T>::Identity(), mat_info,
                                true, false, "Stone", tetset);

        std::vector<std::string> obj_names;
        obj_names.push_back(obj_path + "propella.obj");

        procedural::make_obj_packing<T>(
            engine, EigenVector3<T>(0.0f, large_glass_size * 1.5f + number_of_props_in_y * small_props_size, 0.0f),
            EigenQuaternion<T>::Identity(), obj_names, number_of_props_in_x, number_of_props_in_y, number_of_props_in_z,
            small_props_size, props_spacing, mat_info, tetset);
    }
    if (scene.compare("packing") == 0)
    {
        auto const number_of_spheres = util::to_value<unsigned int>(params.get_value("procedural_param_1", "10"));

        procedural::make_box_container<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(),
                                          5.0f // width x height x depth of space inside box
                                          ,
                                          5.0f, 5.0f, 1.0f // wall thickness
                                          ,
                                          mat_info);

        procedural::make_box_container<T>(engine, EigenVector3<T>(0.0f, -5.0f, 0.0f), EigenQuaternion<T>::Identity(),
                                          25.0f // width x height x depth of space inside box
                                          ,
                                          1.0f, 25.0f, 1.0f // wall thickness
                                          ,
                                          mat_info);

        procedural::make_sphere_packing<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(),
                                           0.1f // minimum sphere radius
                                           ,
                                           0.5f // maximum sphere radius
                                           ,
                                           5.0f // width x height x depth of space to fill up
                                           ,
                                           5.0f, 5.0f, number_of_spheres // total number of sphere objects
                                           ,
                                           mat_info);
    }
    if (scene.compare("sdf_thin_object_drop") == 0)
    {
        auto const scene_size = util::to_value<float>(
            params.get_value("procedural_param_1", "10.0"));

        float const ground_width = scene_size;
        float const ground_height = scene_size / 10.0;
        float const ground_depth = scene_size;

        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0),
                                   EigenQuaternion<T>::Identity(), mat_info,
                                   ground_width, ground_height, ground_depth);

        procedural::make_sdf_thin_object_drop<T>(
            engine, EigenVector3<T>(0.0, 0.0, 0.0),
            EigenQuaternion<T>::Identity(), scene_size, mat_info);
    }
    if (scene.compare("sdf_thin_thin_object_drop") == 0)
    {
        auto const scene_size = util::to_value<float>(
            params.get_value("procedural_param_1", "10.0"));

        float const ground_width = scene_size;
        float const ground_height = scene_size / 10.0;
        float const ground_depth = scene_size;

        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0),
                                   EigenQuaternion<T>::Identity(), mat_info,
                                   ground_width, ground_height, ground_depth);

        procedural::make_sdf_thin_thin_object_drop<T>(
            engine, EigenVector3<T>(0.0, 0.0, 0.0),
            EigenQuaternion<T>::Identity(), scene_size, mat_info);
    }

    if (scene.compare("ccd_thin_thin_object_shoot") == 0)
    {
        auto const scene_size = util::to_value<float>(
            params.get_value("procedural_param_1", "10.0"));

        float const ground_width = scene_size;
        float const ground_height = scene_size / 10.0;
        float const ground_depth = scene_size;

        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0),
                                   EigenQuaternion<T>::Identity(), mat_info,
                                   ground_width, ground_height, ground_depth);

        procedural::make_ccd_thin_thin_object_shoot<T>(
            engine, EigenVector3<T>(0.0, 0.0, 0.0),
            EigenQuaternion<T>::Identity(), scene_size, mat_info);
    }

    if (scene.compare("ccd_cannonball_fast_shoot") == 0)
    {
        auto const scene_size = util::to_value<float>(
            params.get_value("procedural_param_1", "10.0"));

        float const ground_width = scene_size;
        float const ground_height = scene_size / 10.0;
        float const ground_depth = scene_size;

        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0),
                                   EigenQuaternion<T>::Identity(), mat_info,
                                   ground_width, ground_height, ground_depth);

        procedural::make_ccd_cannonball_fast_shoot<T>(
            engine, EigenVector3<T>(0.0, 0.0, 0.0),
            EigenQuaternion<T>::Identity(), scene_size, mat_info);
    }

    if (scene.compare("ccd_drop_fast_thin") == 0)
    {
        auto const scene_size = util::to_value<float>(
            params.get_value("procedural_param_1", "10.0"));

        float const ground_width = scene_size;
        float const ground_height = scene_size / 10.0;
        float const ground_depth = scene_size;

        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0),
                                   EigenQuaternion<T>::Identity(), mat_info,
                                   ground_width, ground_height, ground_depth);

        procedural::make_ccd_drop_fast_thin<T>(
            engine, EigenVector3<T>(0.0, 0.0, 0.0),
            EigenQuaternion<T>::Identity(), scene_size, mat_info);
    }

    if (scene.compare("ccd_cones_fast_towards_eachother") == 0)
    {
        auto const scene_size = util::to_value<float>(
            params.get_value("procedural_param_1", "10.0"));

        float const ground_width = scene_size;
        float const ground_height = scene_size / 10.0;
        float const ground_depth = scene_size;

        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0),
                                   EigenQuaternion<T>::Identity(), mat_info,
                                   ground_width, ground_height, ground_depth);

        procedural::make_ccd_cones_fast_towards_eachother<T>(
            engine, EigenVector3<T>(0.0, 0.0, 0.0),
            EigenQuaternion<T>::Identity(), scene_size, mat_info);
    }

    if (scene.compare("ccd_multiple_bounces_in_one_iteration") == 0)
    {
        auto const scene_size = util::to_value<float>(
            params.get_value("procedural_param_1", "10.0"));

        float const ground_width = scene_size;
        float const ground_height = scene_size / 10.0;
        float const ground_depth = scene_size;

        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0),
                                   EigenQuaternion<T>::Identity(), mat_info,
                                   ground_width, ground_height, ground_depth);

        procedural::make_ccd_multiple_bounces_in_one_iteration<T>(
            engine, EigenVector3<T>(0.0, 0.0, 0.0),
            EigenQuaternion<T>::Identity(), scene_size, mat_info);
    }

    if (scene.compare("ccd_small_space") == 0)
    {
        auto const scene_size = util::to_value<float>(
            params.get_value("procedural_param_1", "10.0"));

        float const ground_width = scene_size;
        float const ground_height = scene_size / 10.0;
        float const ground_depth = scene_size;

        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0),
                                   EigenQuaternion<T>::Identity(), mat_info,
                                   ground_width, ground_height, ground_depth);

        procedural::make_ccd_small_space<T>(
            engine, EigenVector3<T>(0.0, 0.0, 0.0),
            EigenQuaternion<T>::Identity(), scene_size, mat_info);
    }
    if (scene.compare("ccd_large_large_fast") == 0)
    {
        auto const scene_size = util::to_value<float>(
            params.get_value("procedural_param_1", "10.0"));

        float const ground_width = scene_size;
        float const ground_height = scene_size / 10.0;
        float const ground_depth = scene_size;

        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0),
                                   EigenQuaternion<T>::Identity(), mat_info,
                                   ground_width, ground_height, ground_depth);

        procedural::make_ccd_large_large_fast<T>(
            engine, EigenVector3<T>(0.0, 0.0, 0.0),
            EigenQuaternion<T>::Identity(), scene_size, mat_info);
    }

    //======================================
    //|             PORTFOLIO              |
    //======================================

    if (scene.compare("ccd_bunny_bunny") == 0)
    {
        auto const scene_size = util::to_value<float>(
            params.get_value("procedural_param_1", "10.0"));

        float const ground_width = scene_size;
        float const ground_height = scene_size / 10.0;
        float const ground_depth = scene_size;

        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0),
                                   EigenQuaternion<T>::Identity(), mat_info,
                                   ground_width, ground_height, ground_depth);

        T const stone_density
            = get_material_density_eigen<T>(mat_info, "Stone");
        size_t const mid = get_material_id_eigen<T>(mat_info, "Stone");
        procedural::make_ccd_bunny_bunny<T>(
            engine, EigenVector3<T>(0.0, 0.0, 0.0),
            EigenQuaternion<T>::Identity(), scene_size, mat_info);
        {
            GeometryHandleEigen<T> const bunny = create_geometry_handle_obj<T>(
                engine, obj_path + "bunny.obj", 12.0, 12.0, 12.0,
                mesh_array::tetgen_cdt_settings());
            const EigenVector3<T> T_b2m = bunny.Tb2m();
            const EigenQuaternion<T> Q_b2m = bunny.Qb2m();

            const EigenVector3<T> T_m2l = EigenVector3<T>(0, 0.0, 0.0);
            const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

            const EigenVector3<T> T_l2w = EigenVector3<T>(0.0, 5.0, 0.0);
            const EigenQuaternion<T> Q_l2w = EigenQuaternion<T>::Identity();

            EigenVector3<T> T_b2w;
            EigenQuaternion<T> Q_b2w;

            compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l,
                                               T_l2w, Q_l2w, T_b2w, Q_b2w);

            size_t const rid = create_rigid_body<T>(engine, T_b2w, Q_b2w, bunny,
                                                    mid, stone_density);
            engine->set_rigid_body_velocity(rid, 0.0, 20.0, 0.0);
        }
        {
            GeometryHandleEigen<T> const bunny = create_geometry_handle_obj<T>(
                engine, obj_path + "bunny.obj", 12.0, 12.0, 12.0,
                mesh_array::tetgen_cdt_settings());
            const EigenVector3<T> T_b2m = bunny.Tb2m();
            const EigenQuaternion<T> Q_b2m = bunny.Qb2m();

            const EigenVector3<T> T_m2l = EigenVector3<T>(0, 0.0, 0.0);
            const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

            const EigenVector3<T> T_l2w = EigenVector3<T>(0.0, 12.0, 0.0);
            const EigenQuaternion<T> Q_l2w = EigenQuaternion<T>::Identity();

            EigenVector3<T> T_b2w;
            EigenQuaternion<T> Q_b2w;

            compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l,
                                               T_l2w, Q_l2w, T_b2w, Q_b2w);

            size_t const rid = create_rigid_body<T>(engine, T_b2w, Q_b2w, bunny,
                                                    mid, stone_density);
            engine->set_rigid_body_velocity(rid, 0.0, -20.0, 0.0);
        }
    }
    if (scene.compare("ccd_plane_dropped_on_spikes") == 0)
    {
        auto const scene_size = util::to_value<float>(
            params.get_value("procedural_param_1", "10.0"));

        float const ground_width = scene_size;
        float const ground_height = scene_size / 10.0;
        float const ground_depth = scene_size;

        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0),
                                   EigenQuaternion<T>::Identity(), mat_info,
                                   ground_width, ground_height, ground_depth);

        procedural::make_ccd_plane_dropped_on_spikes<T>(
            engine, EigenVector3<T>(0.0, 0.0, 0.0),
            EigenQuaternion<T>::Identity(), scene_size, mat_info);
    }
    if (scene.compare("ccd_spinning_object_shot") == 0)
    {
        auto const scene_size = util::to_value<float>(
            params.get_value("procedural_param_1", "10.0"));

        float const ground_width = scene_size;
        float const ground_height = scene_size / 10.0;
        float const ground_depth = scene_size;

        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0),
                                   EigenQuaternion<T>::Identity(), mat_info,
                                   ground_width, ground_height, ground_depth);

        T const stone_density
            = get_material_density_eigen<T>(mat_info, "Stone");
        size_t const mid = get_material_id_eigen<T>(mat_info, "Stone");
        procedural::make_ccd_bunny_bunny<T>(
            engine, EigenVector3<T>(0.0, 0.0, 0.0),
            EigenQuaternion<T>::Identity(), scene_size, mat_info);
        {
            GeometryHandleEigen<T> const bunny = create_geometry_handle_obj<T>(
                engine, obj_path + "bunny.obj", 12.0, 12.0, 12.0,
                mesh_array::tetgen_cdt_settings());
            const EigenVector3<T> T_b2m = bunny.Tb2m();
            const EigenQuaternion<T> Q_b2m = bunny.Qb2m();

            const EigenVector3<T> T_m2l = EigenVector3<T>(0, 0.0, 0.0);
            const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

            const EigenVector3<T> T_l2w = EigenVector3<T>(-3.7, 3.0, 0.0);
            const EigenQuaternion<T> Q_l2w = EigenQuaternion<T>::Identity();

            EigenVector3<T> T_b2w;
            EigenQuaternion<T> Q_b2w;

            compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l,
                                               T_l2w, Q_l2w, T_b2w, Q_b2w);

            size_t const rid = create_rigid_body<T>(engine, T_b2w, Q_b2w, bunny,
                                                    mid, stone_density);
            //engine->set_rigid_body_velocity(rid, 0.0, 0.0, 0.0);
            engine->set_rigid_body_fixed(rid, true);
        }
        {
            GeometryHandleEigen<T> const bunny = create_geometry_handle_obj<T>(
                engine, obj_path + "knot.obj", 1.0, 1.0, 1.0,
                mesh_array::tetgen_cdt_settings());
            const EigenVector3<T> T_b2m = bunny.Tb2m();
            const EigenQuaternion<T> Q_b2m = bunny.Qb2m();

            const EigenVector3<T> T_m2l = EigenVector3<T>(0, 0.0, 0.0);
            const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

            const EigenVector3<T> T_l2w = EigenVector3<T>(4.0, 4.0, 0.0);
            const EigenQuaternion<T> Q_l2w = EigenQuaternion<T>::Identity();

            EigenVector3<T> T_b2w;
            EigenQuaternion<T> Q_b2w;

            compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l,
                                               T_l2w, Q_l2w, T_b2w, Q_b2w);

            size_t const rid = create_rigid_body<T>(engine, T_b2w, Q_b2w, bunny,
                                                    mid, stone_density);
            engine->set_rigid_body_velocity(rid, -20.0, 1.0, 0.0);
            engine->set_rigid_body_spin(rid, 100.0, 0.0, 0.0);
        }
    }
    if (scene.compare("ccd_tight_space_drop_grazing_wall") == 0)
    {
        auto const scene_size = util::to_value<float>(
            params.get_value("procedural_param_1", "10.0"));

        float const ground_width = scene_size;
        float const ground_height = scene_size / 10.0;
        float const ground_depth = scene_size;

        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0),
                                   EigenQuaternion<T>::Identity(), mat_info,
                                   ground_width, ground_height, ground_depth);

        T const stone_density
            = get_material_density_eigen<T>(mat_info, "Stone");
        size_t const mid = get_material_id_eigen<T>(mat_info, "Stone");
        procedural::make_ccd_tight_space_drop<T>(
            engine, EigenVector3<T>(0.0, 0.0, 0.0),
            EigenQuaternion<T>::Identity(), scene_size, mat_info);
        if (true)
        {
            GeometryHandleEigen<T> const bunny = create_geometry_handle_obj<T>(
                engine, obj_path + "torus.obj", 1.0, 25.0, 1.0,
                mesh_array::tetgen_cdt_settings());
            const EigenVector3<T> T_b2m = bunny.Tb2m();
            const EigenQuaternion<T> Q_b2m = bunny.Qb2m();

            const EigenVector3<T> T_m2l = EigenVector3<T>(0.0, 7.0, 0.0);
            const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

            const EigenVector3<T> T_l2w = EigenVector3<T>(0.0, 0.0, 0.0);
            const EigenQuaternion<T> Q_l2w = EigenQuaternion<T>::Identity();

            EigenVector3<T> T_b2w;
            EigenQuaternion<T> Q_b2w;

            compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l,
                                               T_l2w, Q_l2w, T_b2w, Q_b2w);

            size_t const rid = create_rigid_body<T>(engine, T_b2w, Q_b2w, bunny,
                                                    mid, stone_density);
            //engine->set_rigid_body_velocity(rid, 0.0, 0.0, 0.0);
            engine->set_rigid_body_fixed(rid, true);
        }
        if (true)
        {
            GeometryHandleEigen<T> const bunny = create_geometry_handle_obj<T>(
                engine, obj_path + "sphere.obj", 0.49, 0.49, 0.49,
                mesh_array::tetgen_cdt_settings());
            const EigenVector3<T> T_b2m = bunny.Tb2m();
            const EigenQuaternion<T> Q_b2m = bunny.Qb2m();

            const EigenVector3<T> T_m2l = EigenVector3<T>(0.0, 25.0, 0.0);
            const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

            const EigenVector3<T> T_l2w = EigenVector3<T>(0.0, 0.0, 0.0);
            const EigenQuaternion<T> Q_l2w = EigenQuaternion<T>::Identity();

            EigenVector3<T> T_b2w;
            EigenQuaternion<T> Q_b2w;

            compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l,
                                               T_l2w, Q_l2w, T_b2w, Q_b2w);

            size_t const rid = create_rigid_body<T>(engine, T_b2w, Q_b2w, bunny,
                                                    mid, stone_density);
            engine->set_rigid_body_velocity(rid, 0.0, -20.0, 0.0);
            //engine->set_rigid_body_fixed(rid, true);
        }
    }
    if (scene.compare("ccd_tight_space_drop") == 0)
    {
        auto const scene_size = util::to_value<float>(
            params.get_value("procedural_param_1", "10.0"));

        float const ground_width = scene_size;
        float const ground_height = scene_size / 10.0;
        float const ground_depth = scene_size;

        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0),
                                   EigenQuaternion<T>::Identity(), mat_info,
                                   ground_width, ground_height, ground_depth);

        T const stone_density
            = get_material_density_eigen<T>(mat_info, "Stone");
        size_t const mid = get_material_id_eigen<T>(mat_info, "Stone");
        procedural::make_ccd_tight_space_drop<T>(
            engine, EigenVector3<T>(0.0, 0.0, 0.0),
            EigenQuaternion<T>::Identity(), scene_size, mat_info);
        if (true)
        {
            GeometryHandleEigen<T> const bunny = create_geometry_handle_obj<T>(
                engine, obj_path + "torus.obj", 1.0, 1.0, 1.0,
                mesh_array::tetgen_cdt_settings());
            const EigenVector3<T> T_b2m = bunny.Tb2m();
            const EigenQuaternion<T> Q_b2m = bunny.Qb2m();

            const EigenVector3<T> T_m2l = EigenVector3<T>(0.0, 1.0, 0.0);
            const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

            const EigenVector3<T> T_l2w = EigenVector3<T>(0.0, 0.0, 0.0);
            const EigenQuaternion<T> Q_l2w = EigenQuaternion<T>::Identity();

            EigenVector3<T> T_b2w;
            EigenQuaternion<T> Q_b2w;

            compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l,
                                               T_l2w, Q_l2w, T_b2w, Q_b2w);

            size_t const rid = create_rigid_body<T>(engine, T_b2w, Q_b2w, bunny,
                                                    mid, stone_density);
            //engine->set_rigid_body_velocity(rid, 0.0, 0.0, 0.0);
            engine->set_rigid_body_fixed(rid, true);
        }
        if (true)
        {
            GeometryHandleEigen<T> const bunny = create_geometry_handle_obj<T>(
                engine, obj_path + "sphere.obj", 0.55, 0.55, 0.55,
                mesh_array::tetgen_cdt_settings());
            const EigenVector3<T> T_b2m = bunny.Tb2m();
            const EigenQuaternion<T> Q_b2m = bunny.Qb2m();

            const EigenVector3<T> T_m2l = EigenVector3<T>(0.0, 2.602, 0.0);
            const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

            const EigenVector3<T> T_l2w = EigenVector3<T>(0.0, 0.0, 0.0);
            const EigenQuaternion<T> Q_l2w = EigenQuaternion<T>::Identity();

            EigenVector3<T> T_b2w;
            EigenQuaternion<T> Q_b2w;

            compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l,
                                               T_l2w, Q_l2w, T_b2w, Q_b2w);

            size_t const rid = create_rigid_body<T>(engine, T_b2w, Q_b2w, bunny,
                                                    mid, stone_density);
            engine->set_rigid_body_velocity(rid, 0.0, -20.0, 0.0);
            //engine->set_rigid_body_fixed(rid, true);
        }
    }
    if (scene.compare("ccd_toi_negative") == 0)
    {
        auto const scene_size = util::to_value<float>(
            params.get_value("procedural_param_1", "10.0"));

        float const ground_width = scene_size;
        float const ground_height = scene_size / 10.0;
        float const ground_depth = scene_size;

        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0),
                                   EigenQuaternion<T>::Identity(), mat_info,
                                   ground_width, ground_height, ground_depth);

        procedural::make_ccd_toi_negative<T>(
            engine, EigenVector3<T>(0.0, 0.0, 0.0),
            EigenQuaternion<T>::Identity(), scene_size, mat_info);
    }
    if (scene.compare("ccd_toi_zero_towards_wall") == 0)
    {
        auto const scene_size = util::to_value<float>(
            params.get_value("procedural_param_1", "10.0"));

        float const ground_width = scene_size;
        float const ground_height = scene_size / 10.0;
        float const ground_depth = scene_size;

        /*        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0),
                                   EigenQuaternion<T>::Identity(), mat_info,
                                   ground_width, ground_height, ground_depth);*/

        T const stone_density
            = get_material_density_eigen<T>(mat_info, "Stone");
        size_t const mid = get_material_id_eigen<T>(mat_info, "Stone");
        /*procedural::make_ccd_toi_zero_towards_wall_very_close<T>(
            engine, EigenVector3<T>(0.0, 0.0, 0.0),
            EigenQuaternion<T>::Identity(), scene_size, mat_info);*/
        if (true)
        {
            GeometryHandleEigen<T> const bunny = create_geometry_handle_obj<T>(
                engine, obj_path + "baseSlider.obj", 1.0, 1.0, 1.0,
                mesh_array::tetgen_cdt_settings());
            const EigenVector3<T> T_b2m = bunny.Tb2m();
            const EigenQuaternion<T> Q_b2m = bunny.Qb2m();

            const EigenVector3<T> T_m2l = EigenVector3<T>(0, 0.0, 0.0);
            const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

            const EigenVector3<T> T_l2w = EigenVector3<T>(0.0, 0.0, 0.0);
            const EigenQuaternion<T> Q_l2w = EigenQuaternion<T>::Identity();

            EigenVector3<T> T_b2w;
            EigenQuaternion<T> Q_b2w;

            compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l,
                                               T_l2w, Q_l2w, T_b2w, Q_b2w);

            size_t const rid = create_rigid_body<T>(engine, T_b2w, Q_b2w, bunny,
                                                    mid, stone_density);
            //engine->set_rigid_body_velocity(rid, 0.0, 0.0, 0.0);
            engine->set_rigid_body_fixed(rid, true);
        }
        if (true)
        {
            GeometryHandleEigen<T> const bunny = create_geometry_handle_obj<T>(
                engine, obj_path + "knot.obj", 0.5, 0.5, 0.5,
                mesh_array::tetgen_cdt_settings());
            const EigenVector3<T> T_b2m = bunny.Tb2m();
            const EigenQuaternion<T> Q_b2m = bunny.Qb2m();

            const EigenVector3<T> T_m2l = EigenVector3<T>(0.75, 0.25, 0.0);
            const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

            const EigenVector3<T> T_l2w = EigenVector3<T>(0.0, 0.0, 0.0);
            const EigenQuaternion<T> Q_l2w = EigenQuaternion<T>::Identity();

            EigenVector3<T> T_b2w;
            EigenQuaternion<T> Q_b2w;

            compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l,
                                               T_l2w, Q_l2w, T_b2w, Q_b2w);

            size_t const rid = create_rigid_body<T>(engine, T_b2w, Q_b2w, bunny,
                                                    mid, stone_density);
            engine->set_rigid_body_velocity(rid, -10, 0.0, 0.0);
            //engine->set_rigid_body_fixed(rid, true);
        }
    }
    if (scene.compare("ccd_toi_zero_towards_wall_very_close") == 0)
    {
        auto const scene_size = util::to_value<float>(
            params.get_value("procedural_param_1", "10.0"));

        float const ground_width = scene_size;
        float const ground_height = scene_size / 10.0;
        float const ground_depth = scene_size;

        /*        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0),
                                   EigenQuaternion<T>::Identity(), mat_info,
                                   ground_width, ground_height, ground_depth);*/

        T const stone_density
            = get_material_density_eigen<T>(mat_info, "Stone");
        size_t const mid = get_material_id_eigen<T>(mat_info, "Stone");
        /*procedural::make_ccd_toi_zero_towards_wall_very_close<T>(
            engine, EigenVector3<T>(0.0, 0.0, 0.0),
            EigenQuaternion<T>::Identity(), scene_size, mat_info);*/
        if (true)
        {
            GeometryHandleEigen<T> const bunny = create_geometry_handle_obj<T>(
                engine, obj_path + "baseSlider.obj", 1.0, 1.0, 1.0,
                mesh_array::tetgen_cdt_settings());
            const EigenVector3<T> T_b2m = bunny.Tb2m();
            const EigenQuaternion<T> Q_b2m = bunny.Qb2m();

            const EigenVector3<T> T_m2l = EigenVector3<T>(0, 0.0, 0.0);
            const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

            const EigenVector3<T> T_l2w = EigenVector3<T>(0.0, 0.0, 0.0);
            const EigenQuaternion<T> Q_l2w = EigenQuaternion<T>::Identity();

            EigenVector3<T> T_b2w;
            EigenQuaternion<T> Q_b2w;

            compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l,
                                               T_l2w, Q_l2w, T_b2w, Q_b2w);

            size_t const rid = create_rigid_body<T>(engine, T_b2w, Q_b2w, bunny,
                                                    mid, stone_density);
            //engine->set_rigid_body_velocity(rid, 0.0, 0.0, 0.0);
            engine->set_rigid_body_fixed(rid, true);
        }
        if (true)
        {
            GeometryHandleEigen<T> const bunny = create_geometry_handle_obj<T>(
                engine, obj_path + "knot.obj", 0.5, 0.5, 0.5,
                mesh_array::tetgen_cdt_settings());
            const EigenVector3<T> T_b2m = bunny.Tb2m();
            const EigenQuaternion<T> Q_b2m = bunny.Qb2m();

            const EigenVector3<T> T_m2l = EigenVector3<T>(0.24, 0.285, 0.0);
            const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

            const EigenVector3<T> T_l2w = EigenVector3<T>(0.0, 0.0, 0.0);
            const EigenQuaternion<T> Q_l2w = EigenQuaternion<T>::Identity();

            EigenVector3<T> T_b2w;
            EigenQuaternion<T> Q_b2w;

            compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l,
                                               T_l2w, Q_l2w, T_b2w, Q_b2w);

            size_t const rid = create_rigid_body<T>(engine, T_b2w, Q_b2w, bunny,
                                                    mid, stone_density);
            engine->set_rigid_body_velocity(rid, -15.0, 0.0, 0.0);
            //engine->set_rigid_body_fixed(rid, true);
        }
    }
    if (scene.compare("ccd_two_objects_shot_fast_eachother") == 0)
    {
        auto const scene_size = util::to_value<float>(
            params.get_value("procedural_param_1", "10.0"));

        float const ground_width = scene_size;
        float const ground_height = scene_size / 10.0;
        float const ground_depth = scene_size;

        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0),
                                   EigenQuaternion<T>::Identity(), mat_info,
                                   ground_width, ground_height, ground_depth);

        procedural::make_ccd_two_objects_shot_fast_eachother<T>(
            engine, EigenVector3<T>(0.0, 0.0, 0.0),
            EigenQuaternion<T>::Identity(), scene_size, mat_info);

        T const stone_density
            = get_material_density_eigen<T>(mat_info, "Stone");
        size_t const mid = get_material_id_eigen<T>(mat_info, "Stone");
        procedural::make_ccd_bunny_bunny<T>(
            engine, EigenVector3<T>(0.0, 0.0, 0.0),
            EigenQuaternion<T>::Identity(), scene_size, mat_info);
        {
            GeometryHandleEigen<T> const bunny = create_geometry_handle_obj<T>(
                engine, obj_path + "bunny.obj", 12.0, 12.0, 12.0,
                mesh_array::tetgen_cdt_settings());
            const EigenVector3<T> T_b2m = bunny.Tb2m();
            const EigenQuaternion<T> Q_b2m = bunny.Qb2m();

            const EigenVector3<T> T_m2l = EigenVector3<T>(0, 0.0, 0.0);
            const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

            const EigenVector3<T> T_l2w = EigenVector3<T>(-4.0, 4.0, 0.0);
            const EigenQuaternion<T> Q_l2w = EigenQuaternion<T>::Identity();

            EigenVector3<T> T_b2w;
            EigenQuaternion<T> Q_b2w;

            compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l,
                                               T_l2w, Q_l2w, T_b2w, Q_b2w);

            size_t const rid = create_rigid_body<T>(engine, T_b2w, Q_b2w, bunny,
                                                    mid, stone_density);
            engine->set_rigid_body_velocity(rid, 100.0, 0.0, 0.0);
            //engine->set_rigid_body_fixed(rid, true);
        }
        {
            GeometryHandleEigen<T> const bunny = create_geometry_handle_obj<T>(
                engine, obj_path + "bunny.obj", 12.0, 12.0, 12.0,
                mesh_array::tetgen_cdt_settings());
            const EigenVector3<T> T_b2m = bunny.Tb2m();
            const EigenQuaternion<T> Q_b2m = bunny.Qb2m();

            const EigenVector3<T> T_m2l = EigenVector3<T>(0, 0.0, 0.0);
            const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

            const EigenVector3<T> T_l2w = EigenVector3<T>(4.0, 4.0, 0.0);
            const EigenQuaternion<T> Q_l2w = EigenQuaternion<T>::Identity();

            EigenVector3<T> T_b2w;
            EigenQuaternion<T> Q_b2w;

            compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l,
                                               T_l2w, Q_l2w, T_b2w, Q_b2w);

            size_t const rid = create_rigid_body<T>(engine, T_b2w, Q_b2w, bunny,
                                                    mid, stone_density);
            engine->set_rigid_body_velocity(rid, -100.0, 00.0, 0.0);
        }
    }
    if (scene.compare("ccd_objects_large_spin") == 0)
    {
        auto const scene_size = util::to_value<float>(
            params.get_value("procedural_param_1", "10.0"));

        float const ground_width = scene_size;
        float const ground_height = scene_size / 10.0;
        float const ground_depth = scene_size;

        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0),
                                   EigenQuaternion<T>::Identity(), mat_info,
                                   ground_width, ground_height, ground_depth);

        procedural::make_ccd_two_objects_shot_fast_eachother<T>(
            engine, EigenVector3<T>(0.0, 0.0, 0.0),
            EigenQuaternion<T>::Identity(), scene_size, mat_info);

        T const stone_density
            = get_material_density_eigen<T>(mat_info, "Stone");
        size_t const mid = get_material_id_eigen<T>(mat_info, "Stone");
        procedural::make_ccd_bunny_bunny<T>(
            engine, EigenVector3<T>(0.0, 0.0, 0.0),
            EigenQuaternion<T>::Identity(), scene_size, mat_info);
        {
            GeometryHandleEigen<T> const bunny = create_geometry_handle_obj<T>(
                engine, obj_path + "knot.obj", 4.0, 4.0, 4.0,
                mesh_array::tetgen_cdt_settings());
            const EigenVector3<T> T_b2m = bunny.Tb2m();
            const EigenQuaternion<T> Q_b2m = bunny.Qb2m();

            const EigenVector3<T> T_m2l = EigenVector3<T>(0, 0.0, 0.0);
            const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

            const EigenVector3<T> T_l2w = EigenVector3<T>(-4.0, 4.0, 0.0);
            const EigenQuaternion<T> Q_l2w = EigenQuaternion<T>::Identity();

            EigenVector3<T> T_b2w;
            EigenQuaternion<T> Q_b2w;

            compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l,
                                               T_l2w, Q_l2w, T_b2w, Q_b2w);

            size_t const rid = create_rigid_body<T>(engine, T_b2w, Q_b2w, bunny,
                                                    mid, stone_density);
            engine->set_rigid_body_velocity(rid, 10.0, 0.0, 0.0);
            engine->set_rigid_body_spin(rid, 300.0, 0.0, 0.0);
            //engine->set_rigid_body_fixed(rid, true);
        }
        {
            GeometryHandleEigen<T> const bunny = create_geometry_handle_obj<T>(
                engine, obj_path + "knot.obj", 4.0, 4.0, 4.0,
                mesh_array::tetgen_cdt_settings());
            const EigenVector3<T> T_b2m = bunny.Tb2m();
            const EigenQuaternion<T> Q_b2m = bunny.Qb2m();

            const EigenVector3<T> T_m2l = EigenVector3<T>(0, 0.0, 0.0);
            const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

            const EigenVector3<T> T_l2w = EigenVector3<T>(4.0, 4.0, 0.0);
            const EigenQuaternion<T> Q_l2w = EigenQuaternion<T>::Identity();

            EigenVector3<T> T_b2w;
            EigenQuaternion<T> Q_b2w;

            compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l,
                                               T_l2w, Q_l2w, T_b2w, Q_b2w);

            size_t const rid = create_rigid_body<T>(engine, T_b2w, Q_b2w, bunny,
                                                    mid, stone_density);
            engine->set_rigid_body_velocity(rid, -10.0, 00.0, 0.0);
            engine->set_rigid_body_spin(rid, -300.0, 0.0, 0.0);
        }
    }
    if (scene.compare("ccd_sphere_through_torus") == 0)
    {
        auto const scene_size = util::to_value<float>(
            params.get_value("procedural_param_1", "10.0"));

        float const ground_width = scene_size;
        float const ground_height = scene_size / 10.0;
        float const ground_depth = scene_size;

        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0),
                                   EigenQuaternion<T>::Identity(), mat_info,
                                   ground_width, ground_height, ground_depth);

        T const stone_density
            = get_material_density_eigen<T>(mat_info, "Stone");
        size_t const mid = get_material_id_eigen<T>(mat_info, "Stone");
        procedural::make_ccd_sphere_through_torus<T>(
            engine, EigenVector3<T>(0.0, 0.0, 0.0),
            EigenQuaternion<T>::Identity(), scene_size, mat_info);
        if (true)
        {
            GeometryHandleEigen<T> const bunny = create_geometry_handle_obj<T>(
                engine, obj_path + "torus.obj", 1.0, 1.0, 1.0,
                mesh_array::tetgen_cdt_settings());
            const EigenVector3<T> T_b2m = bunny.Tb2m();
            const EigenQuaternion<T> Q_b2m = bunny.Qb2m();

            const EigenVector3<T> T_m2l = EigenVector3<T>(0.0, 1.0, 0.0);
            const EigenQuaternion<T> Q_m2l
                = Rotatex(std::numbers::pi_v<T> * T(0.5));

            const EigenVector3<T> T_l2w = EigenVector3<T>(0.0, 0.0, 0.0);
            const EigenQuaternion<T> Q_l2w = EigenQuaternion<T>::Identity();

            EigenVector3<T> T_b2w;
            EigenQuaternion<T> Q_b2w;

            compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l,
                                               T_l2w, Q_l2w, T_b2w, Q_b2w);

            size_t const rid = create_rigid_body<T>(engine, T_b2w, Q_b2w, bunny,
                                                    mid, stone_density);
            //engine->set_rigid_body_velocity(rid, 0.0, 0.0, 0.0);
            engine->set_rigid_body_fixed(rid, true);
        }
        if (true)
        {
            GeometryHandleEigen<T> const bunny = create_geometry_handle_obj<T>(
                engine, obj_path + "sphere.obj", 0.49, 0.49, 0.49,
                mesh_array::tetgen_cdt_settings());
            const EigenVector3<T> T_b2m = bunny.Tb2m();
            const EigenQuaternion<T> Q_b2m = bunny.Qb2m();

            const EigenVector3<T> T_m2l = EigenVector3<T>(6.5, 1.0, 0.0);
            const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

            const EigenVector3<T> T_l2w = EigenVector3<T>(0.0, 0.0, 0.0);
            const EigenQuaternion<T> Q_l2w = EigenQuaternion<T>::Identity();

            EigenVector3<T> T_b2w;
            EigenQuaternion<T> Q_b2w;

            compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l,
                                               T_l2w, Q_l2w, T_b2w, Q_b2w);

            size_t const rid = create_rigid_body<T>(engine, T_b2w, Q_b2w, bunny,
                                                    mid, stone_density);
            engine->set_rigid_body_velocity(rid, -250.0, -1.0, 0.0);
            engine->set_gravity_acceleration(0.0);
            //engine->set_rigid_body_fixed(rid, true);
        }
    }
    if (scene.compare("ccd_two_fast_spinning_objects_shot_fast_eachother") == 0)
    {
        auto const scene_size = util::to_value<float>(
            params.get_value("procedural_param_1", "10.0"));

        float const ground_width = scene_size;
        float const ground_height = scene_size / 10.0;
        float const ground_depth = scene_size;

        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0),
                                   EigenQuaternion<T>::Identity(), mat_info,
                                   ground_width, ground_height, ground_depth);

        procedural::make_ccd_two_objects_shot_fast_eachother<T>(
            engine, EigenVector3<T>(0.0, 0.0, 0.0),
            EigenQuaternion<T>::Identity(), scene_size, mat_info);

        T const stone_density
            = get_material_density_eigen<T>(mat_info, "Stone");
        size_t const mid = get_material_id_eigen<T>(mat_info, "Stone");
        procedural::make_ccd_bunny_bunny<T>(
            engine, EigenVector3<T>(0.0, 0.0, 0.0),
            EigenQuaternion<T>::Identity(), scene_size, mat_info);
        {
            GeometryHandleEigen<T> const bunny = create_geometry_handle_obj<T>(
                engine, obj_path + "knot.obj", 4.0, 4.0, 4.0,
                mesh_array::tetgen_cdt_settings());
            const EigenVector3<T> T_b2m = bunny.Tb2m();
            const EigenQuaternion<T> Q_b2m = bunny.Qb2m();

            const EigenVector3<T> T_m2l = EigenVector3<T>(0, 0.0, 0.0);
            const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

            const EigenVector3<T> T_l2w = EigenVector3<T>(-4.0, 4.0, 0.0);
            const EigenQuaternion<T> Q_l2w = EigenQuaternion<T>::Identity();

            EigenVector3<T> T_b2w;
            EigenQuaternion<T> Q_b2w;

            compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l,
                                               T_l2w, Q_l2w, T_b2w, Q_b2w);

            size_t const rid = create_rigid_body<T>(engine, T_b2w, Q_b2w, bunny,
                                                    mid, stone_density);
            engine->set_rigid_body_velocity(rid, 100.0, 0.0, 0.0);
            engine->set_rigid_body_spin(rid, 100.0, 0.0, 0.0);
            //engine->set_rigid_body_fixed(rid, true);
        }
        {
            GeometryHandleEigen<T> const bunny = create_geometry_handle_obj<T>(
                engine, obj_path + "knot.obj", 4.0, 4.0, 4.0,
                mesh_array::tetgen_cdt_settings());
            const EigenVector3<T> T_b2m = bunny.Tb2m();
            const EigenQuaternion<T> Q_b2m = bunny.Qb2m();

            const EigenVector3<T> T_m2l = EigenVector3<T>(0, 0.0, 0.0);
            const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

            const EigenVector3<T> T_l2w = EigenVector3<T>(4.0, 4.0, 0.0);
            const EigenQuaternion<T> Q_l2w = EigenQuaternion<T>::Identity();

            EigenVector3<T> T_b2w;
            EigenQuaternion<T> Q_b2w;

            compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l,
                                               T_l2w, Q_l2w, T_b2w, Q_b2w);

            size_t const rid = create_rigid_body<T>(engine, T_b2w, Q_b2w, bunny,
                                                    mid, stone_density);
            engine->set_rigid_body_velocity(rid, -100.0, 00.0, 0.0);
            engine->set_rigid_body_spin(rid, -100.0, 0.0, 0.0);
        }
    }
    if (scene.compare("ccd_very_large_triangle_pyramid") == 0)
    {
        auto const scene_size = util::to_value<float>(
            params.get_value("procedural_param_1", "10.0"));

        float const ground_width = scene_size;
        float const ground_height = scene_size / 10.0;
        float const ground_depth = scene_size;

        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0),
                                   EigenQuaternion<T>::Identity(), mat_info,
                                   ground_width, ground_height, ground_depth);

        procedural::make_ccd_very_large_triangle_pyramid<T>(
            engine, EigenVector3<T>(0.0, 0.0, 0.0),
            EigenQuaternion<T>::Identity(), scene_size, mat_info);
    }

    // END

    if (scene.compare("sliding_point") == 0)
    {
        auto const scene_size = util::to_value<float>(params.get_value("procedural_param_1", "10.0"));
        bool const use_spike = util::to_value<bool>(params.get_value("procedural_param_2", "true"));
        bool const use_wedge = util::to_value<bool>(params.get_value("procedural_param_3", "true"));

        float const ground_width = scene_size;
        float const ground_height = scene_size / 10.0;
        float const ground_depth = scene_size;

        procedural::make_ground<T>(engine, EigenVector3<T>(0, 0, 0), EigenQuaternion<T>::Identity(), mat_info,
                                   ground_width, ground_height, ground_depth);

        procedural::make_sliding_point<T>(engine, EigenVector3<T>(0.0, 0.0, 0.0), EigenQuaternion<T>::Identity(),
                                          scene_size, mat_info, use_spike, use_wedge);
    }
    if (scene.compare("two_points") == 0)
    {
        auto const scene_size = util::to_value<float>(params.get_value("procedural_param_1", "10.0"));
        bool const use_spike = util::to_value<bool>(params.get_value("procedural_param_2", "true"));
        bool const use_wedge = util::to_value<bool>(params.get_value("procedural_param_3", "true"));
        bool const use_spike_and_wedge = util::to_value<bool>(params.get_value("procedural_param_4", "true"));

        procedural::make_two_points<T>(engine, EigenVector3<T>(0.0, 0.0, 0.0), EigenQuaternion<T>::Identity(),
                                       scene_size, mat_info, use_spike, use_wedge, use_spike_and_wedge);
    }
    if (scene.compare("point_in_crack") == 0)
    {
        auto const scene_size = util::to_value<float>(params.get_value("procedural_param_1", "10.0"));
        bool const use_spike_in_hole = util::to_value<bool>(params.get_value("procedural_param_2", "true"));
        bool const use_spike_in_crack = util::to_value<bool>(params.get_value("procedural_param_3", "true"));
        bool const use_wedge_in_crack = util::to_value<bool>(params.get_value("procedural_param_4", "true"));

        procedural::make_point_in_crack<T>(engine, EigenVector3<T>(0.0, 0.0, 0.0), EigenQuaternion<T>::Identity(),
                                           scene_size, mat_info, obj_path, use_spike_in_hole, use_spike_in_crack,
                                           use_wedge_in_crack);
    }
    if (scene.compare("cliff_edge") == 0)
    {
        auto const scene_size = util::to_value<float>(params.get_value("procedural_param_1", "20.0"));

        procedural::make_cliff_edge<T>(engine, EigenVector3<T>(0.0, 0.0, 0.0), EigenQuaternion<T>::Identity(),
                                       scene_size, mat_info);
    }
    if (scene.compare("internal_edge") == 0)
    {
        auto const scene_size = util::to_value<float>(params.get_value("procedural_param_1", "20.0"));

        procedural::make_internal_edge<T>(engine, EigenVector3<T>(0.0, 0.0, 0.0), EigenQuaternion<T>::Identity(),
                                          scene_size, mat_info, tetset);
    }
  }

}//namespace procedural


// PROCEDURAL_MAKE_SCENE_H
#endif
