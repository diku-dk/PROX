#ifndef RIGID_BODY_GUI_CHALK_MAKE_PROCEDURAL_CONTENT_H
#define RIGID_BODY_GUI_CHALK_MAKE_PROCEDURAL_CONTENT_H

#include <util_config_file.h>

#include <simulators.h>
#include <procedural.h>
#include <mesh_array.h>

#include <cstdlib>

namespace rigid_body
{
  namespace gui
  {
    namespace chalk
    {

      // Units are in voxels
      // Voxel size is either 25nm or 50nm
      // Images from which data is generated from 2056x2056x256 voxels



      template<typename MT>
      inline void make_procedural_content(content::API * engine, util::ConfigFile const & params)
      {
        typedef typename MT::real_type       T;
        typedef typename MT::vector3_type    V;
        typedef typename MT::quaternion_type Q;

        mesh_array::TetGenSettings tetset = mesh_array::tetgen_default_settings();
        tetset.m_quality_ratio      = util::to_value<double>(params.get_value("tetgen_quality_ratio", "2.0"));
        tetset.m_maximum_volume     = util::to_value<double>(params.get_value("tetgen_maximum_volume", "0.0"));
        tetset.m_quiet_output       = util::to_value<bool>(params.get_value("tetgen_quiet_output", "true"));
        tetset.m_suppress_splitting = util::to_value<bool>(params.get_value("tetgen_suppress_splitting", "true"));

        procedural::MaterialInfo<T> mat_info = procedural::create_material_info<MT>(engine);

        std::string  const grain_data_file = params.get_value( "procedural_param_1", "chalk/ellipsoid_50_hard_final.txt");

        // Get voxel size in meters [m], voxel size is the unit of the input data from the file
        float        const voxel_size              =  util::to_value<float>(        params.get_value( "procedural_param_2", "50e-9")  );
        float        const grain_scale             =  util::to_value<float>(        params.get_value( "procedural_param_3", "400000") );
        unsigned int const total_number_of_grains  =  util::to_value<unsigned int>( params.get_value( "procedural_param_4", "1200")   );
        unsigned int const number_of_grains_in_x   =  util::to_value<unsigned int>( params.get_value( "procedural_param_5", "10")     );
        unsigned int const number_of_grains_in_y   =  util::to_value<unsigned int>( params.get_value( "procedural_param_6", "5O")      );
        unsigned int const number_of_grains_in_z   =  util::to_value<unsigned int>( params.get_value( "procedural_param_7", "10")     );
        float        const grain_density           =  util::to_value<float>(        params.get_value( "procedural_param_8", "1000")   );

        float grain_size = 1.0f;

        procedural::make_grain_packing<MT>(
                                           engine
                                           , V::zero()
                                           , Q::identity()
                                           , grain_data_file
                                           , voxel_size
                                           , grain_scale
                                           , number_of_grains_in_x
                                           , number_of_grains_in_z
                                           , total_number_of_grains
                                           , grain_density
                                           , mat_info
                                           , grain_size
                                           );

        float const container_width  = 1.0f * number_of_grains_in_x  * grain_size;
        float const container_height = 1.0f * number_of_grains_in_y  * grain_size;
        float const container_depth  = 1.0f * number_of_grains_in_z  * grain_size;
        float const wall_thickness   = 1.0f * grain_size;

        procedural::make_box_container<MT>(
                                           engine
                                           , V::make( 0.0f, 0.0f, 0.0f )
                                           , Q::identity()
                                           , container_width
                                           , container_height
                                           , container_depth
                                           , wall_thickness
                                           , mat_info
                                           );

        bool         const oscillation_on  =  util::to_value<bool>(    params.get_value( "procedural_param_9",  "false") );
        float        const amplitude_vx  =  util::to_value<float>(     params.get_value( "procedural_param_10",  "1.0")   );  // Amplitude in voxel units
        float        const frequency  =  util::to_value<float>(        params.get_value( "procedural_param_11", "6.0")   );  // Frequency, 2 pi / T where T is the period of the oscillation
        float        const phase      =  util::to_value<float>(        params.get_value( "procedural_param_12", "0")     );  // Phase shift
        float        const dir_x      =  util::to_value<float>(        params.get_value( "procedural_param_13", "0")     );  // x-component of direction of oscilation
        float        const dir_y      =  util::to_value<float>(        params.get_value( "procedural_param_14", "1")     );  // y-component of direction of oscilation
        float        const dir_z      =  util::to_value<float>(        params.get_value( "procedural_param_15", "0")     );  // z-component of direction of oscilation



        if(oscillation_on)
        {
          float const amplitude = amplitude_vx * voxel_size * grain_scale; // Convert amplitude into proper world units.

          std::vector<size_t> rids;

          rids.resize(engine->get_number_of_rigid_bodies());

          engine->get_rigid_body_indices( &rids[0] );

          for(unsigned int i = 0u; i < rids.size(); ++i)
          {
            unsigned int const & rid = rids[i];

            if(!engine->get_rigid_body_fixed(rid))
              continue;

            float ref_x = 0.0f;
            float ref_y = 0.0f;
            float ref_z = 0.0f;

            engine->get_rigid_body_position(rid, ref_x, ref_y, ref_z);

            size_t const oscillation_motion_idx =  engine->create_oscilation_scripted_motion();

            engine->set_scripted_oscilation_paramters(
                                                      oscillation_motion_idx
                                                      , amplitude
                                                      , frequency
                                                      , phase
                                                      , dir_x
                                                      , dir_y
                                                      , dir_z
                                                      , ref_x
                                                      , ref_y
                                                      , ref_z
                                                      );
            
            engine->connect_scripted_motion(rid, oscillation_motion_idx);
          }
        }
        
      }
      
    }//namespace chalk
  }//namespace gui
}//namespace rigid_body

// RIGID_BODY_GUI_CHALK_MAKE_PROCEDURAL_CONTENT_H
#endif
