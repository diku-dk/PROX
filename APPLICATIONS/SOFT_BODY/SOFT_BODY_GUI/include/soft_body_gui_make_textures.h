#ifndef SOFT_BODY_GUI_MAKE_TEXTURES_H
#define SOFT_BODY_GUI_MAKE_TEXTURES_H

#include <soft_body_gui_ucph_colors.h>
#include <soft_body_gui_materials.h>

#include <util_config_file.h>

#include <gl3_texture_manager.h>

#include <vector>

namespace soft_body
{
  namespace gui
  {

    template<typename MT>
    inline void make_textures(
                              std::string const & texture_path
                              , hyper::Engine<MT> & engine
                              , util::ConfigFile const & params
                              , std::vector<MaterialInfo> & materials
                              , gl3::TextureManager & texture_manager
                              , std::vector<glm::vec3> & texture_scale_params
                              )
    {
      gl3::check_errors("make_textures(): invoked");

      std::vector<std::string> default_textures;

      default_textures.push_back("dune.vol");
      default_textures.push_back("jagnow.vol");
      default_textures.push_back("cobblestone.vol");
      default_textures.push_back("520_tor.vol");
      default_textures.push_back("animalskin_closeup.vol");
      default_textures.push_back("fauxfur.vol");
      default_textures.push_back("64bluebrown.vol");
      default_textures.push_back("brickwall2.vol");
      default_textures.push_back("64cracks.vol");
      default_textures.push_back("cabin_firewood.vol");
      default_textures.push_back("prtpca2.vol");
      default_textures.push_back("Nmarble23.vol");
      default_textures.push_back("caustic.vol");
      default_textures.push_back("tomatoes.vol");
      default_textures.push_back("Nmarble24.vol");
      default_textures.push_back("waves.vol");
      default_textures.push_back("RGran_Cool_1.vol");
      default_textures.push_back("corral.vol");
      default_textures.push_back("woodwall.vol");
      default_textures.push_back("Rmetal19.vol");
      default_textures.push_back("zebra.vol");

      std::vector<size_t> visual_indices;
      visual_indices.resize(engine.get_number_of_visual_materials());

      engine.get_visual_material_indices( &visual_indices[0] );

      texture_scale_params.clear();

      for(size_t i = 0u; i < visual_indices.size(); ++i)
      {
        size_t      const visual_index = visual_indices[i];
        std::string const name         = engine.get_visual_material_name( visual_index );

        std::string const texture_file = texture_path + params.get_value(name,  default_textures[i%21] );

        texture_manager.load(texture_file);

        gl3::Texture3D & tex = texture_manager.get(i);

        tex.bind();
        tex.set_repeating(true);
        tex.unbind();

        std::vector<std::string> values = params.get_values(name + "_scale" );

        if(values.empty())
          texture_scale_params.push_back( glm::vec3(1.0,1.0,1.0) );
        else
        {
          float sx = util::to_value<float>( values[0] );
          float sy = util::to_value<float>( values[1] );
          float sz = util::to_value<float>( values[2] );
          texture_scale_params.push_back( glm::vec3(sx,sy,sz) );
        }

        MaterialInfo const material = make_material_from_config_file(name, params);
        
        materials.push_back(material);
      }    
    }
    
  }//namespace gui
}//namespace soft_body

// SOFT_BODY_GUI_MAKE_TEXTURES_H
#endif