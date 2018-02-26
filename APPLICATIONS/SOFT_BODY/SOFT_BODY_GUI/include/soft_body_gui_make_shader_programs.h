#ifndef SOFT_BODY_GUI_MAKE_SHADER_PROGRAMS_H
#define SOFT_BODY_GUI_MAKE_SHADER_PROGRAMS_H

#include <gl3_shader_program_manager.h>

namespace soft_body
{
  namespace gui
  {

    inline void make_shader_programs(
                                     std::string const & shader_path
                                     , gl3::ShaderProgramManager & program_manager
                                     )
    {
      gl3::check_errors("make_shaders(): invoked");

      {
        std::string const vertex_shader_file   = shader_path + "soft_body_gui_solid_vertex.glsl";
        std::string const fragment_shader_file = shader_path + "soft_body_gui_solid_fragment.glsl";

        program_manager.load(vertex_shader_file, fragment_shader_file);
      }

      {
        std::string const vertex_shader_file   = shader_path + "soft_body_gui_wire_vertex.glsl";
        std::string const geometry_shader_file = shader_path + "soft_body_gui_wire_geometry.glsl";
        std::string const fragment_shader_file = shader_path + "soft_body_gui_wire_fragment.glsl";

        program_manager.load(vertex_shader_file, geometry_shader_file, fragment_shader_file);
      }

      {
        std::string const vertex_shader_file   = shader_path + "soft_body_gui_bluish_vertex.glsl";
        std::string const fragment_shader_file = shader_path + "soft_body_gui_bluish_fragment.glsl";

        program_manager.load(vertex_shader_file, fragment_shader_file);
      }

      {
        std::string const vertex_shader_file   = shader_path + "soft_body_gui_aabb_vertex.glsl";
        std::string const fragment_shader_file = shader_path + "soft_body_gui_aabb_fragment.glsl";

        program_manager.load(vertex_shader_file, fragment_shader_file);
      }

      gl3::check_errors("make_shaders(): leaving");
    }
    
  }//namespace gui
}//namespace soft_body

// SOFT_BODY_GUI_MAKE_SHADER_PROGRAMS_H
#endif
