#ifndef SOFT_BODY_GUI_MAKE_SCENE_H
#define SOFT_BODY_GUI_MAKE_SCENE_H

#include <soft_body_gui_ucph_colors.h>
#include <soft_body_gui_make_textures.h>
#include <soft_body_gui_make_shader_programs.h>
#include <soft_body_gui_scene_manager.h>
#include <soft_body_gui_update_scene.h>
#include <soft_body_gui_materials.h>

#include <util_config_file.h>
#include <gl3_texture_manager.h>
#include <gl3_shader_program_manager.h>
#include <gl3_vbo.h>
#include <tiny_math_types.h>
#include <mesh_array_t3mesh.h>
#include <mesh_array_vertex_attribute.h>
#include <hyper_engine.h>

namespace soft_body
{
  namespace gui
  {

    template<typename MT>
    inline void make_scene(
                           util::ConfigFile const & params
                           , hyper::Engine<MT> & engine
                           , std::vector<MaterialInfo> & materials
                           , std::string const & texture_path
                           , gl3::TextureManager & texture_manager
                           , std::string const & shader_path
                           , gl3::ShaderProgramManager & program_manager
                           , SceneManager & scene_manager
                           )
    {
      gl3::check_errors("make_scene(): invoked");

      make_shader_programs(shader_path, program_manager);

      std::vector<glm::vec3> texture_scale_params;

      make_textures(texture_path, engine, params, materials, texture_manager, texture_scale_params);

      gl3::Program & solid_program = program_manager.get(0);
      gl3::Program & wire_program  = program_manager.get(1);

      typename hyper::Engine<MT>::body_iterator body = engine.body_begin();
      typename hyper::Engine<MT>::body_iterator end  = engine.body_end();

      for(;body!=end;++body)
      {
        SceneObject object;

        object.m_body_idx = body->get_idx();

        size_t      const visual_idx    = engine.get_body(object.m_body_idx).get_visual_idx();
        std::string const name          = engine.get_visual_material_name(  visual_idx );

        object.m_visual_idx = visual_idx;

        {
          object.m_red   = UCPHGreenColor::red();
          object.m_green = UCPHGreenColor::green();
          object.m_blue  = UCPHGreenColor::blue();
        }
        //      {
        //        object.m_red   = UCPHBlueColor::red();
        //        object.m_green = UCPHBlueColor::green();
        //        object.m_blue  = UCPHBlueColor::blue();
        //      }
        //      {
        //        object.m_red   = UCPHRedColor::red();
        //        object.m_green = UCPHRedColor::green();
        //        object.m_blue  = UCPHRedColor::blue();
        //      }

        object.m_texture_matrix = glm::scale( texture_scale_params[visual_idx] );

        mesh_array::T3Mesh mesh;
        mesh_array::VertexAttribute<float, mesh_array::T3Mesh> X;
        mesh_array::VertexAttribute<float, mesh_array::T3Mesh> Y;
        mesh_array::VertexAttribute<float, mesh_array::T3Mesh> Z;

        mesh_array::make_t3mesh( body->m_mesh, body->m_X, body->m_Y, body->m_Z, mesh, X, Y, Z );

        object.m_vbo       = gl3::make_vbo(mesh, X, Y, Z, gl3::DYNAMIC_VBO() );
        object.m_solid_vao = gl3::make_vao(object.m_vbo, solid_program, "position", "normal");
        object.m_wire_vao  = gl3::make_vao(object.m_vbo, wire_program, "position", "normal");
        
        scene_manager.add( object );
      }
      
      update_scene( scene_manager, engine );
    }
    
  }//namespace gui
}//namespace soft_body

// SOFT_BODY_GUI_MAKE_SCENE_H
#endif
