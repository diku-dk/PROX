#ifndef SOFT_BODY_GUI_DRAW_SCENE_H
#define SOFT_BODY_GUI_DRAW_SCENE_H

#include <soft_body_gui_lighting.h>
#include <soft_body_gui_materials.h>
#include <soft_body_gui_scene_manager.h>

#include <gl3_texture_manager.h>
#include <gl3_shader_program_manager.h>
#include <gl3_vbo.h>
#include <gl3_glm.h>

namespace soft_body
{
  namespace gui
  {

    inline void draw_scene(
                           bool const & draw_wireframe
                           , glm::mat4 const & projection_matrix
                           , glm::mat4 const & view_matrix
                           , std::vector<LightInfo> & lights
                           , std::vector<MaterialInfo> & materials
                           , gl3::TextureManager const & texture_manager
                           , gl3::ShaderProgramManager & program_manager
                           , SceneManager const & scene_manager
                           )
    {
      gl3::check_errors("draw_scene() invoked");

      gl3::Program & solid_program = program_manager.get(0);
      gl3::Program & wire_program  = program_manager.get(1);

      gl3::Program & program = draw_wireframe ? wire_program : solid_program;

      program.use();

      program.set_uniform("projection_matrix", projection_matrix );
      program.set_uniform("view_matrix", view_matrix );

      if(draw_wireframe)
        program.set_uniform("wire_color", 0.9f, 0.1f, 0.1f );  // 2014-10-6 Kenny: Consider having this as hardwired constant?

      set_lighting_uniforms(program, lights);

      std::vector<SceneObject>::const_iterator object = scene_manager.m_objects.begin();
      std::vector<SceneObject>::const_iterator end    = scene_manager.m_objects.end();

      for(;object!=end;++object)
      {
        gl3::Texture3D const & texture = texture_manager.get( object->m_visual_idx );

        glm::mat4 const model_matrix      = glm::mat4(1.0f);   // 2015-03-19 Kenny: T4Meshes live in world coordinate system, no model transform needed.

        glm::mat4 const model_view_matrix =  view_matrix * model_matrix;

        glActiveTexture( GL_TEXTURE0 );       // Active texture unit 0
        texture.bind();                       // Bind this texture to texture unit 0
        program.set_uniform( "solid_texture", 0 ); // solid_texture sampler should use texture unit 0

        program.set_uniform( "model_view_matrix", model_view_matrix);
        program.set_uniform( "texture_matrix", object->m_texture_matrix );

        set_material_uniforms(program, materials[object->m_visual_idx]);

        gl3::VAO const & vao =  draw_wireframe ?  object->m_wire_vao : object->m_solid_vao;

        vao.bind();
        object->m_vbo.draw();
        vao.unbind();
        
        texture.unbind();
      }
      
      program.stop();
    }
    
  }//namespace gui
}//namespace soft_body

// SOFT_BODY_GUI_DRAW_SCENE_H
#endif