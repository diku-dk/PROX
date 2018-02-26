#ifndef SOFT_BODY_GUI_DRAW_AABBS_H
#define SOFT_BODY_GUI_DRAW_AABBS_H

#include <soft_body_gui_scene_manager.h>

#include <gl3_shader_program_manager.h>
#include <gl3_glm.h>

#include <hyper_engine.h>

namespace soft_body
{
  namespace gui
  {

    inline void make_aabb_geometry(
                                   gl3::ShaderProgramManager & program_manager
                                   , SceneObject & aabb_geometry
                                   )
    {
      gl3::Program & program = program_manager.get(3);

      aabb_geometry.m_vbo        = gl3::make_unit_aabb_vbo();
      aabb_geometry.m_solid_vao  = gl3::make_vao(aabb_geometry.m_vbo, program, "position", "normal");
    }

    template<typename MT>
    inline void draw_aabbs(
                           glm::mat4 const & projection_matrix
                           , glm::mat4 const & view_matrix
                           , gl3::ShaderProgramManager  & program_manager
                           , SceneObject & aabb_geometry
                           , SceneManager  & scene_manager
                           , hyper::Engine<MT> & engine
                           )
    {
      gl3::check_errors("draw_aabbs() invoked");

      gl3::Program & program = program_manager.get(3);

      program.use();

      program.set_uniform("projection_matrix", projection_matrix );

      std::vector<SceneObject>::const_iterator object = scene_manager.m_objects.begin();
      std::vector<SceneObject>::const_iterator end    = scene_manager.m_objects.end();

      for(;object!=end;++object)
      {
        typename hyper::Engine<MT>::body_type const & body = engine.get_body(object->m_body_idx);

        float const center_x = 0.5f * (body.m_max_point(0) + body.m_min_point(0));
        float const center_y = 0.5f * (body.m_max_point(1) + body.m_min_point(1));
        float const center_z = 0.5f * (body.m_max_point(2) + body.m_min_point(2));

        glm::mat4 const model_matrix = glm::translate(
                                                      glm::mat4(1.0)
                                                      , glm::vec3(
                                                                  center_x
                                                                  , center_y
                                                                  , center_z
                                                                  )
                                                      );

        float const aabb_width  = (body.m_max_point(0) - body.m_min_point(0));
        float const aabb_height = (body.m_max_point(1) - body.m_min_point(1));
        float const aabb_depth  = (body.m_max_point(2) - body.m_min_point(2));

        glm::mat4 const model_view_matrix =  view_matrix * model_matrix;

        program.set_uniform( "scale", aabb_width, aabb_height, aabb_depth );

        program.set_uniform( "color", object->m_red, object->m_green, object->m_blue );
        program.set_uniform( "model_view_matrix", model_view_matrix);
        
        aabb_geometry.m_solid_vao.bind();
        aabb_geometry.m_vbo.draw();
        aabb_geometry.m_solid_vao.unbind();
      }
      
      program.stop();
    }
    
  }//namespace gui
}//namespace soft_body

// SOFT_BODY_GUI_DRAW_AABBS_H
#endif