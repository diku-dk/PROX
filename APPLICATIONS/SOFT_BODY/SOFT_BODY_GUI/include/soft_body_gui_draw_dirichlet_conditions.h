#ifndef SOFT_BODY_GUI_DRAW_DIRICHLET_CONDITIONS_H
#define SOFT_BODY_GUI_DRAW_DIRICHLET_CONDITIONS_H

#include <soft_body_gui_lighting.h>
#include <soft_body_gui_scene_object.h>

#include <gl3_shader_program_manager.h>
#include <gl3_glm.h>

#include <mesh_array.h>
#include <tiny.h>
#include <hyper_engine.h>

#include <cmath> // std::acos

namespace soft_body
{
  namespace gui
  {

    inline void make_dirichlet_geometry(
                                        gl3::ShaderProgramManager & program_manager
                                        , SceneObject & dirichlet_geometry
                                        )
    {
      typedef tiny::MathTypes<float>       MT;
      typedef typename MT::vector3_type    V;
      typedef typename MT::quaternion_type Q;

      gl3::Program & program = program_manager.get(2);

      float radius = 1.00f;

      mesh_array::T3Mesh mesh;
      mesh_array::VertexAttribute<float, mesh_array::T3Mesh> meshX;
      mesh_array::VertexAttribute<float, mesh_array::T3Mesh> meshY;
      mesh_array::VertexAttribute<float, mesh_array::T3Mesh> meshZ;

      mesh_array::make_sphere<MT>(radius, 12, 12, mesh, meshX, meshY, meshZ);

      dirichlet_geometry.m_vbo        = gl3::make_vbo(mesh, meshX, meshY, meshZ);
      dirichlet_geometry.m_solid_vao  = gl3::make_vao(dirichlet_geometry.m_vbo, program, "position", "normal");
    }

    template<typename MT>
    inline void draw_dirichlet_conditions(
                                          hyper::Engine<MT> const & engine
                                          , util::ConfigFile const & params
                                          , glm::mat4 const & projection_matrix
                                          , glm::mat4 const & view_matrix
                                          , std::vector<LightInfo>  & lights
                                          , gl3::ShaderProgramManager  & program_manager
                                          , SceneObject & dirichlet_geometry
                                          )
    {
      typedef typename MT::real_type    T;
      typedef typename MT::vector3_type V;

      using std::acos;
      using std::fabs;

      gl3::check_errors("draw_dirichlet_conditions() invoked");

      gl3::Program & program = program_manager.get(2);

      program.use();

      program.set_uniform("projection_matrix", projection_matrix );
      program.set_uniform("view_matrix", view_matrix );

      set_lighting_uniforms(program, lights);

      typename hyper::Engine<MT>::const_body_iterator body = engine.body_begin();
      typename hyper::Engine<MT>::const_body_iterator end  = engine.body_end();

      for(;body!=end;++body)
      {
        typename std::vector<hyper::DirichletInfo<MT> >::const_iterator c   = body->m_dirichlet_conditions.begin();
        typename std::vector<hyper::DirichletInfo<MT> >::const_iterator end = body->m_dirichlet_conditions.end();

        for( ; c != end; ++ c)
        {
          unsigned int const idx = c->idx();

          T const x = body->m_X[idx];
          T const y = body->m_Y[idx];
          T const z = body->m_Z[idx];
          glm::vec3 const p = glm::vec3( x, y, z );

          float const scale = util::to_value<float>(params.get_value("draw_dirichlet_scale", "1.0"));

          glm::mat4 const scale_matrix = glm::scale(glm::mat4(1.0), glm::vec3(scale,scale,scale));
          glm::mat4 const translation_matrix = glm::translate( glm::mat4(1.0), p);
          glm::mat4 const model_view_matrix =  view_matrix * translation_matrix * scale_matrix;

          program.set_uniform( "model_view_matrix", model_view_matrix);
          
          dirichlet_geometry.m_solid_vao.bind();
          dirichlet_geometry.m_vbo.draw();
          dirichlet_geometry.m_solid_vao.unbind();
        }
      }
      
      program.stop();
    }
    
  }//namespace gui
}//namespace soft_body

// SOFT_BODY_GUI_DRAW_DIRICHLET_CONDITIONS_H
#endif