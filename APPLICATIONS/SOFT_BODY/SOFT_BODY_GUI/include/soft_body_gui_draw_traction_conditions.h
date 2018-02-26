#ifndef SOFT_BODY_GUI_DRAW_TRACTION_CONDITIONS_H
#define SOFT_BODY_GUI_DRAW_TRACTION_CONDITIONS_H

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

    inline void make_traction_geometry(
                                       gl3::ShaderProgramManager & program_manager
                                       , SceneObject & traction_geometry
                                       )
    {
      typedef tiny::MathTypes<float>       MT;
      typedef typename MT::vector3_type    V;
      typedef typename MT::quaternion_type Q;

      gl3::Program & program = program_manager.get(2);

      float height  = 0.5f;
      float radius1 = 0.25f;
      float radius2 = 0.15f;

      mesh_array::T3Mesh base;
      mesh_array::VertexAttribute<float, mesh_array::T3Mesh> baseX;
      mesh_array::VertexAttribute<float, mesh_array::T3Mesh> baseY;
      mesh_array::VertexAttribute<float, mesh_array::T3Mesh> baseZ;

      mesh_array::T3Mesh shaft;
      mesh_array::VertexAttribute<float, mesh_array::T3Mesh> shaftX;
      mesh_array::VertexAttribute<float, mesh_array::T3Mesh> shaftY;
      mesh_array::VertexAttribute<float, mesh_array::T3Mesh> shaftZ;

      mesh_array::T3Mesh head;
      mesh_array::VertexAttribute<float, mesh_array::T3Mesh> headX;
      mesh_array::VertexAttribute<float, mesh_array::T3Mesh> headY;
      mesh_array::VertexAttribute<float, mesh_array::T3Mesh> headZ;

      mesh_array::make_sphere<MT>(radius1, 12, 12, base, baseX, baseY, baseZ);
      mesh_array::make_cone<MT>(radius1, height, 12, head, headX, headY, headZ);
      mesh_array::make_cylinder<MT>(radius2, height, 12, shaft, shaftX, shaftY, shaftZ);

      mesh_array::translate<MT>(V::make(0, height/2, 0), shaft, shaftX, shaftY, shaftZ);
      mesh_array::translate<MT>(V::make(0, height, 0), head, headX, headY, headZ);

      mesh_array::T3Mesh mesh;
      mesh_array::VertexAttribute<float, mesh_array::T3Mesh> X;
      mesh_array::VertexAttribute<float, mesh_array::T3Mesh> Y;
      mesh_array::VertexAttribute<float, mesh_array::T3Mesh> Z;

      mesh_array::concatenation<MT>( base, baseX, baseY, baseZ, mesh, X, Y, Z );
      mesh_array::concatenation<MT>( shaft, shaftX, shaftY, shaftZ, mesh, X, Y, Z );
      mesh_array::concatenation<MT>( head, headX, headY, headZ, mesh, X, Y, Z );

      traction_geometry.m_vbo        = gl3::make_vbo(mesh, X, Y, Z);
      traction_geometry.m_solid_vao  = gl3::make_vao(traction_geometry.m_vbo, program, "position", "normal");
    }

    template<typename MT>
    inline void draw_traction_conditions(
                                         hyper::Engine<MT> const & engine
                                         , util::ConfigFile const & params
                                         , glm::mat4 const & projection_matrix
                                         , glm::mat4 const & view_matrix
                                         , std::vector<LightInfo>  & lights
                                         , gl3::ShaderProgramManager  & program_manager
                                         , SceneObject & traction_geometry
                                         )
    {
      typedef typename MT::real_type    T;
      typedef typename MT::vector3_type V;
      typedef typename MT::value_traits VT;

      using std::acos;
      using std::fabs;

      gl3::check_errors("draw_traction_conditions() invoked");

      gl3::Program & program = program_manager.get(2);

      program.use();

      program.set_uniform("projection_matrix", projection_matrix );
      program.set_uniform("view_matrix", view_matrix );

      set_lighting_uniforms(program, lights);

      glm::vec3 const up =glm::vec3(0,1,0);

      typename hyper::Engine<MT>::const_body_iterator body   = engine.body_begin();
      typename hyper::Engine<MT>::const_body_iterator b_end  = engine.body_end();

      for(;body!=b_end;++body)
      {
        typename std::vector<hyper::TractionInfo<MT> >::const_iterator c     = body->m_traction_conditions.begin();
        typename std::vector<hyper::TractionInfo<MT> >::const_iterator c_end = body->m_traction_conditions.end();

        for( ; c != c_end; ++ c)
        {
          unsigned int const i = c->i();
          unsigned int const j = c->j();
          unsigned int const k = c->k();

          V const l = c->traction();

          T const x = (body->m_X[i] + body->m_X[j] + body->m_X[k]) / VT::three();
          T const y = (body->m_Y[i] + body->m_Y[j] + body->m_Y[k]) / VT::three();
          T const z = (body->m_Z[i] + body->m_Z[j] + body->m_Z[k]) / VT::three();

          glm::vec3 const p = glm::vec3( x, y, z );
          glm::vec3 const n = glm::vec3( l(0), l(1), l(2) );

          float const cos_theta = glm::dot(up, n);

          glm::vec3 const axis = (fabs(cos_theta) < 0.99) ? glm::cross(up, n) : glm::vec3(1,0,0);

          float const radians = acos(cos_theta);

          float const scale = util::to_value<float>(params.get_value("draw_traction_scale", "1.0"));

          glm::mat4 const scale_matrix       = glm::scale(glm::mat4(1.0), glm::vec3(scale,scale,scale));
          glm::mat4 const translation_matrix = glm::translate( glm::mat4(1.0), p);
          glm::mat4 const rotation_matrix    = glm::rotate( glm::mat4(1.0), radians, axis);
          glm::mat4 const model_view_matrix  =  view_matrix * translation_matrix * rotation_matrix * scale_matrix;
          
          program.set_uniform( "model_view_matrix", model_view_matrix);
          
          traction_geometry.m_solid_vao.bind();
          traction_geometry.m_vbo.draw();
          traction_geometry.m_solid_vao.unbind();
        }
      }
      
      program.stop();
    }
    
  }//namespace gui
}//namespace soft_body

// SOFT_BODY_GUI_DRAW_TRACTION_CONDITIONS_H
#endif
