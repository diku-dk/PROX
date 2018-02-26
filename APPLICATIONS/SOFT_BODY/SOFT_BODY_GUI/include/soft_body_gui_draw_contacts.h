#ifndef SOFT_BODY_GUI_DRAW_CONTACTS_H
#define SOFT_BODY_GUI_DRAW_CONTACTS_H

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

    inline void make_contact_geometry(
                                      gl3::ShaderProgramManager & program_manager
                                      , SceneObject & contact_geometry
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

      contact_geometry.m_vbo        = gl3::make_vbo(mesh, X, Y, Z);
      contact_geometry.m_solid_vao  = gl3::make_vao(contact_geometry.m_vbo, program, "position", "normal");
    }

    template<typename MT>
    inline void draw_contacts(
                              hyper::Engine<MT> const & engine
                              , util::ConfigFile const & params
                              , glm::mat4 const & projection_matrix
                              , glm::mat4 const & view_matrix
                              , std::vector<LightInfo>  & lights
                              , gl3::ShaderProgramManager  & program_manager
                              , SceneObject & contact_geometry
                              )
    {
      typedef typename MT::vector3_type V;

      using std::acos;
      using std::fabs;

      gl3::check_errors("draw_contacts() invoked");

      gl3::Program & program = program_manager.get(2);

      program.use();

      program.set_uniform("projection_matrix", projection_matrix );
      program.set_uniform("view_matrix", view_matrix );

      set_lighting_uniforms(program, lights);

      glm::vec3 const up =glm::vec3(0,1,0);

      typename hyper::Engine<MT>::const_contact_iterator c   = engine.contact_begin();
      typename hyper::Engine<MT>::const_contact_iterator end = engine.contact_end();

      for( ; c != end; ++ c)
      {
        V const pp = c->get_position();
        V const nn = c->get_normal();

        glm::vec3 const p = glm::vec3( pp(0), pp(1), pp(2) );
        glm::vec3 const n = glm::normalize( glm::vec3( nn(0), nn(1), nn(2)) );

        float const cos_theta = glm::dot(up, n);

        glm::vec3 const axis = (fabs(cos_theta) < 0.99) ? glm::cross(up, n) : glm::vec3(1,0,0);

        float const radians = acos(cos_theta);

        float const scale = util::to_value<float>(params.get_value("draw_contacts_scaling", "1.0"));

        glm::mat4 const scale_matrix = glm::scale(glm::mat4(1.0), glm::vec3(scale,scale,scale));
        glm::mat4 const rotation_matrix = glm::rotate( glm::mat4(1.0), radians, axis);
        glm::mat4 const translation_matrix = glm::translate( glm::mat4(1.0), p);
        glm::mat4 const model_view_matrix =  view_matrix * translation_matrix * rotation_matrix * scale_matrix;
        
        program.set_uniform( "model_view_matrix", model_view_matrix);
        
        contact_geometry.m_solid_vao.bind();
        contact_geometry.m_vbo.draw();
        contact_geometry.m_solid_vao.unbind();
      }
      
      program.stop();
    }
    
  }//namespace gui
}//namespace soft_body

// SOFT_BODY_GUI_DRAW_CONTACTS_H
#endif