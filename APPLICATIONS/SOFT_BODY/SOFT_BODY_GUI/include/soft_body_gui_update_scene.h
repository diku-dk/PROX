#ifndef SOFT_BODY_GUI_UPDATE_SCENE_H
#define SOFT_BODY_GUI_UPDATE_SCENE_H

#include <soft_body_gui_scene_object.h>
#include <soft_body_gui_scene_manager.h>

#include <gl3_glm.h>

#include <hyper_engine.h>

namespace soft_body
{
  namespace gui
  {

    template<typename MT>
    inline void update_scene(
                             SceneManager & scene_manager
                             , hyper::Engine<MT> & engine
                             )
    {
      std::vector<SceneObject>::iterator object = scene_manager.m_objects.begin();
      std::vector<SceneObject>::iterator end    = scene_manager.m_objects.end();

      for(;object!=end;++object)
      {
        typename hyper::Engine<MT>::body_type const & body = engine.get_body(object->m_body_idx);

        mesh_array::T3Mesh mesh;
        mesh_array::VertexAttribute<float, mesh_array::T3Mesh> X;
        mesh_array::VertexAttribute<float, mesh_array::T3Mesh> Y;
        mesh_array::VertexAttribute<float, mesh_array::T3Mesh> Z;

        mesh_array::make_t3mesh( body.m_mesh, body.m_X, body.m_Y, body.m_Z, mesh, X, Y, Z );

        gl3::update_vbo(object->m_vbo,mesh, X, Y, Z);
      }
    }

  }//namespace gui
}//namespace soft_body

// SOFT_BODY_GUI_UPDATE_SCENE_H
#endif
