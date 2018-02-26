#ifndef SOFT_BODY_GUI_SCENE_OBJECT_H
#define SOFT_BODY_GUI_SCENE_OBJECT_H

#include <gl3_glm.h>
#include <gl3_vbo.h>
#include <gl3_vao.h>

namespace soft_body
{
  namespace gui
  {

    class SceneObject
    {
    public:

      size_t       m_body_idx;       ///< Body index for looking up simulation data
      size_t       m_visual_idx;     ///< Visual index that identifies texture/material information

      glm::mat4    m_texture_matrix;  ///< texture scaling

      float        m_red;             ///< Object red color intensity
      float        m_green;           ///< Object green color intensity
      float        m_blue;            ///< Object blue color intensity

      gl3::VBO     m_vbo;              ///< Vertex buffer object.
      gl3::VAO     m_solid_vao;        ///< Vertex array object.
      gl3::VAO     m_wire_vao;         ///< Vertex array object.
    };

  }//namesapce gui
}//namesapce soft_body

// SOFT_BODY_GUI_SCENE_OBJECT_H
#endif