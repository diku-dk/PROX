#ifndef SOFT_BODY_GUI_SCENE_MANAGER_H
#define SOFT_BODY_GUI_SCENE_MANAGER_H

#include <soft_body_gui_scene_object.h>

namespace soft_body
{
  namespace gui
  {

    class SceneManager
    {
    public:

      std::vector<SceneObject> m_objects;

    public:

      void clear()
      {
        m_objects.clear();
      }

      void add(SceneObject const & object)
      {
        m_objects.push_back(object);
      }

    };

  }//namespace gui
}//namespace soft_body

// SOFT_BODY_GUI_SCENE_MANAGER_H
#endif


