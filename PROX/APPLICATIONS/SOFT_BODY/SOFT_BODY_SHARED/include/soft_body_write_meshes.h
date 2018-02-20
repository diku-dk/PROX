#ifndef SOFT_BODY_WRITE_MESHES_H
#define SOFT_BODY_WRITE_MESHES_H

#include <hyper.h>


#include <mesh_array.h>

#include <tiny.h>

#include <util_string_helper.h>
#include <util_config_file.h>
#include <util_log.h>

#include <iomanip>      // std::setw
#include <sstream>      // std::stringstream
#include <string>


namespace soft_body
{

  template<typename MT>
  inline void write_meshes(
                           hyper::Engine<MT> & engine
                           , std::string const & output_path
                           , unsigned int const & frame_number
                           )
  {
    typedef typename MT::real_type                T;
    typedef typename MT::vector3_type             V;
    typedef typename MT::value_traits             VT;
    typedef          hyper::Engine<MT>                   engine_type;
    typedef typename engine_type::body_iterator   body_iterator;

    //util::Log logging;

    for (body_iterator body = engine.body_begin();body != engine.body_end(); ++body)
    {
      std::string const & name = body->get_name();

      std::stringstream filename;
      filename << output_path << std::setw(6) << std::setfill('0') << frame_number << "_" << name << ".obj";

      mesh_array::T3Mesh mesh;
      mesh_array::VertexAttribute<T,mesh_array::T3Mesh> X;
      mesh_array::VertexAttribute<T,mesh_array::T3Mesh> Y;
      mesh_array::VertexAttribute<T,mesh_array::T3Mesh> Z;

      make_t3mesh(
                  body->m_mesh
                  , body->m_X
                  , body->m_Y
                  , body->m_Z
                  , mesh
                  , X
                  , Y
                  , Z
                  );

      write_obj(
                filename.str()
                , mesh
                , X
                , Y
                , Z
                );

    }

  }

} // end of namespace soft_body

// SOFT_BODY_WRITE_MESHES_H
#endif
