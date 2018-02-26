#include <tiny.h>
#include <mesh_array.h>
#include <content.h>
#include <simulators.h>

#include <util_string_helper.h>
#include <util_config_file.h>
#include <util_log.h>
#include <util_timestamp.h>
#include <util_file_locator.h>

#include <boost/algorithm/string.hpp>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <string>

namespace rigid_body
{
  namespace xml2mel
  {

    class Application
    {
    public:

      typedef tiny::MathTypes<float> MT;
      typedef MT::real_type          T;
      typedef MT::vector3_type       V;
      typedef MT::matrix3x3_type     M;
      typedef MT::quaternion_type    Q;
      typedef MT::coordsys_type      C;
      typedef MT::value_traits       VT;


    protected:

      std::string   m_xml_scene_file;
      std::string   m_xml_channel_file;

      std::string   m_mel_geometry_file;
      std::string   m_mel_bodies_file;
      std::string   m_mel_keyframes_file;
      std::string   m_mel_run_file;
      std::string   m_maya_binary_file;

      std::string   m_output_path;
      std::string   m_working_directory;

      util::ConfigFile                    m_config_file;
      content::ChannelStorage             m_channel_storage;
      simulators::ProxEngine              m_engine;

    protected:

      std::string make_mel_keyframe(  unsigned int const & body_idx
                                    , float const & time
                                    , float const & x
                                    , float const & y
                                    , float const & z
                                    , float const & qs
                                    , float const & qx
                                    , float const & qy
                                    , float const & qz
                                    )
      {
        std::stringstream stream;

        Q const q = Q(qs,qx,qy,qz);
        M const R = tiny::make(q);

        T xangle = VT::zero();
        T yangle = VT::zero();
        T zangle = VT::zero();

        tiny::XYZ_euler_angles(R,xangle,yangle,zangle);

        xangle = VT::convert_to_degrees( xangle );
        yangle = VT::convert_to_degrees( yangle );
        zangle = VT::convert_to_degrees( zangle );

        stream << "move -a -ws  -xyz "
        << x
        << " "
        << y
        << " "
        << z
        << " body" << body_idx
        << ";"
        << std::endl;


        stream << "rotate -a "
        << xangle
        << " "
        << yangle
        << " "
        << zangle
        << " "
        << " body" << body_idx
        << ";"
        << std::endl;

        stream << "setKeyframe -t "
        << time
        << "sec "
        << " body" << body_idx
        << ";"
        << std::endl;

        return stream.str();
      }

      std::string make_mel_euler_filter(unsigned int const & body_idx)
      {
        std::stringstream stream;

        stream << "filterCurve body"
        << body_idx
        << ".rotateX body"
        << body_idx
        << ".rotateY body"
        << body_idx
        <<".rotateZ;"
        << std::endl;

        return stream.str();
      }

      std::string make_mel_instance(unsigned int const & body_idx, unsigned int const & geometry_idx)
      {
        std::stringstream stream;

        stream << "instance -name body"
        << body_idx
        << "  geometry"
        << geometry_idx
        << ";"
        << std::endl;

        return stream.str();
      }

      std::string make_mel_hide(unsigned int const & geometry_idx)
      {
        std::stringstream stream;

        stream << "hide "
        << "  geometry"
        << geometry_idx
        << ";"
        << std::endl;

        return stream.str();
      }

      std::string make_mel_polygon_geometry(
                                            unsigned int const & geometry_idx
                                            , mesh_array::T3Mesh const & mesh
                                            , mesh_array::VertexAttribute<float, mesh_array::T3Mesh> const & X
                                            , mesh_array::VertexAttribute<float, mesh_array::T3Mesh> const & Y
                                            , mesh_array::VertexAttribute<float, mesh_array::T3Mesh> const & Z
                                            )
      {
        std::stringstream stream;

        unsigned int N = mesh.vertex_size();
        unsigned int T = mesh.triangle_size();

        std::vector<bool>         free_vertex(N, true);
        std::vector<unsigned int> vertex_idx(N, N);

        unsigned int next_vertex_idx = 0u;

        bool first = true;

        for(unsigned int idx = 0u; idx < T; ++idx)
        {
          mesh_array::Triangle const & tri = mesh.triangle(idx);

          unsigned int const & i = tri.i();
          unsigned int const & j = tri.j();
          unsigned int const & k = tri.k();

          float const & xi = X[i];
          float const & yi = Y[i];
          float const & zi = Z[i];

          float const & xj = X[j];
          float const & yj = Y[j];
          float const & zj = Z[j];

          float const & xk = X[k];
          float const & yk = Y[k];
          float const & zk = Z[k];

          if( free_vertex[i] && free_vertex[j]  && free_vertex[k])
          {

            if(first)
            {
              stream << "polyCreateFacet -ch off"
              << " -p "     << xi << " " << yi << " "<< zi
              << " -p "     << xj << " " << yj << " "<< zj
              << " -p "     << xk << " " << yk << " "<< zk
              << " -n geometry" << geometry_idx
              << ";"
              << std::endl;

              first = false;
            }
            else
            {
              stream << "polyAppendVertex -ch off"
              << " -p "     << xi << " " << yi << " "<< zi
              << " -p "     << xj << " " << yj << " "<< zj
              << " -p "     << xk << " " << yk << " "<< zk
              << ";"        << std::endl;
            }

            vertex_idx[i] = next_vertex_idx++;
            vertex_idx[j] = next_vertex_idx++;
            vertex_idx[k] = next_vertex_idx++;

            free_vertex[i] = false;
            free_vertex[j] = false;
            free_vertex[k] = false;
          }
          else if( !free_vertex[i] && !free_vertex[j]  && !free_vertex[k])
          {
            stream << "polyAppendVertex -ch off"
            << " -v "     << vertex_idx[i]
            << " -v "     << vertex_idx[j]
            << " -v "     << vertex_idx[k]
            << ";"        << std::endl;
          }
          else if( !free_vertex[i] && free_vertex[j]  && free_vertex[k])
          {
            stream << "polyAppendVertex -ch off"
            << " -v "     << vertex_idx[i]
            << " -p "     << xj << " " << yj << " "<< zj
            << " -p "     << xk << " " << yk << " "<< zk
            << ";"
            << std::endl;

            vertex_idx[j] = next_vertex_idx++;
            vertex_idx[k] = next_vertex_idx++;

            free_vertex[j] = false;
            free_vertex[k] = false;
          }
          else if( free_vertex[i] && !free_vertex[j]  && free_vertex[k])
          {
            stream << "polyAppendVertex -ch off"
            << " -p "     << xi << " " << yi << " "<< zi
            << " -v "     << vertex_idx[j]
            << " -p "     << xk << " " << yk << " "<< zk
            << ";"
            << std::endl;

            vertex_idx[i] = next_vertex_idx++;
            vertex_idx[k] = next_vertex_idx++;

            free_vertex[i] = false;
            free_vertex[k] = false;
          }
          else if( free_vertex[i] && free_vertex[j]  && !free_vertex[k])
          {
            stream << "polyAppendVertex -ch off"
            << " -p "     << xi << " " << yi << " "<< zi
            << " -p "     << xj << " " << yj << " "<< zj
            << " -v "     << vertex_idx[k]
            << ";"
            << std::endl;

            vertex_idx[i] = next_vertex_idx++;
            vertex_idx[j] = next_vertex_idx++;

            free_vertex[i] = false;
            free_vertex[j] = false;
          }
          else if( !free_vertex[i] && !free_vertex[j]  && free_vertex[k])
          {
            stream << "polyAppendVertex -ch off"
            << " -v "     << vertex_idx[i]
            << " -v "     << vertex_idx[j]
            << " -p "     << xk << " " << yk << " "<< zk
            << ";"
            << std::endl;

            vertex_idx[k] = next_vertex_idx++;

            free_vertex[k] = false;
          }
          else if( free_vertex[i] && !free_vertex[j]  && !free_vertex[k])
          {
            stream << "polyAppendVertex -ch off"
            << " -p "     << xi << " " << yi << " "<< zi
            << " -v "     << vertex_idx[j]
            << " -v "     << vertex_idx[k]
            << ";"
            << std::endl;

            vertex_idx[i] = next_vertex_idx++;

            free_vertex[i] = false;
          }
          else if( !free_vertex[i] && free_vertex[j]  && !free_vertex[k])
          {
            stream << "polyAppendVertex -ch off"
            << " -v "     << vertex_idx[i]
            << " -p "     << xj << " " << yj << " "<< zj
            << " -v "     << vertex_idx[k]
            << ";"
            << std::endl;

            vertex_idx[j] = next_vertex_idx++;

            free_vertex[j] = false;
          }
        }


        stream << std::endl;
        stream << "polySetToFaceNormal geometry" << geometry_idx << ";" << std::endl;
        stream << std::endl;

        return stream.str();
      }

      std::string make_mel_new_scene( std::string const & maya_binary_filename)
      {
        std::stringstream stream;

        stream << "file -f -new;" << std::endl;
        stream << "file -rename \""<< maya_binary_filename << "\";" << std::endl;

        return stream.str();
      }

      std::string make_mel_save_file( )
      {
        std::stringstream stream;

        stream << "file -f -save  -options \"v=0;\" -type \"mayaBinary\";" << std::endl;

        return stream.str();
      }

      std::string make_mel_source( std::string const & mel_filename )
      {
        std::stringstream stream;

        stream << "source  \""<< mel_filename <<"\";" << std::endl;

        return stream.str();
      }

      std::string create_indexed_filename( std::string const & filename, unsigned int const &idx )
      {
        std::stringstream stream;

        std::size_t const found = filename.find_last_of(".");
        std::string const base  = filename.substr(0, found);
        std::string const ext   = filename.substr(found+1);

        stream << base
        <<  "_"
        << std::setw(6)
        << std::setfill('0')
        << idx
        << "."
        << ext;

        return stream.str();
      }


    public:

      int convert(std::string const & cfg_file)
      {
        util::Log logging;

        using std::max;

        m_engine.clear();
        m_channel_storage.clear();
        m_config_file.clear();

        if(cfg_file.empty())
        {

          logging << "Application::convert(): Missing cfg? " << util::Log::newline();

          return 1;
        }

        logging << "Application::convert(): loading cfg: " << cfg_file << util::Log::newline();

        if( ! m_config_file.load(cfg_file) )
        {
          logging << "Application::convert(): could not load cfg-file = " << cfg_file << util::Log::newline();

          return 1;
        }

        util::LogInfo::on()            = util::to_value<bool>(m_config_file.get_value("logging","true"));
        util::LogInfo::console()       = util::to_value<bool>(m_config_file.get_value("console","true"));
        util::LogInfo::filename()      = m_config_file.get_value("log_file","log.txt");

        {
          util::Log logging;

          logging << "### " << util::timestamp() << util::Log::newline();
        }

        m_working_directory      = m_config_file.get_value( "working_directory",       ""                   );
        m_output_path            = m_config_file.get_value( "output_path",             ""                   );

        m_xml_scene_file    = m_config_file.get_value( "xml_scene_file",   "scene.xml"    );
        m_xml_channel_file  = m_config_file.get_value( "xml_channel_file", "channels.xml" );

        m_mel_geometry_file  = m_config_file.get_value( "mel_geometry_file",   "geometry.mel"    );
        m_mel_bodies_file    = m_config_file.get_value( "mel_bodies_file",     "bodies.mel"      );
        m_mel_keyframes_file = m_config_file.get_value( "mel_keyframes_file",  "keyframes.mel"   );
        m_mel_run_file       = m_config_file.get_value( "mel_run_file",        "run.mel"         );
        m_maya_binary_file   = m_config_file.get_value( "maya_bondary_file",   "test.mb"         );

        content::xml_read( util::get_full_filename(m_working_directory, m_xml_scene_file), &m_engine);

        logging << "Application::convert(): read  xml scene file   = " << m_xml_scene_file   << util::Log::newline();

        content::xml_read( util::get_full_filename(m_working_directory, m_xml_channel_file), m_channel_storage);

        logging << "Application::convert(): read  xml channel file = " << m_xml_channel_file << util::Log::newline();

        std::vector<std::string>  mel_filenames;  // Container of all mel scripts that have been generated

        size_t max_keys = 0u;

        for(size_t channel_idx = 0u; channel_idx < m_channel_storage.get_number_of_channels(); ++channel_idx )
        {
          max_keys = max( max_keys,  m_channel_storage.get_number_of_keys(channel_idx) );
        }

        logging << "Application::convert(): Maximum number of keys =  " << max_keys << util::Log::newline();

        //--- Genereate the raw geometry -----------------------------------------
        {
          std::vector<size_t> gids;
          gids.resize(m_engine.get_number_of_geometries());
          m_engine.get_geometry_indices( &gids[0] );

          for(size_t i = 0u; i < gids.size(); ++i)
          {
            size_t const gid = gids[i];

            std::ofstream mel_file;

            std::string const filename = create_indexed_filename(
                                                                 m_output_path + m_mel_geometry_file
                                                                 , gid
                                                                 );

            mel_filenames.push_back(filename);

            mel_file.open(filename.c_str(),std::ios::out);

            if(!mel_file.is_open())
            {
              logging << "Application::convert(): Error could not open file = " << (m_output_path + m_mel_geometry_file) << util::Log::newline();
              return 1;
            }

            mesh_array::T3Mesh mesh;
            mesh_array::VertexAttribute<float, mesh_array::T3Mesh> X;
            mesh_array::VertexAttribute<float, mesh_array::T3Mesh> Y;
            mesh_array::VertexAttribute<float, mesh_array::T3Mesh> Z;

            if( m_engine.get_number_of_tetrameshes(gid) > 0)
            {
              content::add_visual_mesh_of_tetrameshes( gid, &m_engine, mesh, X, Y, Z );
            }
            else
            {
              content::add_visual_mesh_of_boxes(      gid, &m_engine, mesh, X, Y, Z );
              content::add_visual_mesh_of_capsules(   gid, &m_engine, mesh, X, Y, Z );
              content::add_visual_mesh_of_cones(      gid, &m_engine, mesh, X, Y, Z );
              content::add_visual_mesh_of_convexes(   gid, &m_engine, mesh, X, Y, Z );
              content::add_visual_mesh_of_cylinders(  gid, &m_engine, mesh, X, Y, Z );
              content::add_visual_mesh_of_ellipsoids( gid, &m_engine, mesh, X, Y, Z );
              content::add_visual_mesh_of_spheres(    gid, &m_engine, mesh, X, Y, Z );
            }

            std::string const geometry = make_mel_polygon_geometry(gid, mesh, X, Y, Z);

            mel_file << geometry << std::endl;

            mel_file.flush();
            mel_file.close();

          }


        }
        //--- Genereate the bodies -----------------------------------------------
        {
          std::ofstream mel_file;

          std::string const filename = m_output_path + m_mel_bodies_file;

          mel_filenames.push_back(filename);

          mel_file.open(filename.c_str(),std::ios::out);

          if(!mel_file.is_open())
          {
            logging << "Application::convert(): Error could not open file = " << filename << util::Log::newline();
            return 1;
          }

          std::vector<size_t> rids;
          rids.resize(m_engine.get_number_of_rigid_bodies());
          m_engine.get_rigid_body_indices( &rids[0] );

          for(size_t i = 0u; i < rids.size(); ++i)
          {
            size_t const rid = rids[i];
            size_t const gid = m_engine.get_rigid_body_collision_geometry( rid );

            //          size_t      const mid           = m_engine.get_rigid_body_material( rid );
            //          std::string const material_name = m_engine.get_material_name( mid );

            if (m_engine.get_rigid_body_fixed(rid))
            {
            }
            else if (m_engine.has_scripted_motion(rid))
            {
            }
            else
            {
            }

            std::string const body = make_mel_instance(rid, gid);

            mel_file << body << std::endl;
          }

          std::vector<size_t> gids;
          gids.resize(m_engine.get_number_of_geometries());
          m_engine.get_geometry_indices( &gids[0] );

          for(size_t i = 0u; i < gids.size(); ++i)
          {
            size_t const gid = gids[i];

            std::string const hide_geometry = make_mel_hide(gid);

            mel_file << hide_geometry << std::endl;
          }

          mel_file.flush();
          mel_file.close();
        }
        //--- Generate keyframe data ---------------------------------------------
        {
          for(size_t channel_idx = 0u; channel_idx < m_channel_storage.get_number_of_channels(); ++ channel_idx )
          {
            size_t const body_idx = m_channel_storage.get_channel_id( channel_idx);

            std::ofstream mel_file;

            std::string const filename = create_indexed_filename(m_output_path + m_mel_keyframes_file, channel_idx);

            mel_filenames.push_back(filename);

            mel_file.open(filename.c_str(),std::ios::out);

            if(!mel_file.is_open())
            {
              logging << "Application::convert(): Error could not open file = " << filename << util::Log::newline();
              return 1;
            }

            for(unsigned int key_idx = 0u;key_idx< max_keys;++key_idx)
            {
              float x= 0.0f;
              float y= 0.0f;
              float z= 0.0f;

              float qs = 0.0f;
              float qx = 0.0f;
              float qy = 0.0f;
              float qz = 0.0f;

              float time = 0.0f;

              time = m_channel_storage.get_key_time(channel_idx, key_idx);
              m_channel_storage.get_key_position( channel_idx, key_idx, x, y, z);
              m_channel_storage.get_key_orientation( channel_idx, key_idx, qs, qx, qy, qz);

              std::string keyframe = make_mel_keyframe(body_idx, time, x, y, z, qs, qx, qy, qz);
              
              mel_file << keyframe << std::endl;
            }
            
            mel_file << make_mel_euler_filter(body_idx)  << std::endl;
            
            mel_file.flush();
            mel_file.close();
          }
          
        }
        //--- Generate run script ---------------------------------------------
        {
          std::ofstream mel_file;

          std::string const filename = m_output_path + m_mel_run_file;
          
          mel_file.open(filename.c_str(),std::ios::out);
          
          if(!mel_file.is_open())
          {
            logging << "Application::convert(): Error could not open file = " << filename << util::Log::newline();
            return 1;
          }
          
          mel_file << make_mel_new_scene(m_maya_binary_file) << std::endl;
          
          unsigned int number_of_script_files = mel_filenames.size();
          
          for(unsigned int i = 0u; i < number_of_script_files; ++ i)
          {
            mel_file << make_mel_source( mel_filenames[i]) << std::endl;
            
            // Too annoying with student version, got to click on a continue button everytime:-(((
            // mel_file << make_mel_save_file() << std::endl;
          }
          mel_file << make_mel_save_file() << std::endl;
          
          mel_file.flush();
          mel_file.close();
        }
        
        logging << "Application::convert(): Finished converting from xml to mel" << util::Log::newline();
        return 0;
      }
      
    };
    
  }// end namespace xml2mel
}// end namespace rigid_body

int main(int argc, char **argv)
{
  rigid_body::xml2mel::Application app;
  
  std::string cfg_file = "";
  
  if (argv[1])
    cfg_file = argv[1];
  else
    cfg_file = "default.cfg";
  
  return app.convert(cfg_file);
}
