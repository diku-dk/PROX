#ifndef SOFT_BODY_GUI_APPLICATION_H
#define SOFT_BODY_GUI_APPLICATION_H

#include <soft_body_make_procedural_scenes.h>
#include <soft_body_write_meshes.h>
#include <soft_body_gui_lighting.h>
#include <soft_body_gui_scene_manager.h>
#include <soft_body_gui_draw_scene.h>
#include <soft_body_gui_draw_contacts.h>
#include <soft_body_gui_draw_aabbs.h>
#include <soft_body_gui_draw_dirichlet_conditions.h>
#include <soft_body_gui_draw_traction_conditions.h>
#include <soft_body_gui_make_scene.h>
#include <soft_body_gui_update_scene.h>
//#include <soft_body_gui_widget.h>

#include <gl3.h>
#include <gl3_glm.h>
#include <gl3_camera.h>
#include <gl3_frustum.h>
#include <gl3_texture_manager.h>
#include <gl3_shader_program_manager.h>

#include <util_string_helper.h>
#include <util_config_file.h>
#include <util_log.h>
#include <util_timestamp.h>

#include <tiny.h>
#include <geometry.h>
#include <mesh_array.h>
#include <hyper.h>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

#include <cstdio>
#include <iomanip>      // std::setw
#include <sstream>      // std::stringstream
#include <fstream>
#include <string>
#include <cfloat>
#include <vector>
#include <cmath>

typedef hyper::MathPolicy<float> MT;
typedef MT::real_type          T;
typedef MT::vector3_type       V;
typedef MT::matrix3x3_type     M;
typedef MT::quaternion_type    Q;
typedef MT::coordsys_type      C;
typedef MT::value_traits       VT;

namespace soft_body
{
  namespace gui
  {

    class Application
    {
    protected:

      float m_fovy;              ///< Field of view in the y-direction.
      float m_z_near;            ///< Near clipping plane.
      float m_z_far;             ///< Far clipping plane.

      float m_dolly_sensitivity; ///< The dolly mouse sensitivity.
      float	m_pan_sensitivity;   ///< The panning mouse sensitivity.

      bool  m_dolly_mode;        ///< Boolean flag indicating whether the application is currently doing a dolly operation.
      bool	m_pan_mode;          ///< Boolean flag indicating whether the application is currently doing a pan operation.
      bool	m_trackball_mode;    ///< Boolean flag indicating whether the application is currently doing a trackball operation.
      bool	m_selection_mode;    ///< Boolean flag indicating whether the application is currently doing a selection operation.

      float	m_begin_x;           ///< The starting x-pixel coordinate when doing a mouse operation.
      float	m_begin_y;           ///< The starting y-pixel coordinate when doing a mouse operation.

      gl3::Camera	  m_camera;  ///< Camera class taking care of the view transformation (from worlds-space to camera space).
      gl3::Frustum  m_frustum; ///< Frustum class, can be used for optimizing the rendering process by offering simple frustum culling.

      float m_clear_color[4];    ///< RGBA in 0.0 to 1.0 range

      T             m_time;
      T             m_total_time;
      T             m_time_step;
      unsigned int  m_frame_number;

      bool          m_draw_contacts;
      bool          m_draw_wireframe;
      bool          m_draw_aabbs;
      bool          m_draw_dirichlet_conditions;
      bool          m_draw_traction_conditions;

      bool          m_framegrab;
      bool          m_profiling;
      bool          m_did_auto_save;
      bool          m_write_obj_files;

      std::string m_framegrab_file;
      std::string m_matlab_file;

      std::string m_obj_path;
      std::string m_texture_path;
      std::string m_shader_path;
      std::string m_working_directory;
      std::string m_output_path;

      util::ConfigFile                    m_config_file;
      hyper::Engine<MT>                   m_engine;

      gl3::TextureManager                 m_texture_manager;
      gl3::ShaderProgramManager           m_shader_manager;

      SceneManager                        m_scene_manager;
      SceneObject                         m_contact_geometry;
      SceneObject                         m_aabb_geometry;
      SceneObject                         m_dirichlet_geometry;
      SceneObject                         m_traction_geometry;

      std::vector<LightInfo>              m_lights;
      std::vector<MaterialInfo>           m_materials;

      Widget * m_gl_widget;

    protected:

      void framegrab(unsigned int cnt = 0)
      {
        util::Log logging;

        std::stringstream filename;

        // compute output width for filename
        static int width = std::ceil( std::log10( std::ceil( m_total_time / m_time_step ) ) ) + 1;
        filename << m_output_path
                 << m_framegrab_file
                 << std::setw(width)
                 << std::setfill('0')
                 << ++cnt
                 << ".png";

        // create directories if necessary
        boost::filesystem::path framegrab_file( filename.str() );
        if ( framegrab_file.has_parent_path() &&
            ! boost::filesystem::exists( framegrab_file.parent_path() ) &&
            ! boost::filesystem::create_directories( framegrab_file.parent_path() ) )
        {
          logging << "Could not create directories for file '"
          << filename.str()
          << "'"
          << util::Log::newline();
        }

        m_gl_widget->grabFrameBuffer().save( filename.str().c_str() );

        logging << "framegrab(...): Wrote file " << filename.str() << util::Log::newline();
      }

      void reset()
      {
        m_obj_path                 = "";
        m_texture_path             = "";
        m_shader_path              = "";
        m_working_directory        = "";
        m_output_path              = "";

        m_engine.clear();

        m_time      = VT::zero();
        m_time_step = VT::numeric_cast(0.01f);

        m_draw_contacts             = false;
        m_draw_wireframe            = false;
        m_draw_aabbs                = false;
        m_draw_dirichlet_conditions = false;
        m_draw_traction_conditions  = false;
        m_framegrab                 = false;
        m_did_auto_save             = false;
        m_profiling                 = false;
        m_write_obj_files           = false;

        m_framegrab_file           = "screen_";
        m_matlab_file              = "out.m";
        m_frame_number             = 0u;

        m_config_file.clear();

        m_lights.clear();
        m_materials.clear();
        m_scene_manager.clear();
        m_shader_manager.clear();
        m_texture_manager.clear();

        make_default_lighting(m_lights);

        m_dolly_mode             = false;
        m_pan_mode               = false;
        m_trackball_mode         = false;
        m_camera.target_locked() = true;

        {
          util::Log logging;

          logging << "Rotate : hold down left mouse button" << util::Log::newline();
          logging << "Dolly  : hold down middle mouse button" << util::Log::newline();
          logging << "Pan    : press shift and hold down left mouse button" << util::Log::newline();
        }

        glm::vec3 const position = glm::vec3(0,0,20);
        glm::vec3 const target   = glm::vec3(0,0,0);
        glm::vec3 const up       = glm::vec3(0,1,0);

        m_camera.init(position,target,up);

        m_clear_color[0] = 0.2f;
        m_clear_color[1] = 0.2f;
        m_clear_color[2] = 0.2f;
        m_clear_color[3] = 1.0f;
      }

      glm::mat4 compute_view_matrix() const
      {
        glm::mat4 const view_matrix = glm::lookAt(m_camera.position(),m_camera.target(),m_camera.up());

        return view_matrix;
      }

      glm::mat4 compute_projection_matrix() const
      {
        GLfloat viewport[4];
        glGetFloatv(GL_VIEWPORT,viewport);

        float const width  = viewport[2]*1.0f;
        float const height = viewport[3]*1.0f;

        float const aspect = ( 1.0f * width ) / height;

        glm::mat4 const projection_matrix = glm::perspective( glm::radians(m_fovy), aspect, m_z_near, m_z_far );

        return projection_matrix;
      }

    public:

      void load_config_file(std::string const & cfg_file)
      {
        reset();

        m_config_file.load(cfg_file);

        //------------------------------------------------------------------------
        //--- Get Application Params ---------------------------------------------
        //------------------------------------------------------------------------
        util::LogInfo::on()            = util::to_value<bool>(m_config_file.get_value("logging","true"));
        util::LogInfo::console()       = util::to_value<bool>(m_config_file.get_value("console","true"));
        util::LogInfo::filename()      = m_config_file.get_value("log_file","log.txt");

        {
          util::Log logging;

          logging << "### " << util::timestamp() << util::Log::newline();
        }

        m_obj_path               = m_config_file.get_value( "obj_path",                       ""                 );
        m_texture_path           = m_config_file.get_value( "texture_path",                   ""                 );
        m_shader_path            = m_config_file.get_value( "shader_path",                    ""                 );
        m_working_directory      = m_config_file.get_value( "working_directory",              ""                 );
        m_output_path            = m_config_file.get_value( "output_path",                    ""                 );

        m_draw_contacts              = util::to_value<bool>( m_config_file.get_value("draw_contacts", "false"));
        m_draw_wireframe             = util::to_value<bool>( m_config_file.get_value("draw_wireframe", "false"));
        m_draw_aabbs                 = util::to_value<bool>( m_config_file.get_value("draw_aabbs", "false"));
        m_draw_dirichlet_conditions  = util::to_value<bool>( m_config_file.get_value("draw_dirichlet_conditions", "false"));
        m_draw_traction_conditions   = util::to_value<bool>( m_config_file.get_value("draw_traction_conditions", "false"));

        m_total_time             = util::to_value<float>( m_config_file.get_value("total_time", "3.0")       );
        m_framegrab              = util::to_value<bool>( m_config_file.get_value("framegrab", "false")       );
        m_profiling              = util::to_value<bool>( m_config_file.get_value("profiling", "false")       );
        m_framegrab_file         = m_config_file.get_value( "framegrab_file", "screen_" );
        m_matlab_file            = m_config_file.get_value( "matlab_file", "output.m" );
        m_write_obj_files        = util::to_value<bool>(m_config_file.get_value("write_obj_files", "false"));

        //------------------------------------------------------------------------
        //--- Get Engine Params --------------------------------------------------
        //------------------------------------------------------------------------

        m_time_step              = util::to_value<float>( m_config_file.get_value("time_step", "0.01")       );
        m_engine.params().set_time_step(m_time_step);

        bool const use_cfl_condition  = util::to_value<bool>( m_config_file.get_value("use_cfl_condition", "false")  );
        m_engine.params().set_use_cfl_condition(use_cfl_condition);

        std::string time_step_method = m_config_file.get_value("time_step_method", "semi_implicit");
        m_engine.params().set_time_step_method_type(time_step_method);

        float const adaptive_min_dt            = util::to_value<float>( m_config_file.get_value("adaptive_min_dt", "0.001")               );
        float const adaptive_max_dt            = util::to_value<float>( m_config_file.get_value("adaptive_max_dt", "0.01")                );
        float const adaptive_halving_tolerance = util::to_value<float>( m_config_file.get_value("adaptive_halving_tolerance", "0.000001") );
        size_t const adaptive_doubling_count   = util::to_value<size_t>( m_config_file.get_value("adaptive_doubling_count", "5")          );

        m_engine.params().set_adaptive_min_dt( adaptive_min_dt );
        m_engine.params().set_adaptive_max_dt( adaptive_max_dt );
        m_engine.params().set_adaptive_halving_tolerance( adaptive_halving_tolerance );
        m_engine.params().set_adaptive_doubling_count( adaptive_doubling_count );

        float const tetgen_quality_ratio       = util::to_value<float>( m_config_file.get_value("tetgen_quality_ratio", "2.0")       );
        float const tetgen_maximum_volume      = util::to_value<float>( m_config_file.get_value("tetgen_maximum_volume", "0.1")      );
        bool  const tetgen_quiet_output        = util::to_value<bool>(  m_config_file.get_value("tetgen_quiet_output", "true")       );
        bool  const tetgen_suppress_splitting  = util::to_value<bool>(  m_config_file.get_value("tetgen_suppress_splitting", "true") );

        m_engine.params().set_tetgen_quality_ratio( tetgen_quality_ratio );
        m_engine.params().set_tetgen_maximum_volume( tetgen_maximum_volume );
        m_engine.params().set_tetgen_quiet_output( tetgen_quiet_output );
        m_engine.params().set_tetgen_suppress_splitting( tetgen_suppress_splitting );

        float  const envelope          = util::to_value<float>(  m_config_file.get_value("envelope", "0.01") );

        m_engine.params().set_envelope( envelope );

        //------------------------------------------------------------------------
        //---Get Configuration/Scene Params --------------------------------------
        //---Create the configuration in engine ----------------------------------
        //------------------------------------------------------------------------

        std::string const scene_name = m_config_file.get_value("procedural_scene", "default");

        if(scene_name.compare("default") == 0)
        {
          soft_body::make_default_scene(m_engine, m_config_file);
        }
        else if(scene_name.compare("cantilever_tower") == 0)
        {
          soft_body::make_cantilever_tower_scene(m_engine, m_config_file);
        }
        else if(scene_name.compare("plate_stack") == 0)
        {
          soft_body::make_plate_stack_scene(m_engine, m_config_file);
        }
        else
        {
          util::Log logging;

          logging << "ERROR no such scene exist : " << scene_name << util::Log::newline();
        }

        //------------------------------------------------------------------------
        //--- Now we got all simulation content in the world  so now we
        //--- create the corresponding visual content and set visualization params
        //------------------------------------------------------------------------
        make_scene(
                   m_config_file
                   , m_engine
                   , m_materials
                   , m_texture_path
                   , m_texture_manager
                   , m_shader_path
                   , m_shader_manager
                   , m_scene_manager
                   );

        make_contact_geometry(m_shader_manager, m_contact_geometry);
        make_aabb_geometry(m_shader_manager,m_aabb_geometry);
        make_dirichlet_geometry(m_shader_manager,m_dirichlet_geometry);
        make_traction_geometry(m_shader_manager,m_traction_geometry);
        make_config_file_lighting(m_config_file, m_lights);

        if(m_lights.empty())
        {
          make_default_lighting(m_lights);
        }

        float const cam_eye_x = util::to_value<float>( m_config_file.get_value("camera_position_x", "0.0"));
        float const cam_eye_y = util::to_value<float>( m_config_file.get_value("camera_position_y", "4.0"));
        float const cam_eye_z = util::to_value<float>( m_config_file.get_value("camera_position_z", "10.0"));

        float const cam_target_x = util::to_value<float>( m_config_file.get_value("camera_target_x", "0.0"));
        float const cam_target_y = util::to_value<float>( m_config_file.get_value("camera_target_y", "0.0"));
        float const cam_target_z = util::to_value<float>( m_config_file.get_value("camera_target_z", "0.0"));

        float const cam_up_x = util::to_value<float>( m_config_file.get_value("camera_up_x", "0.0"));
        float const cam_up_y = util::to_value<float>( m_config_file.get_value("camera_up_y", "1.0"));
        float const cam_up_z = util::to_value<float>( m_config_file.get_value("camera_up_z", "0.0"));

        glm::vec3 const eye    = glm::vec3(cam_eye_x, cam_eye_y, cam_eye_z);
        glm::vec3 const center = glm::vec3(cam_target_x, cam_target_y, cam_target_z);
        glm::vec3 const up     = glm::vec3(cam_up_x, cam_up_y, cam_up_z);

        m_z_near = util::to_value<float>( m_config_file.get_value("z_near", "0.1"));
        m_z_far  = util::to_value<float>( m_config_file.get_value("z_far", "700.0"));

        m_camera.init(eye,center,up);

        m_clear_color[0] = util::to_value<float>( m_config_file.get_value("clear_color_r", "0.2"));
        m_clear_color[1] = util::to_value<float>( m_config_file.get_value("clear_color_g", "0.2"));
        m_clear_color[2] = util::to_value<float>( m_config_file.get_value("clear_color_b", "0.2"));
        m_clear_color[3] = util::to_value<float>( m_config_file.get_value("clear_color_a", "1.0"));

      }

    public:

      Application()
      : m_fovy(30.0f)
      , m_z_near(0.1f)
      , m_z_far(700.0f)
      , m_dolly_sensitivity(0.25f)
      , m_pan_sensitivity(0.25f)
      , m_dolly_mode(false)
      , m_pan_mode(false)
      , m_trackball_mode(false)
      , m_selection_mode(false)
      , m_begin_x(0.0f)
      , m_begin_y(0.0f)
      , m_camera()
      , m_frustum()
      , m_gl_widget(0)
      {
        m_clear_color[0] = 0.2f;
        m_clear_color[1] = 0.2f;
        m_clear_color[2] = 0.2f;
        m_clear_color[3] = 1.0f;
      }

      virtual ~Application()
      {
        reset();
      }

    public:

      void paint()
      {
        glClearColor(  m_clear_color[0], m_clear_color[1], m_clear_color[2] // rgb
                     , m_clear_color[3]);                                   // a
        gl3::check_errors("paint(): glClearColor");

        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        gl3::check_errors("paint(): glClear");

        if(m_scene_manager.m_objects.empty())  // nothing to draw
          return;

        glm::mat4 const view_matrix       = compute_view_matrix();
        glm::mat4 const projection_matrix = compute_projection_matrix();

        draw_scene(
                                  m_draw_wireframe
                                  , projection_matrix
                                  , view_matrix
                                  , m_lights
                                  , m_materials
                                  , m_texture_manager
                                  , m_shader_manager
                                  , m_scene_manager
                                  );

        if (m_draw_contacts)
        {
          draw_contacts(
                                       m_engine
                                       , m_config_file
                                       , projection_matrix
                                       , view_matrix
                                       , m_lights
                                       , m_shader_manager
                                       , m_contact_geometry
                                       );
        }
        if (m_draw_aabbs)
        {
          draw_aabbs(
                                    projection_matrix
                                    , view_matrix
                                    , m_shader_manager
                                    , m_aabb_geometry
                                    , m_scene_manager
                                    , m_engine
                                    );
        }
        if (m_draw_dirichlet_conditions)
        {
          draw_dirichlet_conditions(
                                                   m_engine
                                                   , m_config_file
                                                   , projection_matrix
                                                   , view_matrix
                                                   , m_lights
                                                   , m_shader_manager
                                                   , m_dirichlet_geometry
                                                   );
        }
        if (m_draw_traction_conditions)
        {
          draw_traction_conditions(
                                                  m_engine
                                                  , m_config_file
                                                  , projection_matrix
                                                  , view_matrix
                                                  , m_lights
                                                  , m_shader_manager
                                                  , m_traction_geometry
                                                  );
        }
      }

      void key_down(unsigned char choice)
      {
        util::Log logging;

        switch(choice)
        {
          case '0':  load_config_file(m_working_directory + "default.cfg");                                 break;
          case '1':  load_config_file(m_working_directory + "cantilever_tower_adaptive_time_stepping.cfg"); break;
          case '2':  load_config_file(m_working_directory + "cantilever_tower_fixed_time_stepping.cfg");    break;
          case '3':  load_config_file(m_working_directory + "plate_stack.cfg");                             break;

          case 'W':

            this->m_draw_wireframe = ! this->m_draw_wireframe;

            if (this->m_draw_wireframe)
              logging << "draw wireframe on" << util::Log::newline();
            else
              logging << "draw wireframe off " << util::Log::newline();
            break;

          case 'C':

            this->m_draw_contacts = ! this->m_draw_contacts;

            if (this->m_draw_contacts)
              logging << "draw contacts on" << util::Log::newline();
            else
              logging << "draw contacts off " << util::Log::newline();
            break;

          case 'D':

            this->m_draw_dirichlet_conditions = ! this->m_draw_dirichlet_conditions;

            if (this->m_draw_dirichlet_conditions)
              logging << "draw dirichlet on" << util::Log::newline();
            else
              logging << "draw dirichlet off " << util::Log::newline();
            break;

          case 'T':

            this->m_draw_traction_conditions = ! this->m_draw_traction_conditions;

            if (this->m_draw_traction_conditions)
              logging << "draw traction on" << util::Log::newline();
            else
              logging << "draw traction off " << util::Log::newline();
            break;

          case 'A':

            this->m_draw_aabbs = ! this->m_draw_aabbs;

            if (this->m_draw_aabbs)
              logging << "draw AABBs on" << util::Log::newline();
            else
              logging << "draw AABBs off " << util::Log::newline();
            break;

          case 'F':

            this->m_framegrab = ! this->m_framegrab;

            if (this->m_framegrab)
              logging << "frame grab on" << util::Log::newline();
            else
              logging << "frame grab off " << util::Log::newline();
            break;

            break;

          case 'M': hyper::write_profiling(m_output_path + m_matlab_file); break;

          case '+': run(); break;

          case 'O':

            this->m_camera.orbit_mode() = ! this->m_camera.orbit_mode();

            if(this->m_camera.orbit_mode())
              logging << "orbit mode on" << util::Log::newline();
            else
              logging << "orbit mode off " << util::Log::newline();
            break;

          case 'I':

            this->m_camera.target_locked() = ! this->m_camera.target_locked();

            if(this->m_camera.target_locked())
              logging << "target is locked" << util::Log::newline();
            else
              logging << "target is free " << util::Log::newline();
            break;
        };
      }

      void init(Widget * widget)
      {
        m_gl_widget = widget;
        reset();
      }

      /**
       * Runs a simulation step.
       *
       * @return false if the simulation is finished, true otherwise.
       */
      bool run()
      {
        if(m_engine.empty())
          return true;

        if (m_time >= m_total_time)
        {
          if(! m_did_auto_save)
          {
            if(m_profiling)
            {
              hyper::write_profiling(m_output_path + m_matlab_file);
            }

            m_did_auto_save = true;
          }

          {
            util::Log logging;

            logging << "Application::run(): Total time reached = "
            << m_total_time
            << " seconds, no more simulation is done... reload to restart"
            << util::Log::newline();
          }

          return false;
        }

        {
          util::Log logging;

          logging << "Application::run(): time = " << m_time << util::Log::newline();
          logging.flush();
        }

        if (m_framegrab)
        {
          framegrab(m_frame_number);
        }

        simulate(m_engine, m_time_step);

        m_time = m_time + m_time_step;

        if (m_write_obj_files)
          soft_body::write_meshes(m_engine, m_output_path, m_frame_number );

        ++m_frame_number;

        update_scene(m_scene_manager, m_engine);
        
        return true;
      }

      void mouse_down(double cur_x,double cur_y,bool shift,bool ctrl,bool alt,bool left,bool middle,bool right)
      {
        if (middle || (alt && left))  // 2008-08-13 micky: not all mice have a "normal" middle button!
          m_dolly_mode = true;
        
        if ( shift && left )
          m_pan_mode = true;
        
        if(!middle && !right && !ctrl && !alt && !shift && left)// only left button allowed
        {
          m_camera.mouse_down( cur_x, cur_y );
          m_trackball_mode = true;
        }
        
        m_begin_x = cur_x;
        m_begin_y = cur_y;
      }
      
      void mouse_up(double cur_x,double cur_y,bool shift,bool ctrl,bool alt,bool left,bool middle,bool right)
      {
        if (m_dolly_mode )
        {
          m_camera.dolly( m_dolly_sensitivity * (cur_y - m_begin_y) );
          m_dolly_mode = false;
        }
        else if ( m_pan_mode )
        {
          m_camera.pan(
                       m_pan_sensitivity*( m_begin_x - cur_x)
                       , m_pan_sensitivity*(cur_y - m_begin_y)
                       );
          m_pan_mode = false;
        }
        else if (m_trackball_mode )
        {
          m_camera.mouse_up( cur_x, cur_y );
          m_trackball_mode = false;
        }
        
        m_begin_x = cur_x;
        m_begin_y = cur_y;
      }
      
      void mouse_move(double cur_x,double cur_y)
      {
        if (m_dolly_mode )
        {
          m_camera.dolly( m_dolly_sensitivity*(cur_y - m_begin_y) );
        }
        else if ( m_pan_mode )
        {
          m_camera.pan(
                       m_pan_sensitivity*(m_begin_x - cur_x)
                       , m_pan_sensitivity*(cur_y - m_begin_y)
                       );
        }
        else if ( m_trackball_mode )
        {
          m_camera.mouse_move( cur_x, cur_y);
        }
        
        m_begin_x = cur_x;
        m_begin_y = cur_y;
      }
      
    };
    
    
  } // end of namespace gui
} // end of namespace soft_body

class Instance
{
public:
  
  static soft_body::gui::Application & app()
  {
    static soft_body::gui::Application my_instance;
    
    return my_instance;
  }

};

// SOFT_BODY_GUI_APPLICATION_H
#endif
