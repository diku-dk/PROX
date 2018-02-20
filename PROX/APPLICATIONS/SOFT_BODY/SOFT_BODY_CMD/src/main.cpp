#include <hyper.h>

#include <soft_body_make_procedural_scenes.h>
#include <soft_body_write_meshes.h>

#include <util_config_file.h>
#include <util_profiling.h>
#include <util_string_helper.h>
#include <util_log.h>
#include <util_timestamp.h>


#include <cassert>
#include <string>

namespace soft_body
{
  namespace cmd
  {

    class Application
    {
    public:

      typedef hyper::MathPolicy<float> MT;
      typedef MT::real_type             T;
      typedef MT::value_traits         VT;

    protected:

      T m_time;
      T m_total_time;
      T m_time_step;

      bool m_profiling;
      bool m_did_auto_save;
      bool m_write_obj_files;

      std::string       m_obj_path;
      std::string       m_output_path;
      std::string       m_working_directory;
      std::string       m_matlab_file;
      util::ConfigFile  m_config_file;
      hyper::Engine<MT> m_engine;
      unsigned int      m_frame_number;

    protected:

      void clear()
      {
        m_obj_path                 = "";
        m_output_path              = "";
        m_working_directory        = "";

        m_engine.clear();

        m_time          = VT::zero();
        m_time_step     = VT::numeric_cast(0.01f);
        m_did_auto_save = false;
        m_profiling     = false;
        m_matlab_file   = "out.m";
        m_write_obj_files = false;

        m_frame_number  = 0u;

        m_config_file.clear();
      }

    public:

      Application() {}

      virtual ~Application()
      {
        clear();
      }

      void load_config_file(std::string const & cfg_file)
      {
        clear();

        {
          util::Log logging;

          logging << "Application::load_config_file(): loading cfg: " << cfg_file << util::Log::newline();
        }

        bool cfg_loaded = m_config_file.load(cfg_file);

        assert(cfg_loaded || !"load_config_file(): cfg fiel could not be loaded");

        util::LogInfo::on()            = util::to_value<bool>(m_config_file.get_value("logging","true"));
        util::LogInfo::console()       = util::to_value<bool>(m_config_file.get_value("console","true"));
        util::LogInfo::filename()      = m_config_file.get_value("log_file","log.txt");

        {
          util::Log logging;

          logging << "### " << util::timestamp() << util::Log::newline();
        }

        m_obj_path               = m_config_file.get_value( "obj_path",                ""                   );
        m_working_directory      = m_config_file.get_value( "working_directory",       ""                   );
        m_output_path            = m_config_file.get_value( "output_path",             ""                   );

        m_total_time       = util::to_value<float>(m_config_file.get_value("total_time", "3.0"));
        m_time_step        = util::to_value<float>(m_config_file.get_value("time_step", "0.01"));
        m_profiling        = util::to_value<bool>(m_config_file.get_value("profiling", "false"));
        m_matlab_file      = m_config_file.get_value("matlab_file", "output.m");
        m_write_obj_files = util::to_value<bool>(m_config_file.get_value("write_obj_files", "false"));

        m_engine.params().set_time_step(m_time_step);

        float const tetgen_quality_ratio = util::to_value<float>(
                                                                 m_config_file.get_value("tetgen_quality_ratio", "2.0")
                                                                 );
        float const tetgen_maximum_volume = util::to_value<float>(
                                                                  m_config_file.get_value("tetgen_maximum_volume", "0.1")
                                                                  );
        bool const tetgen_quiet_output = util::to_value<bool>(
                                                              m_config_file.get_value("tetgen_quiet_output", "true")
                                                              );
        bool const tetgen_suppress_splitting = util::to_value<bool>(
                                                                    m_config_file.get_value("tetgen_suppress_splitting", "true")
                                                                    );

        m_engine.params().set_tetgen_quality_ratio(tetgen_quality_ratio);
        m_engine.params().set_tetgen_maximum_volume(tetgen_maximum_volume);
        m_engine.params().set_tetgen_quiet_output(tetgen_quiet_output);
        m_engine.params().set_tetgen_suppress_splitting(tetgen_suppress_splitting);

        bool const use_open_cl = util::to_value<bool>(
                                                      m_config_file.get_value("use_open_cl", "false")
                                                      );
        size_t const open_cl_platform = util::to_value<size_t>(
                                                               m_config_file.get_value("open_cl_platform", "0")
                                                               );
        size_t const open_cl_device = util::to_value<size_t>(
                                                             m_config_file.get_value("open_cl_device", "0")
                                                             );
        bool const use_gproximity = util::to_value<bool>(
                                                         m_config_file.get_value("use_gproximity", "false")
                                                         );
        bool const use_batching = util::to_value<bool>(
                                                       m_config_file.get_value("use_batching", "true")
                                                       );
        float const envelope = util::to_value<float>(
                                                     m_config_file.get_value("envelope", "0.01")
                                                     );

        m_engine.params().set_use_open_cl(use_open_cl);
        m_engine.params().set_open_cl_device(open_cl_platform);
        m_engine.params().set_open_cl_device(open_cl_device);
        m_engine.params().set_use_gproximity(use_gproximity);
        m_engine.params().set_use_batching(use_batching);
        m_engine.params().set_envelope(envelope);

        std::string const scene_name =
        m_config_file.get_value("procedural_scene", "default");

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

          logging << "Application::load_config_file(): ERROR no such scene exist : " << scene_name << util::Log::newline();
        }
      }

      int run()
      {
        if(m_engine.empty())
          return 0;

        while (m_time < m_total_time)
        {
          {
            util::Log logging;

            logging << "Application::run() time = " << m_time << util::Log::newline();
            logging.flush();
          }

          hyper::simulate(m_engine, m_time_step);

          m_time = m_time + m_time_step;

          if (m_write_obj_files)
            soft_body::write_meshes(m_engine, m_output_path, m_frame_number );

          ++m_frame_number;
        }

        if(! m_did_auto_save)
        {
          if(m_profiling)
          {
            hyper::write_profiling( m_output_path + m_matlab_file );
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

        return 0;
      }

    };

  }// end of namespace cmd

}// end of namespace soft_body

int main(int argc, char** argv)
{
  soft_body::cmd::Application app;

  std::string cfg_file = "";

  if (argv[1])
    cfg_file = argv[1];
  else
    cfg_file = "default.cfg";

  app.load_config_file(cfg_file);

  return app.run();
}
