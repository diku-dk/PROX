#ifndef UTIL_PYTHON_WRITE_PROFILING_H
#define UTIL_PYTHON_WRITE_PROFILING_H

#include <util_profiling.h>
#include <util_python_write_vector.h>

#include <sstream>
#include <string>
#include <vector>
#include <iomanip>

namespace util
{

  /**
   * Convenience function used to generate unique array names for python
   * arrays when outputting all profiling data.
   */
  inline std::string generate_name_py( std::string const & base_name, unsigned int const & idx)
  {
    std::ostringstream name;

    name << base_name << "_" << std::setw(6) << std::setfill('0') << idx;

    name.flush();

    return name.str();
  }

  /**
   * This function will extract all profiling information from util::Profiling
   * and convert the raw data into python arrays and write these into a
   * string using python syntax.
   *
   * This provides end-users with an easy get it all into a file
   * in one-go functionality.
   *
   * @return      Upon return this string holds all the raw profiling data
   *              in python syntax.
   */
  inline std::string python_write_profiling()
  {
    std::stringstream output;

    {
      Profiling::monitors_container const * monitors = Profiling::get_monitors_instance();

      Profiling::monitors_container::const_iterator m   = monitors->begin();
      Profiling::monitors_container::const_iterator end = monitors->end();

      for( ; m != end; ++m)
      {
        std::string              const & name    = m->first;
        Profiling::Monitor       const & monitor = m->second;
        std::vector<float>       const & values  = monitor.get_values();

        output <<  python_write_vector( name, values ) << std::endl;
      }
    }

    {
      Profiling::timer_monitors_container const * monitors = Profiling::get_timer_monitors_instance();

      Profiling::timer_monitors_container::const_iterator m   = monitors->begin();
      Profiling::timer_monitors_container::const_iterator end = monitors->end();

      for( ; m != end; ++m)
      {
        std::string                   const & name    = m->first;
        Profiling::TimerMonitor       const & monitor = m->second;
        std::vector<float>            const & values  = monitor.get_values();

        output <<  python_write_vector( name, values ) << std::endl;
      }

    }

    {
      Profiling::vector_monitors_container const * monitors = Profiling::get_vector_monitors_instance();

      Profiling::vector_monitors_container::const_iterator m   = monitors->begin();
      Profiling::vector_monitors_container::const_iterator end = monitors->end();

      for( ; m != end; ++m)
      {
        std::string                    const & name    = m->first;
        Profiling::VectorMonitor       const & monitor = m->second;
        unsigned int                   const   N       = monitor.get_number_of_entries();

        output << name << " = [" << std::endl;
        for (unsigned int idx = 0u; idx < N; ++idx)
        {
          std::vector<float> const & values = monitor.get_values(idx);
          std::string        const   name_idx = generate_name(name, idx);

            if (idx != N-1)
            {
                output <<  python_write_vector( values ) << ", " << "\n";
            }
            else
            {
                output <<  python_write_vector( values ) << "\n";
            }
        }
        output << "]" << std::endl;

      }

    }

    output.flush();

    return output.str();
  }

}// namespace util

// UTIL_PYTHON_WRITE_PROFILING_H
#endif
