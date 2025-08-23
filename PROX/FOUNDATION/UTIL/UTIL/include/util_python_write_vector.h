#ifndef UTIL_PYTHON_WRITE_VECTOR_H
#define UTIL_PYTHON_WRITE_VECTOR_H

#include <sstream>
#include <string>
#include <vector>
#include <tiny_vector.h>
#include <prox_rigid_body.h>


namespace util
{


  /**
   * This version fo the write_vector is usefull for writing indices into
   * python scripts.
   */
  template<typename T>
  inline std::string python_write_vector(std::vector<T> const & values, size_t const & size)
  {
    std::stringstream output;

    output << "[";


    for(size_t i = 0u; i< size; ++i)
    {
        if (i == size-1)
        {
            output <<  values[i];
        }
        else
        {
            output <<  values[i] << ", ";
        }

    }

    output << "]";

    output.flush();

    return output.str();
  }

  template<typename T>
  inline std::string python_write_vector(std::vector<T> const & values )
  {
    return python_write_vector(values, values.size());
  }

  inline std::string python_write_quaternion(std::vector<tiny::MathTypes<float>::quaternion_type> QUAT)
  {
      std::stringstream output;

      output << "[";

      for(size_t i = 0u; i < QUAT.size(); ++i)
      {
          output << "[" << QUAT[i][0] << ", " << QUAT[i][1] << ", " << QUAT[i][2] << ", " << QUAT[i][3] << "]";
          if (i != QUAT.size()-1)
          {
              output << ", ";
          }
      }
      output << "]";

      output.flush();

      return output.str();
  }

  template<typename T>
  inline std::string python_write_vector(std::string const & name
                                         , std::vector<T> const & values)
  {
    std::stringstream output;

    output << name << " = [";

    for(size_t i = 0u; i< values.size(); ++i)
    {
      output << values[i] << ", ";
    }
    output << "];";

    output.flush();

    return output.str();
  }

}// namespace util

// UTIL_PYTHON_WRITE_VECTOR_H
#endif
