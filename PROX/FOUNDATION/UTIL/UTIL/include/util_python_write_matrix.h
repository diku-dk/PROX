#ifndef UTIL_PYTHON_WRITE_MATRIX_H
#define UTIL_PYTHON_WRITE_MATRIX_H

#include <util_python_write_vector.h>

#include <sstream>
#include <string>
#include <vector>

namespace util
{

  template<typename T, typename I>
  inline std::string python_write_matrix(
                                         std::vector<I> const & rows
                                         , std::vector<I> const & columns
                                         , std::vector<T> const & values
                                         , std::size_t const & nzeros
                                         , std::size_t const & m
                                         , std::size_t const & n
                                         )
  {
    std::stringstream output;

    output << "[";
    output << python_write_vector( rows, nzeros);
    output << ",";
    output << python_write_vector( columns, nzeros);
    output << ",";
    output << python_write_vector( values, nzeros );
    output << ",";
    output << m;
    output << ",";
    output << n;
    output << ",";
    output << nzeros;
    output << "]";

    output.flush();

    return output.str();
  }

  template<typename T, typename I>
  inline std::string python_write_matrix(
                                           std::string const & name
                                         , std::vector<I> const & rows
                                         , std::vector<I> const & columns
                                         , std::vector<T> const & values
                                         , std::size_t const & nzeros
                                         , std::size_t const & m
                                         , std::size_t const & n
                                         )
  {
    std::stringstream output;

    output << name << " = " << python_write_matrix(rows,columns, values, nzeros, m ,n) << ";";

    output.flush();

    return output.str();
  }

  template<int N, typename MatType>
  inline std::string python_write_matrix_vector(const std::vector<MatType>& mats)
  {
      std::stringstream ss;
      ss << "[";

      for (size_t k = 0; k < mats.size(); ++k)
      {
          const MatType& m = mats[k];
          ss << "["; // start matrix
          for (int r = 0; r < N; ++r)
          {
              ss << "[";
              for (int c = 0; c < N; ++c)
              {
                  std::size_t logical_idx = static_cast<std::size_t>(r * N + c);
                  ss << m.data()[logical_idx];
                  if (c + 1 < N) ss << ", ";
              }
              ss << "]";
              if (r + 1 < N) ss << ", ";
          }
          ss << "]";
          if (k + 1 < mats.size()) ss << ", ";
      }

      ss << "]";
      return ss.str();
  }

}

// UTIL_PYTHON_WRITE_MATRIX_H
#endif
