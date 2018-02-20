#ifndef BENCHMARK_MATRIX_H
#define BENCHMARK_MATRIX_H

#include <tiny_matrix.h>
#include <tiny_matrix_functions.h>
#include <tiny_type_traits.h>

#include "measure.h"

namespace matrix_benchmark
{
  
  template<size_t I, size_t J, typename matrix_type, typename matrix_rot_type, typename vector_type>
  inline void benchmark(Measure& measures)
  {      
    matrix_type a(matrix_type::value_traits::zero()),b(matrix_type::value_traits::zero());
    matrix_rot_type m( matrix_rot_type::value_traits::zero());
    
    vector_type v( vector_type::value_traits::one() );
    
    size_t const loops = 10000;
    
    //Timer<double> timer;
    
    //timer.start();
    for (size_t i = 0; i<loops; ++i)
      a + b;
    //timer.stop();
    measures.SelectOperation(std::string("addition"));
    //measures[J]=timer()/double(loops);
    
    //timer.start();
    for (size_t i = 0; i<loops; ++i)
      a - b;
    //timer.stop();
    measures.SelectOperation(std::string("subtraction"));
    //measures[J]=timer()/double(loops);
            
    //timer.start();
    for (size_t i = 0; i<loops; ++i)
      a * v;
    //timer.stop();
    measures.SelectOperation(std::string("vm-mult"));
    //measures[J]=timer()/double(loops);
    
    //timer.start();
    for (size_t i = 0; i<loops; ++i)
      a * m;
    //timer.stop();
    measures.SelectOperation(std::string("mm-mult"));
    //measures[J]=timer()/double(loops);
  }
  
  template <size_t I,size_t J>
  struct Unroller
  {
    inline static void unroll(Measure& measures)
    {      
      measures.SelectType(std::string("float"));
      benchmark<I,J, tiny::Matrix<I,J,tiny::float_traits>, tiny::Matrix<J, I,tiny::float_traits> , tiny::Vector<J,tiny::float_traits> >(measures);
      
      measures.SelectType(std::string("double"));
      benchmark<I,J, tiny::Matrix<I,J,tiny::double_traits>, tiny::Matrix<J, I,tiny::double_traits> ,tiny::Vector<J,tiny::double_traits> >(measures);
      
//      measures.SelectType(std::string("SSE"));
//      benchmark<I,J, tiny::Matrix<I,J,tiny::sse_float_traits>, tiny::Matrix<J, I, tiny::sse_float_traits> ,tiny::Vector<J,tiny::sse_float_traits> >(measures);
            
      Unroller<I-1, J-1>::unroll(measures);
    }
  };
  
  
  template<>
  struct Unroller<0,0>
  {
    static void unroll(Measure& measures) {}
  };
  
  template <size_t I,size_t J>
  void run()
  {
    Measure measures;
    Unroller<I,J>::unroll(measures);
    measures.dump();
    measures.clear();
  }
}// namespace matrix_benchmark

// BENCHMARK_MATRIX_H
#endif 
