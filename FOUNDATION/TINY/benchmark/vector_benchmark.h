#ifndef BENCHMARK_VECTOR_H
#define BENCHMARK_VECTOR_H

#include <tiny_vector.h>
#include <tiny_vector_functions.h>
#include <tiny_type_traits.h>

#include "measure.h"

namespace vector_benchmark
{
  template<size_t N, typename V>
  inline void benchmark (Measure& measures)
  {    
    size_t const loops = 100000;
    
    V A( V::value_traits::one() );
    V B( V::value_traits::two() );
    
    //Timer<double> timer;
    
    //timer.start();
    for (size_t i = 0; i<loops; ++i)
      inner_prod(A, B);
    //timer.stop();
    measures.SelectOperation(std::string("inner_prod"));
    //measures[N]=timer()/double(loops);
        
    //timer.start();
    for (size_t i = 0; i<loops; ++i)
      norm(B);
    //timer.stop();
    measures.SelectOperation(std::string("norm"));
    //measures[N]=timer()/double(loops);
    
    //timer.start();
    for (size_t i = 0; i<loops; ++i)
      unit(B);
    //timer.stop();
    measures.SelectOperation(std::string("unit"));
    //measures[N]=timer()/double(loops);
  }
  
  template <size_t N>
  struct Unroller
  {
    static void unroll(Measure & measures)
    {
      measures.SelectType(std::string("Float"));
      benchmark<N, tiny::Vector<N,tiny::float_traits> >(measures);
   
      measures.SelectType(std::string("Double"));
      benchmark<N, tiny::Vector<N,tiny::double_traits> >(measures);
      
//      measures.SelectType(std::string("SSE"));
//      benchmark<N, tiny::Vector<N,tiny::sse_float_traits> >(measures);
      
      Unroller<N-1>::unroll(measures);
    }
  };
  
  template<>
  struct Unroller<0>
  {
    static void unroll(Measure& measures)  {}
  };
  
  template <size_t N>
  void run()
  {
    Measure measure;
    Unroller<N>::unroll(measure);
    measure.dump();
  }
  
}//namespace vector_benchmark

// BENCHMARK_VECTOR_H
#endif 

