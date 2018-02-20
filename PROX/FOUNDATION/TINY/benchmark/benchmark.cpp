#include "vector_benchmark.h"
#include "matrix_benchmark.h"

int main()
{
  vector_benchmark::run<8>();
  matrix_benchmark::run<6,6>();  
  return 0;
}


