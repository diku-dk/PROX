#include <mesh_array.h>

int main()
{
  using namespace mesh_array;
  
  T3Mesh  M;

  t3_vertex_float_attribute X(M);
  t3_vertex_float_attribute Y(M);
  t3_vertex_float_attribute Z(M);
  
  read_tetgen("../bin/resources/tetgen/tri", M, X, Y, Z);
  
  return 0;
}
