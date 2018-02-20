#include <tiny.h>

#include <iostream>

using namespace tiny;

int main()
{   
  std::cout << "Small-sized  Matrix Library Demo" << std::endl;
  
  // Convenience Tool for creating common data types
  typedef tiny::MathTypes<float>  math_types;
  
  // Using a 3D vector
  typedef math_types::vector3_type  V;
  
  V a(1.0f);
  V b = V::make(2.0f, 3.0f, 5.0f );
  V c = a + b * 3.0f;

  // Using a 3-by-3 matrix
  typedef math_types::matrix3x3_type  M;
  
  M A(1.0f);
  M B = M::make_diag( b );
  M C = A*B + B*3.0f;
  
  // Dealing with rotations
  typedef math_types::quaternion_type  Q;
  
  Q q = Q::Ru( 2.0f, V::i() );
  b = tiny::rotate( q, a);
  
  M R = M::Ru( 2.0f, V::i() );
  b = R*a;
  
  // Quaternion/Matrix conversions
  q = tiny::make( R );
  R = tiny::make( q );
  
  // Or use specialized trais for user-specified data types
  tiny::Matrix<6,8,tiny::double_traits>  M6x8d;
  
  std::cout << "Goodbye..." << std::endl;
	return 0;
}
