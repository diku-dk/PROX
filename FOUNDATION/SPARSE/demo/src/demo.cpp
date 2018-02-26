#include <sparse.h>

#include <iostream>

using namespace sparse;

int main()
{
  std::cout << "Sparse Block Matrix Library Demo" << std::endl;
  
  typedef Block<4,6,float> block_type;  
  block_type b;
  b(0,0) = 0.5f;
  b(1,0) = 1.5f;
  
  typedef Vector<block_type> vector_type;
  vector_type v;
  v.resize(1);
  typename vector_type::iterator iter = v.begin();
  ++iter;
  v.clear();
  
  typedef CompressedRowMatrix<block_type> compressed_type;
  compressed_type C;
  C.resize(12,12,20);
  C.clear();
  
  typedef DiagonalMatrix<block_type> diagonal_type;
  diagonal_type D;
  D.resize(12);
  D.clear();
  
  typedef TwoColumnMatrix<block_type> tow_column_type;
  tow_column_type T;
  T.resize(12,12,24);
  T.clear();
    
  std::cout << "Goodbye..." << std::endl;
	return 0;
}
