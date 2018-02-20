#include <big_prod_trans.h>

#include <cassert>

namespace big
{
  template<typename T>
  void prod_trans(
                         boost::numeric::ublas::compressed_matrix<T> const & A
                         , boost::numeric::ublas::vector<T> const & x
                         , boost::numeric::ublas::vector<T>       & y
                         )
  {
    typedef boost::numeric::ublas::vector<T> vector_type;
    typedef typename vector_type::size_type  size_type;
    
    assert(A.size1()>0            || !"prod_trans(): A was empty"            );
    assert(A.size2()>0            || !"prod_trans(): A was empty"            );
    assert(A.size1() ==  x.size() || !"prod_trans(): incompatible dimensions");
    assert(A.size2() ==  y.size() || !"prod_trans(): incompatible dimensions");
    
    //
    //  Example of compressed matrix format:
    //
    //    0 x y 0
    //    0 0 z 0
    //    w 0 0 0
    //
    //   value_data: [ x y z w ]
    //
    //    Note this is simply all non-zero elements of the sparse matrix stored in a row-by-row wise manner.
    //
    //   index1_data: [ 0 2 3 ]
    //
    //    Note that the number of elements is equal to the number of rows in the
    //    sparse matrix. Each element points to the index in value_data that
    //    corresponds to the starting element of the corresponding row.
    //
    //   index2_data: [ 1 2 2 0 ]
    //
    //    Note this array have the same dimension as value_data. Each element in index2_data
    //    stores the corresponind column index of the matching element in value_data.
    //
    y.clear();
    
    size_type const row_end = A.filled1 () - 1;
    for (size_type i = 0u; i < row_end; ++ i)
    {
      size_type const begin = A.index1_data()[i];
      size_type const end   = A.index1_data()[i + 1];
      for (size_type j = begin; j < end; ++ j)
        y( A.index2_data()[j] ) += A.value_data()[j] * x( i );
    }
  }

  template
  void prod_trans<float>(
                  boost::numeric::ublas::compressed_matrix<float> const & A
                  , boost::numeric::ublas::vector<float> const & x
                  , boost::numeric::ublas::vector<float>       & y
                  );

  template
  void prod_trans<double>(
                  boost::numeric::ublas::compressed_matrix<double> const & A
                  , boost::numeric::ublas::vector<double> const & x
                  , boost::numeric::ublas::vector<double>       & y
                  );

  template<typename T>
  void prod_trans(
                         boost::numeric::ublas::compressed_matrix<T> const & A
                         , boost::numeric::ublas::vector<T> const & x
                         , boost::numeric::ublas::vector<T> const & b
                         , boost::numeric::ublas::vector<T>       & y
                         )
  {
    typedef boost::numeric::ublas::vector<T> vector_type;
    typedef typename vector_type::size_type  size_type;
    
    assert(A.size1()>0            || !"prod_trans(): A was empty"            );
    assert(A.size2()>0            || !"prod_trans(): A was empty"            );
    assert(A.size1() ==  x.size() || !"prod_trans(): incompatible dimensions");
    assert(A.size2() ==  y.size() || !"prod_trans(): incompatible dimensions");
    assert(b.size()  ==  y.size() || !"prod_trans(): incompatible dimensions");
    
    //
    //  Example of compressed matrix format:
    //
    //    0 x y 0
    //    0 0 z 0
    //    w 0 0 0
    //
    //   value_data: [ x y z w ]
    //
    //    Note this is simply all non-zero elements of the sparse matrix stored in a row-by-row wise manner.
    //
    //   index1_data: [ 0 2 3 ]
    //
    //    Note that the number of elements is equal to the number of rows in the
    //    sparse matrix. Each element points to the index in value_data that
    //    corresponds to the starting element of the corresponding row.
    //
    //   index2_data: [ 1 2 2 0 ]
    //
    //    Note this array have the same dimension as value_data. Each element in index2_data
    //    stores the corresponind column index of the matching element in value_data.
    //
    y.assign( b );
    
    size_type const row_end = A.filled1 () - 1;
    for (size_type i = 0u; i < row_end; ++ i)
    {
      size_type const begin = A.index1_data()[i];
      size_type const end   = A.index1_data()[i + 1];
      for (size_type j = begin; j < end; ++ j)
        y( A.index2_data()[j] ) += A.value_data()[j] * x( i );
    }
  }
  
  template
  void prod_trans<float>(
                  boost::numeric::ublas::compressed_matrix<float> const & A
                  , boost::numeric::ublas::vector<float> const & x
                  , boost::numeric::ublas::vector<float> const & b
                  , boost::numeric::ublas::vector<float>       & y
                  );
  
  template
  void prod_trans<double>(
                  boost::numeric::ublas::compressed_matrix<double> const & A
                  , boost::numeric::ublas::vector<double> const & x
                  , boost::numeric::ublas::vector<double> const & b
                  , boost::numeric::ublas::vector<double>       & y
                  );
  
} // end namespace big
