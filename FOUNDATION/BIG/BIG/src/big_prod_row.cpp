#include <big_prod_row.h>

#include <big_is_number.h>

#include <cassert>

namespace big
{
  
  template<typename T>
  T prod_row(
             boost::numeric::ublas::compressed_matrix<T> const & A
             , boost::numeric::ublas::vector<T> const & x
             , typename boost::numeric::ublas::vector<T>::size_type i
             )
  {
    typedef boost::numeric::ublas::vector<T> vector_type;
    typedef typename vector_type::size_type  size_type;
    typedef typename vector_type::value_type real_type;
    
    assert(A.size1()>0            || !"prod_row(): A was empty"            );
    assert(A.size2()>0            || !"prod_row(): A was empty"            );
    assert(A.size2() ==  x.size() || !"prod_row(): incompatible dimensions");
    assert(i < A.size1()          || !"prod_row(): incompatible dimensions");
    
    real_type value = real_type();
    
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
    
    size_type const row_end = A.filled1 () - 1;
    if(i>=row_end)
      return value;
    
    size_type const begin = A.index1_data () [i];
    size_type const end   = A.index1_data()[i + 1];
    for (size_type j = begin; j < end; ++j)
      value += A.value_data()[j] * x(  A.index2_data()[j]  );
    
    assert(is_number(value) || !"prod_row(): not a number encountered");
    
    return value;
  }
  
  template
  float prod_row<float>(
                        boost::numeric::ublas::compressed_matrix<float> const & A
                        , boost::numeric::ublas::vector<float> const & x
                        , boost::numeric::ublas::vector<float>::size_type i
                        );
  
  template
  double prod_row<double>(
                          boost::numeric::ublas::compressed_matrix<double> const & A
                          , boost::numeric::ublas::vector<double> const & x
                          , boost::numeric::ublas::vector<double>::size_type i
                          );
  
} // end namespace big
