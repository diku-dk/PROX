#include <big_prod_sub.h>

#include <cassert>

namespace big
{
  
  template<typename T>
  void prod_sub(
                boost::numeric::ublas::compressed_matrix<T> const & A
                , boost::numeric::ublas::vector<T> const & x
                , boost::numeric::ublas::vector<T>       & y
                )
  {
    typedef boost::numeric::ublas::vector<T> vector_type;
    typedef typename vector_type::size_type  size_type;
    typedef typename vector_type::value_type real_type;
    
    assert(A.size1()>0            || !"prod_sub(): A was empty"            );
    assert(A.size2()>0            || !"prod_sub(): A was empty"            );
    assert(A.size2() ==  x.size() || !"prod_sub(): incompatible dimensions");
    assert(A.size1() ==  y.size() || !"prod_sub(): incompatible dimensions");
    
    size_type const row_end = A.filled1 () - 1;
    for (size_type i = 0u; i < row_end; ++ i)
    {
      size_type const begin = A.index1_data () [i];
      size_type const end   = A.index1_data () [i + 1];
      real_type t     = real_type();
      for (size_type j = begin; j < end; ++j)
        t += A.value_data()[j] * x(  A.index2_data()[j]   );
      y(i) -= t;
    }
  }
  
  template
  void prod_sub<float>(
                       boost::numeric::ublas::compressed_matrix<float> const & A
                       , boost::numeric::ublas::vector<float> const & x
                       , boost::numeric::ublas::vector<float>       & y
                       );
  
  template
  void prod_sub<double>(
                        boost::numeric::ublas::compressed_matrix<double> const & A
                        , boost::numeric::ublas::vector<double> const & x
                        , boost::numeric::ublas::vector<double>       & y
                        );
  
} // end namespace big
