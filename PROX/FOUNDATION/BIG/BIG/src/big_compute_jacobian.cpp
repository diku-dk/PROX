#include <min_map/big_compute_jacobian.h>

#include <big_return_codes.h>

namespace big
{
  
  template<typename T>
  void  compute_jacobian(
                         ublas::compressed_matrix<T> const & A
                         , ublas::vector<size_t> const & bitmask
                         , ublas::compressed_matrix<T>  & J
                         )
  {
    using std::min;
    
    size_t const m = A.size1();
    size_t const n = A.size2();
    
    J.resize(m,n,false);
    
    size_t const row_end = A.filled1() - 1;
    
    for (size_t i = 0u; i < row_end; ++ i)
    {
      if(bitmask(i) == IN_NON_ACTIVE)
      {
        J(i,i) = 1;
      }
      if(bitmask(i) == IN_ACTIVE)
      {
        // J_i* = A_i*
        size_t const begin = A.index1_data()[i];
        size_t const end   = A.index1_data()[i+1];
        for (size_t k = begin; k < end; ++ k)
        {
          size_t const j = A.index2_data()[k];
          J(i,j) = A.value_data()[k];
        }
      }
    }
  }
  
  template
  void  compute_jacobian<float>(
                                ublas::compressed_matrix<float> const & A
                                , ublas::vector<size_t> const & bitmask
                                , ublas::compressed_matrix<float>  & J
                                );
  
  template
  void  compute_jacobian<double>(
                                 ublas::compressed_matrix<double> const & A
                                 , ublas::vector<size_t> const & bitmask
                                 , ublas::compressed_matrix<double>  & J
                                 );
  
  
} // namespace big
