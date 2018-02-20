#include <min_map/big_compute_partitioned_jacobian.h>

#include <big_return_codes.h>

namespace big
{
  
  template<typename T>
  void  compute_partitioned_jacobian(
                                     ublas::compressed_matrix<T> const & A
                                     , ublas::vector<size_t> const & bitmask
                                     , ublas::vector<size_t> const & old2new
                                     , size_t const & cnt_active
                                     , size_t const & cnt_inactive
                                     , ublas::compressed_matrix<T> & A_aa
                                     , ublas::compressed_matrix<T> & A_ab
                                     )
  {
    A_aa.resize( cnt_active, cnt_active, false );
    A_ab.resize( cnt_active, cnt_inactive, false );
    
    size_t const row_end = A.filled1() - 1;
    for (size_t i_old = 0u; i_old < row_end; ++i_old)
    {
      if( bitmask(i_old) == IN_ACTIVE )
      {
        size_t const i_new    = old2new( i_old );
        size_t const begin    = A.index1_data()[i_old];
        size_t const end      = A.index1_data()[i_old + 1];
        
        for (size_t j = begin; j < end; ++ j)
        {
          size_t const j_old = A.index2_data()[j];
          size_t const j_new = old2new( j_old );
          if( bitmask(j_old) == IN_ACTIVE )
            A_aa( i_new, j_new            ) = A.value_data()[j];
          else
            A_ab( i_new, j_new-cnt_active ) = A.value_data()[j];
        }
      }
    }
    
  }
  
  
  template
  void  compute_partitioned_jacobian<float>(
                                            ublas::compressed_matrix<float> const & A
                                            , ublas::vector<size_t> const & bitmask
                                            , ublas::vector<size_t> const & old2new
                                            , size_t const & cnt_active
                                            , size_t const & cnt_inactive
                                            , ublas::compressed_matrix<float> & A_aa
                                            , ublas::compressed_matrix<float> & A_ab
                                            );
  
  template
  void  compute_partitioned_jacobian<double>(
                                             ublas::compressed_matrix<double> const & A
                                             , ublas::vector<size_t> const & bitmask
                                             , ublas::vector<size_t> const & old2new
                                             , size_t const & cnt_active
                                             , size_t const & cnt_inactive
                                             , ublas::compressed_matrix<double> & A_aa
                                             , ublas::compressed_matrix<double> & A_ab
                                             );
  
  
} // namespace big
