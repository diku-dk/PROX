#include <min_map/big_partition_vector.h>

#include <big_return_codes.h>

namespace big
{
  template<typename T >
  void  partition_vector(
                         ublas::vector<T> const & x
                         , ublas::vector<size_t> const & bitmask
                         , ublas::vector<size_t> const & old2new
                         , size_t const & cnt_active
                         , size_t const & cnt_inactive
                         , ublas::vector<T>  & x_a
                         , ublas::vector<T>  & x_b
                         )
  {
    x_a.resize( cnt_active );
    x_b.resize( cnt_inactive );
    
    size_t const n = x.size();
    
    for(size_t i_old = 0u; i_old < n; ++i_old)
    {
      size_t i_new = old2new( i_old );
      if( bitmask(i_old) == IN_ACTIVE )
        x_a( i_new ) = x( i_old);
      else
        x_b( i_new-cnt_active ) = x( i_old);
    }
  }
  
  template
  void  partition_vector<float>(
                                ublas::vector<float> const & x
                                , ublas::vector<size_t> const & bitmask
                                , ublas::vector<size_t> const & old2new
                                , size_t const & cnt_active
                                , size_t const & cnt_inactive
                                , ublas::vector<float>  & x_a
                                , ublas::vector<float>  & x_b
                                );
  template
  void  partition_vector<double>(
                                 ublas::vector<double> const & x
                                 , ublas::vector<size_t> const & bitmask
                                 , ublas::vector<size_t> const & old2new
                                 , size_t const & cnt_active
                                 , size_t const & cnt_inactive
                                 , ublas::vector<double>  & x_a
                                 , ublas::vector<double>  & x_b
                                 );
  
} // namespace big
