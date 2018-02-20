#include <min_map/big_compute_index_sets.h>

#include <big_return_codes.h>

#include <cassert>

namespace big
{
  
  template < typename T>
  void compute_index_sets(
                          ublas::vector<T> const & y
                          , ublas::vector<T> const & x
                          , ublas::vector<size_t> & bitmask
                          , size_t  & cnt_active
                          , size_t  & cnt_inactive
                          )
  {
    size_t const n = x.size();
    
    bitmask.resize(n);
    
    cnt_active   = 0u;
    cnt_inactive = 0u;
    
    for (size_t i = 0u; i < n; ++i)
    {
      if( y(i) >= x(i) )
      {
        bitmask(i) = IN_NON_ACTIVE;
        ++cnt_inactive;
      }
      else
      {
        bitmask(i) = IN_ACTIVE;
        ++cnt_active;
      }
    }
    assert( (cnt_active + cnt_inactive) == n        || !"compute_index_sets(): index sets were inconsistent");
  }
  
  template
  void compute_index_sets<float>(
                                 ublas::vector<float> const & y
                                 , ublas::vector<float> const & x
                                 , ublas::vector<size_t> & bitmask
                                 , size_t  & cnt_active
                                 , size_t  & cnt_inactive
                                 );
  
  template
  void compute_index_sets<double>(
                                  ublas::vector<double> const & y
                                  , ublas::vector<double> const & x
                                  , ublas::vector<size_t> & bitmask
                                  , size_t  & cnt_active
                                  , size_t  & cnt_inactive
                                  );
  
  
} // namespace big
