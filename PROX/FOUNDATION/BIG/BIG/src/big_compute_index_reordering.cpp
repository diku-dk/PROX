#include <min_map/big_compute_index_reordering.h>

#include <big_return_codes.h>

#include <cassert>

namespace big
{
  void  compute_index_reordering(
                                 ublas::vector<size_t> const & bitmask
                                 , ublas::vector<size_t>    & old2new
                                 , ublas::vector<size_t>    & new2old
                                 )
  {
    size_t const n = bitmask.size();
    
    old2new.resize(n);
    new2old.resize(n);
    
    size_t r = 0u;
    
    for (size_t i = 0u; i < n; ++ i)
    {
      if(bitmask(i) == IN_ACTIVE)
      {
        old2new( i ) = r;
        new2old( r ) = i;
        ++r;
      }
    }
    
    for (size_t i = 0u; i < n; ++ i)
    {
      if(bitmask(i) == IN_NON_ACTIVE)
      {
        old2new( i ) = r;
        new2old( r ) = i;
        ++r;
      }
    }
    
    assert( r==n || !"compute_index_reordering(): something went wrong");
  }
    
} // namespace big
