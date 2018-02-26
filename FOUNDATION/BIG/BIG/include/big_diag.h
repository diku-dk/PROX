#ifndef BIG_DIAG_H
#define BIG_DIAG_H

#include <big_types.h>

namespace big
{
  
  /**
   * Generate Diagonal Matrix.
   * This function is a convenience function that is usefull
   * for quick-and-dirty initialization.
   *
   * @param v          A vector which holds the diagonal values of the resulting matrix.
   * @param A          A square diagonal matrix holding the values of the vector along the diagonal.
   */
  template<typename value_type, typename matrix_type>
  inline void diag( ublas::vector<value_type> const & v, matrix_type & D  )
  {
    size_t const n = v.size();
    
    assert( n > 0u         || !"diag(): n was out of range");
    
    D.resize(n,n,false);
    D.clear();
    
    for(size_t i=0u;i<n;++i)
      D(i,i) = v(i);
  }
  
} // end of namespace big

// BIG_DIAG_H
#endif
