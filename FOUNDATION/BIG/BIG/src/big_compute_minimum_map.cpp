#include <min_map/big_compute_minimum_map.h>

#include <cmath>

namespace big
{
  
  template < typename T >
  void compute_minimum_map(
                           ublas::vector<T> const & y
                           , ublas::vector<T> const & x
                           , ublas::vector<T> & H
                           )
  {
    using std::min;
    
    size_t const m = y.size();
    
    if(m==0u)
      return;
    
    H.resize(m);
    
    for (size_t i = 0u; i < m; ++ i)
    {
      T const y_i = y(i);
      T const x_i = x(i);
      H(i) = min( x_i , y_i );
    }
  }
  
  template
  void compute_minimum_map<float>(
                                  ublas::vector<float> const & y
                                  , ublas::vector<float> const & x
                                  , ublas::vector<float> & H
                                  );
  
  template
  void compute_minimum_map<double>(
                                   ublas::vector<double> const & y
                                   , ublas::vector<double> const & x
                                   , ublas::vector<double> & H
                                   );
  
} // namespace big
