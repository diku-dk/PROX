#include <min_map/big_agglomerate_vector.h>

#include <stdexcept>

namespace big
{
  template<typename T>
  void agglomerate_vector(
                          ublas::vector<T> const & x_a
                          , ublas::vector<T> const & x_b
                          , ublas::vector<size_t> const & new2old
                          , ublas::vector<T> & x
                          )
  {
    size_t const n   = x.size();
    size_t const n_a = x_a.size();
    size_t const n_b = x_b.size();
    
    if ( (n_a + n_b) != n )
      throw std::invalid_argument( "Incompatible dimensions" );
    
    for(size_t i = 0u; i< n_a;++i)
      x( new2old( i ) ) = x_a( i );
    
    for(size_t i = 0u; i< n_b; ++i)
      x( new2old( n_a + i ) ) = x_b( i );
  }
  
  template
  void agglomerate_vector<float>(
                                 ublas::vector<float> const & x_a
                                 , ublas::vector<float> const & x_b
                                 , ublas::vector<size_t> const & new2old
                                 , ublas::vector<float> & x
                                 );
  
  template
  void agglomerate_vector<double>(
                                  ublas::vector<double> const & x_a
                                  , ublas::vector<double> const & x_b
                                  , ublas::vector<size_t> const & new2old
                                  , ublas::vector<double> & x
                                  );
  
  
} // namespace big
