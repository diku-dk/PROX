#include <big_project.h>

#include <big_is_number.h>

#include <stdexcept>
#include <cassert>

namespace big
{
  
  template < typename T >
  void project(
               boost::numeric::ublas::vector<T> const & x
               , boost::numeric::ublas::vector<T> & new_x
               )
  {
    size_t const m = x.size();
    
    if(x.size()<=0)
      throw std::invalid_argument("size of x-vector were zero");
    
    if(new_x.size()!=m)
      throw std::invalid_argument("size of new_x-vector were incompatible");
        
    //--- Here we project x_k onto the box defined by l(x_k) and u(x_k)
    for (size_t i = 0u; i < m; ++ i)
    {
      T const x_i = x(i);
      assert( is_number( x_i ) || !"project: x_i was not a number");
      new_x(i) = (x_i < 0.0) ? 0.0 : x_i;
    }
  }
  
  
  template
  void project<float>(
                      boost::numeric::ublas::vector<float> const & x
                      , boost::numeric::ublas::vector<float> & new_x
                      );
  
  template
  void project<double>(
                       boost::numeric::ublas::vector<double> const & x
                       , boost::numeric::ublas::vector<double> & new_x
                       );
  
  
  
} // namespace big
