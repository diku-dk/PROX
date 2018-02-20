#include <big_value_traits.h>

#include <boost/cast.hpp>

namespace big
{
  namespace detail
  {
        
    template <>
    float zero<float>() { return 0.0f; }
    
    template <>
    double zero<double>() { return 0.0; }
    
    template <>
    int zero<int>() { return 0; }
    
    template <>
    long unsigned int zero<long unsigned int>() { return 0; }
    
    template <>
    unsigned int zero<unsigned int>() { return 0u; }
    
    template <>
    float one<float>() { return 1.0f; }
    
    template <>
    double one<double>() { return 1.0; }
    
    template <>
    int one<int>() { return 1; }
    
    template <>
    unsigned int one<unsigned int>() { return 1u; }
        
    template <>
    float two<float>() { return 2.0f; }
    
    template <>
    double two<double>() { return 2.0; }
    
    template <>
    int two<int>() { return 2; }
    
    template <>
    unsigned int two<unsigned int>() { return 2u; }
    
    template <>
    float three<float>() { return 3.0f; }
    
    template <>
    double three<double>() { return 3.0; }
    
    template <>
    int three<int>() { return 3; }
    
    template <>
    unsigned int three<unsigned int>() { return 3u; }
        
    template <>
    float four<float>() { return 4.0f; }
    
    template <>
    double four<double>() { return 4.0; }
    
    template <>
    int four<int>() { return 4; }
    
    template <>
    unsigned int four<unsigned int>() { return 4u; }
    
    template <>
    float eight<float>() { return 8.0f; }
    
    template <>
    double eight<double>() { return 8.0; }
    
    template <>
    int eight<int>() { return 8; }
    
    template <>
    unsigned int eight<unsigned int>() { return 8u; }
    
    template <>
    float half<float>() { return 0.5f; }
    
    template <>
    double half<double>() { return 0.5; }

    template<>
    float pi<float>() { return boost::numeric_cast<float>(M_PI); }
    
    template<>
    double pi<double>() { return boost::numeric_cast<double>(M_PI); }
    
    template<>
    float pi_half<float>() {  return boost::numeric_cast<float>(M_PI_2); }

    template<>
    double pi_half<double>() {  return boost::numeric_cast<double>(M_PI_2); }

    template<>
    float pi_quarter<float>() {  return boost::numeric_cast<float>(M_PI_4); }
    
    template<>
    double pi_quarter<double>() {  return boost::numeric_cast<double>(M_PI_4); }
    
    template <>
    double highest<double>() { return boost::numeric::bounds<double>::highest(); }
        
    
    template <>
    float highest<float>() { return boost::numeric::bounds<float>::highest(); }
    
    template <>
    double lowest<double>() { return boost::numeric::bounds<double>::lowest(); }

    template <>
    float lowest<float>() { return boost::numeric::bounds<float>::lowest(); }
    
    template <>
    double infinity<double>(){ return highest<double>(); }

    template <>
    float infinity<float>(){ return highest<float>(); }

    // one degree in radians
    template<>
    double degree<double>()
    {
      return boost::numeric_cast<double>(0.017453292519943295769236907684886);
    }
    
    template<>
    float degree<float>()
    {
      return boost::numeric_cast<float>(0.017453292519943295769236907684886);
    }
    
    // one radian in degrees
    template<>
    double radian<double>()
    {
      return boost::numeric_cast<double>(57.295779513082320876798154814105);
    }

    template<>
    float radian<float>()
    {
      return boost::numeric_cast<float>(57.295779513082320876798154814105);
    }
    
    
  } // end namespace detail
    
  template<>
  class ValueTraits<float>;
  
  template<>
  class ValueTraits<double>;
  
}  // end namespace big
