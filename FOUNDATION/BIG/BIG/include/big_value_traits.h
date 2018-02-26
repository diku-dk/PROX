#ifndef BIG_VALUE_TRAITS_H
#define BIG_VALUE_TRAITS_H

namespace big
{
  namespace detail
  {
    template <typename T>
    T zero();
    
    template <>
    float zero<float>();
    
    template <>
    double zero<double>();
    
    template <>
    int zero<int>();
    
    template <>
    long unsigned int zero<long unsigned int>();
    
    template <>
    unsigned int zero<unsigned int>();
    
    template <typename T>
    T one();
    
    template <>
    float one<float>();
    
    template <>
    double one<double>();
    
    template <>
    int one<int>();
    
    template <>
    unsigned int one<unsigned int>();
    
    template <typename T>
    T two();
    
    template <>
    float two<float>();
    
    template <>
    double two<double>();
    
    template <>
    int two<int>();
    
    template <>
    unsigned int two<unsigned int>();
    
    template <typename T>
    T three();
    
    template <>
    float three<float>();
    
    template <>
    double three<double>();
    
    template <>
    int three<int>();
    
    template <>
    unsigned int three<unsigned int>();
    
    template <typename T>
    T four();
    
    template <>
    float four<float>();
    
    template <>
    double four<double>();
    
    template <>
    int four<int>();
    
    template <>
    unsigned int four<unsigned int>();
    
    template <typename T>
    T eight();
    
    template <>
    float eight<float>();
    
    template <>
    double eight<double>();
    
    template <>
    int eight<int>();
    
    template <>
    unsigned int eight<unsigned int>();
    
    template <typename T>
    T half();
    
    template <>
    float half<float>();
    
    template <>
    double half<double>();
    
    template<typename T>
    T pi();
    
    template<typename T>
    T pi_half();
    
    template<typename T>
    T pi_quarter();
    
    template <typename T>
    T highest();
    
    template <typename T>
    T lowest();
    
    template <typename T>
    T infinity();
    
    // one degree in radians
    template<typename T>
    T degree();
    
    // one radian in degrees
    template<typename T>
    T radian();
    
  } // end namespace detail
  
  template <typename T>
  class ValueTraits
  {
  public:
    
    static T zero()      { return detail::zero<T>();      }
    static T one()       { return detail::one<T>();       }
    static T two()       { return detail::two<T>();       }
    static T three()     { return detail::three<T>();     }
    static T four()      { return detail::four<T>();      }
    static T eight()     { return detail::eight<T>();     }
    static T infinity()  { return detail::infinity<T>();  }
    static T half()      { return detail::half<T>();      }
    static T pi()        { return detail::pi<T>();        }
    static T pi_2()      { return detail::pi_half<T>();   }
    static T pi_half()   { return detail::pi_half<T>();   }
    static T pi_quarter(){ return detail::pi_quarter<T>();}
    static T pi_4()      { return detail::pi_quarter<T>();}
    static T degree()    { return detail::degree<T>();    }
    static T radian()    { return detail::radian<T>();    }
    
  };
  
} // end namespace big

//BIG_VALUE_TRAITS_H
#endif
