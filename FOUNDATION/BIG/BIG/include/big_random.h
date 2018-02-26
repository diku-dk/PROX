#ifndef BIG_RANDOM_H
#define BIG_RANDOM_H

#include <big_value_traits.h>

#ifdef BOOST_VERSION  //--- FIXME: Nicer way too see if we have boost?
#  include <boost/random.hpp>
#else
#  include <cstdlib> // for std::srand() and std::rand(), defaults to this if no boost!
#endif

#include <ctime>   // for std::time()

#include <boost/cast.hpp> // Needed for boost::numeric_cast

namespace big
{
  
  
#ifdef BOOST_RANDOM_HPP
  
  
  template<typename value_type>
  class Random
  {
  public:
    
    typedef boost::minstd_rand          generator_type;
    
  protected:
    
    static generator_type &  generator()
    {
      static generator_type tmp(static_cast<unsigned int>(std::time(0)));
      return tmp;
    }
    
  public:
    
    typedef value_type                                                     T;
    typedef boost::uniform_real<T>                                         distribution_type;
    typedef boost::variate_generator<generator_type&, distribution_type >  random_type;
    
    distribution_type       m_distribution;
    random_type             m_random;
    
  public:
    
    Random()
    : m_distribution(big::detail::zero<T>(), big::detail::one<T>())
    , m_random(generator(), m_distribution)
    {}
    
    Random(T lower,T upper)
    : m_distribution(lower,upper)
    , m_random(generator(), m_distribution)
    {}
    
  private:
    
    Random(Random const & rnd){}
    
    Random & operator=(Random const & rnd){return *this;}
    
  public:
    
    T operator()() { return m_random();  }
    
    bool operator==(Random const & rnd) const { return m_distribution == rnd.m_distribution; }
    
  };
  
#else
  
  template <typename value_type>
  class Random
  {
  protected:
    
    typedef value_type  T;
    typedef Random<T>   self;
    
    T m_lower;
    T m_upper;
    
  protected:
    
    static bool & is_initialized()
    {
      static bool initialized = false;
      return initialized;
    }
    
  public:
    
    Random()
    : m_lower(big::detail::zero<T>())
    , m_upper(big::detail::one<T>())
    {
      using std::time;
      if(!is_initialized())
      {
        std::srand(static_cast<unsigned int>(std::time(0)));
        is_initialized() = true;
      }
    }
    
    Random(T lower,T upper)
    : m_lower(lower)
    , m_upper(upper)
    {
      self();
    }
    
  private:
    
    Random(Random const & rnd){}
    Random & operator=(Random const & rnd){return *this;}
    
  public:
    
    T operator()() const
    {
      double rnd = rand()/(1.0*RAND_MAX);
      return boost::numeric_cast<T>(m_lower+(m_upper-m_lower)*rnd);
    }
    
    bool operator==(Random const & rnd) const { return (m_lower==rnd.m_lower && m_upper==rnd.m_upper);  }
    
  };
  
#endif
  
  
} // end namespace big

//BIG_RANDOM_H
#endif
