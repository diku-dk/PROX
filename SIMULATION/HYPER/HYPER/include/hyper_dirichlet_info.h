#ifndef HYPER_DIRICHLET_INFO_H
#define HYPER_DIRICHLET_INFO_H

namespace hyper
{

  template<typename MT>
  class DirichletInfo
  {
  public:

    typedef typename MT::vector3_type V;

  protected:

    unsigned int m_idx;       ///< The nodal index on where to apply the Dirichlet Condition
    V            m_value;     ///< The positional value to enforce, ie. x(idx) = value for all times

  public:

    unsigned int const & idx()   const { return this->m_idx;   }
    V            const & value() const { return this->m_value; }
    unsigned int       & idx()         { return this->m_idx;   }
    V                  & value()       { return this->m_value; }

  public:

    DirichletInfo()
    : m_idx()
    , m_value(V::zero())
    {}

    DirichletInfo(unsigned int const & idx, V const & value)
    : m_idx(idx)
    , m_value(value)
    {}

    virtual ~DirichletInfo(){}

    DirichletInfo(DirichletInfo const & info)
    {
      *this = info;
    }

    DirichletInfo const & operator=(DirichletInfo const & info)
    {
      if( this != &info)
      {
        this->m_idx   = info.m_idx;
        this->m_value = info.m_value;
      }
      return *this;
    }

    bool operator==(DirichletInfo const & info) const { return this->m_idx == info.m_idx; }
    bool operator!=(DirichletInfo const & info) const { return this->m_idx != info.m_idx; }
    bool operator<(DirichletInfo const & info)  const { return this->m_idx < info.m_idx;  }
    bool operator>(DirichletInfo const & info)  const { return this->m_idx > info.m_idx;  }
    bool operator<=(DirichletInfo const & info) const { return this->m_idx <= info.m_idx; }
    bool operator>=(DirichletInfo const & info) const { return this->m_idx >= info.m_idx; }

  };

  template<typename MT>
  inline DirichletInfo<MT> make_dirichlet_info( unsigned int const & idx, typename MT::vector3_type const & value)
  {
    return DirichletInfo<MT>(idx, value);
  }


}// namespace hyper

// HYPER_DIRICHLET_INFO_H
#endif
