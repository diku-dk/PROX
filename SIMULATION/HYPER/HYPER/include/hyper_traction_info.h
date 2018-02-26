#ifndef HYPER_TRACTION_INFO_H
#define HYPER_TRACTION_INFO_H

namespace hyper
{

  /**
   * This class encapsulates information about how
   * traction (force per area) is applied to a surface triangle.
   *
   * The surface triangle is specified by the nodal indices
   * corresponding to the triangle vertices and it is assumed
   * that a constant traction (load) is applied over the whole
   * surface of the triangle.
   *
   */
  template<typename MT>
  class TractionInfo
  {
  public:

    typedef typename MT::matrix3x3_type          M;
    typedef typename MT::vector3_type            V;
    typedef typename MT::real_type               T;
    typedef typename MT::value_traits            VT;

  protected:

    unsigned int m_i;
    unsigned int m_j;
    unsigned int m_k;

    V m_traction;

  public:

    unsigned int const & i() const { return this->m_i; }
    unsigned int const & j() const { return this->m_j; }
    unsigned int const & k() const { return this->m_k; }
    unsigned int       & i()       { return this->m_i; }
    unsigned int       & j()       { return this->m_j; }
    unsigned int       & k()       { return this->m_k; }

    V const & traction() const { return this->m_traction; }
    V       & traction()       { return this->m_traction; }

  public:

    TractionInfo()
    : m_i()
    , m_j()
    , m_k()
    , m_traction( V::zero() )
    {}

    virtual ~TractionInfo()
    {}

    TractionInfo(TractionInfo const & info)
    {
      *this = info;
    }

    TractionInfo(
                 unsigned int const & i
                 , unsigned int const & j
                 , unsigned int const & k
                 , V const & traction
                 )
    : m_i(i)
    , m_j(j)
    , m_k(k)
    , m_traction( traction )
    {}

    TractionInfo const & operator=(TractionInfo const & info)
    {
      if( this != &info)
      {
        this->m_i = info.m_i;
        this->m_j = info.m_j;
        this->m_k = info.m_k;
        this->m_traction = info.m_traction;
      }
      return *this;
    }

  };

  template<typename MT>
  inline TractionInfo<MT> make_traction_info(
                                  unsigned int const & i
                                , unsigned int const & j
                                , unsigned int const & k
                                , typename MT::vector3_type const & traction
                                )
  {
    return TractionInfo<MT>(i,j,k,traction);
  }

}// namespace hyper

// HYPER_TRACTION_INFO_H
#endif
