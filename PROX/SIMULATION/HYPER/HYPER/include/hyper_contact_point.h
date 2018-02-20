#ifndef HYPER_CONTACT_POINT_H
#define HYPER_CONTACT_POINT_H

#include <hyper_body.h>

namespace hyper
{
  
  template< typename MT >
  class ContactPoint
  {
  public:
    
    typedef typename MT::real_type       real_type;
    typedef typename MT::vector3_type    vector3_type;
    typedef          Body<MT>            body_type;
    
  protected:
    
    vector3_type     m_position;
    vector3_type     m_normal;
    real_type        m_depth;
    
    body_type *      m_body_i;
    body_type *      m_body_j;
    
  public:

    ContactPoint()
    : m_position(vector3_type::make(real_type(0),real_type(0),real_type(0)))
    , m_normal(vector3_type::make(real_type(0),real_type(0),real_type(0)))
    , m_depth(real_type(0))
    , m_body_i(0)
    , m_body_j(0)
    {}
    
    virtual ~ContactPoint(){}
    
    ContactPoint (ContactPoint const & point)
    {
      *this = point;
    }
    
    ContactPoint & operator= (ContactPoint const & point)
    {
      if( this != &point)
      {
        this->m_position    = point.m_position;
        this->m_normal      = point.m_normal;
        this->m_depth       = point.m_depth;
        this->m_body_i      = point.m_body_i;
        this->m_body_j      = point.m_body_j;
      }
      return *this;
    }
    
    vector3_type const & get_position() const  {  return this->m_position; }
    vector3_type const & get_normal()   const  {  return this->m_normal;   }
    real_type    const & get_depth()    const  {  return this->m_depth;    }
    body_type    const * get_body_i()    const {  return this->m_body_i;   }
    body_type    const * get_body_j()    const {  return this->m_body_j;   }
    
    void set_position(vector3_type const & p)  {  this->m_position = p;    }
    void set_normal(vector3_type const & n)    {  this->m_normal = n;      }
    void set_depth(real_type const & d)        {  this->m_depth = d;       }
    void set_body_i(body_type const * body_i)  {  this->m_body_i = const_cast<body_type*>(body_i); }  // 2009-11-25 Kenny: hmm can we not get rid of const casts?
    void set_body_j(body_type const * body_j)  {  this->m_body_j = const_cast<body_type*>(body_j); }  // 2009-11-25 Kenny: hmm can we not get rid of const casts?
    
  };

}// namespace hyper

// HYPER_CONTACT_POINT_H
#endif 
