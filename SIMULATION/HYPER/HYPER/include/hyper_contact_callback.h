#ifndef HYPER_CONTACT_CALLBACK_H
#define HYPER_CONTACT_CALLBACK_H

#include <hyper_body.h>
#include <hyper_contact_point.h>

#include <contacts/geometry_contacts_callback.h>


#include <cassert>
#include <vector>

namespace hyper
{

  namespace detail
  {

    /**
     * Contact Point Generation Callback Functor.
     * This functor acts as the glue between the contact point generation algorithm and the Prox contact point data types and containers.
     *
     * The functor hides the nasty details and data types from the contact point generation library.
     */
    template<typename MT>
    class ContactCallbackFunctor
    : public geometry::ContactsCallback<typename MT::vector3_type>
    {
    public:

      typedef          hyper::Body<MT>            body_type;
      typedef          hyper::ContactPoint<MT>    contact_type;
      typedef typename MT::real_type              T;
      typedef typename MT::vector3_type           V;

    protected:

      body_type                   * m_body_i;   ///< A pointer to body i of the contact.
      body_type                   * m_body_j;   ///< A pointer to body j of the contact.
      std::vector< contact_type > * m_results;  ///< A pointer to a contact point container where all generated contacts should be added to.

    public:

      ContactCallbackFunctor()
      : m_body_i(0)
      , m_body_j(0)
      , m_results(0)
      {}

      ContactCallbackFunctor(body_type * A, body_type * B, std::vector< hyper::ContactPoint<MT> > & results)
      : m_body_i(A)
      , m_body_j(B)
      , m_results(&results)
      {}

      ~ContactCallbackFunctor(){}

      ContactCallbackFunctor& operator=( ContactCallbackFunctor const & callback)
      {
        if (this != &callback)
        {
          this->m_body_i  = callback.m_body_i;
          this->m_body_j  = callback.m_body_j;
          this->m_results = callback.m_results;
        }
        return *this;
      }

      /*
       * Add new contact callback function.
       * The contact point generation library invokes this function call whenever it has discovered a new contact that should be stored.
       *
       * @param p    The position of the new contact point.
       * @param n    The unit normal of the new contact point.
       * @param d    The penetration distance meassure of the contact point.
       */
      void operator()( V const & p, V const & n, T const & d)
      {
        assert( this->m_body_i  || !"operator(): body i is null");
        assert( this->m_body_j  || !"operator(): body j is null");
        assert( this->m_results || !"operator(): results is null");

        contact_type contact;

        contact.set_position( p );
        contact.set_depth( d );
        contact.set_normal( n );
        contact.set_body_i( this->m_body_i );
        contact.set_body_j( this->m_body_j );

        this->m_results->push_back( contact );
      }
    };

  }// namespace detail

} // namespace hyper

// HYPER_CONTACT_CALLBACK_H
#endif
