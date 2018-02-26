#ifndef BIG_IDENTITY_PRECONDITIONER_H
#define BIG_IDENTITY_PRECONDITIONER_H

namespace big
{
  
  /**
   * An identity preconditioner.
   * This preconditioner does not do anything.
   */
  template<typename T>
  inline void identity(
                       boost::numeric::ublas::compressed_matrix<T> const & /*P*/
                       , boost::numeric::ublas::vector<T>       & e
                       , boost::numeric::ublas::vector<T> const & r
                       )
  {
    ublas::noalias( e ) = r;
  }

} // end of namespace big

// BIG_IDENTITY_PRECONDITIONER_H
#endif
