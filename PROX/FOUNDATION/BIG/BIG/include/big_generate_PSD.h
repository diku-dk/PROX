#ifndef BIG_GENERATE_PSD_H
#define BIG_GENERATE_PSD_H

#include <big_random.h>
#include <big_types.h>
#include <big_generate_random.h>
#include <big_gram_schmidt.h>
#include <big_diag.h>
#include <big_value_traits.h>

namespace big
{
  
  /**
   * Generate Symmetric Positive Semi-Definite (PSD) Matrix.
   * This function is a convenience function that is usefull
   * for quick-and-dirty initialization.
   *
   * Notice, The eigenvalues of the resulting matrix is generated from
   * random values between zero and one.
   *
   * @param n          The number of rows and columns in the resulting matrix.
   * @param fraction   The (stocastic) fraction of zero-valued eigenvalues. If
   *                   this is set to zero then the matrix is a positive definite
   *                   matrix. If the fraction is set close to one then the number
   *                   of zero-eigenvalues grow very large.
   * @param A          Upon return this argument holds an n-by-n PSD matrix.
   */
  template<typename matrix_type>
  inline void generate_PSD(
                           size_t const & n
                           , matrix_type & A
                           , typename matrix_type::value_type const & fraction
                           )
  {
    typedef typename matrix_type::value_type        value_type;
    typedef          big::ValueTraits<value_type>   value_traits;
    typedef          ublas::vector<value_type>      vector_type;
    
    assert( fraction >= value_traits::zero() || !"generate_PSD(): invalid fraction specified");
    assert( fraction <= value_traits::one()  || !"generate_PSD(): invalid fraction specified");
    assert( n > 0                            || !"generate_PSD(): invalid problem size specified");
    
    Random<value_type> value(value_traits::zero(),value_traits::one());
    
    matrix_type Q;
    matrix_type D;
    matrix_type M;
    vector_type d;
    
    big::generate_random( n, d );
    
    if(fraction>value_traits::zero())
    {
      for(size_t i = 0;i< n;++i)
      {
        if(value()<fraction)
          d(i) = value_traits::zero();
      }
    }
    
    big::diag( d, D );
    big::generate_random(n, n, Q);
    big::gram_schmidt(Q);
    
    A.resize(n,n,false);
    M.resize(n,n,false);
    M = ublas::prod(D, ublas::trans(Q) );
    A = ublas::prod( Q, M );
  }
  
} // end of namespace big

// BIG_GENERATE_PSD_H
#endif
