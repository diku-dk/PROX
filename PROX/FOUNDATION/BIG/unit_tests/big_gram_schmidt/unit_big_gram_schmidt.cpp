#include <big_random.h>
#include <big_types.h>
#include <big_generate_random.h>
#include <big_is_orthonormal.h>
#include <big_gram_schmidt.h>

#include <big_generate_PD.h>
#include <big_generate_PSD.h>
#include <big_is_symmetric.h>

#include <big_matlab_write.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(big_modified_gram_schmidt);

BOOST_AUTO_TEST_CASE(random_test_case)
{
  typedef ublas::compressed_matrix<double> matrix_type;

  matrix_type A;

  for(size_t tst=0;tst<5;++tst)
  {
    big::generate_random(10, 10, A);

    bool not_ortho = !big::is_orthonormal( A );
    BOOST_CHECK(not_ortho);

    big::gram_schmidt(A);

    bool did_it = big::is_orthonormal( A );
    BOOST_CHECK(did_it);
  }


  {
    using namespace big;
    big::fast_generate_PD( 10, A );

    bool is_ok = big::is_symmetric( A );
    BOOST_CHECK(is_ok);

    big::generate_PSD( 10, A, 0.5 );

    bool is_also_ok = big::is_symmetric( A );
    BOOST_CHECK(is_also_ok);
  }
}

BOOST_AUTO_TEST_SUITE_END();
