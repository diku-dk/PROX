#include <big_random.h>
#include <big_types.h>
#include <big_moore_penrose_pseudoinverse.h>
#include <big_generate_random.h>


#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>


template<typename matrix_type,typename vector_type>
void test(matrix_type const & A, vector_type  & x, vector_type const & b, vector_type const & y)
{
  using std::min;
  using std::max;

  typedef typename matrix_type::value_type real_type;
  typedef typename matrix_type::size_type  size_type;

  real_type const tol = boost::numeric_cast<real_type>(0.01);

  matrix_type invA;  
  big::svd_moore_penrose_pseudoinverse(A,invA);
  {
    ublas::noalias(x) = ublas::prod( invA, b);
    vector_type tst;
    tst.resize( b.size(), false);
    ublas::noalias(tst) = ublas::prod( A, x);
    for(size_type i = 0; i < b.size(); ++i)
      BOOST_CHECK_CLOSE( real_type( tst(i) ), real_type( b(i) ), tol );
  }

  big::lu_moore_penrose_pseudoinverse(A,invA);
  {
    ublas::noalias(x) = ublas::prod( invA, b);
    vector_type tst;
    tst.resize( b.size(), false);
    ublas::noalias(tst) = ublas::prod( A, x);
    for(size_type i = 0; i < b.size(); ++i)
      BOOST_CHECK_CLOSE( real_type( tst(i) ), real_type( b(i) ), tol );
  }
}

BOOST_AUTO_TEST_SUITE(big_svd);


BOOST_AUTO_TEST_CASE(logic_and_valid_arguments_testing)
{
  typedef ublas::compressed_matrix<double> matrix_type;
  typedef ublas::vector<double>            vector_type;
  typedef vector_type::size_type           size_type;

  {
    matrix_type A;
    matrix_type invA;
    BOOST_CHECK_THROW( big::lu_moore_penrose_pseudoinverse(A, invA), std::invalid_argument );
    BOOST_CHECK_THROW( big::svd_moore_penrose_pseudoinverse(A, invA), std::invalid_argument );
  }
}

BOOST_AUTO_TEST_CASE(random_test_case)
{

  typedef ublas::compressed_matrix<double> matrix_type;
  typedef ublas::vector<double>            vector_type;
  typedef vector_type::size_type           size_type;

  for(size_type tst=0;tst<10;++tst)
  {

    size_t M = 10;
    size_t N = 5;

    matrix_type A;
    vector_type x;
    vector_type b;
    vector_type y;

    // square system
    big::generate_random( M, M, A); // A may not be invertibel
    big::generate_random( M, x);
    b.resize(M,false);
    y.resize(M,false);
    ublas::noalias(b) = ublas::prod(A,x);
    y.assign(x);
    x.clear();
    test(A,x,b,y);

    // more rows than columns
    big::generate_random( M, N, A);  // A must have full column rank
    big::generate_random( N, x);
    b.resize(M,false);
    y.resize(N,false);
    ublas::noalias(b) = ublas::prod(A,x);
    y.assign(x);
    x.clear();
    test(A,x,b,y);

    // more columns than rows
    big::generate_random( N, M, A);  // A must have full row rank
    big::generate_random( M, x);
    b.resize(N,false);
    y.resize(M,false);
    ublas::noalias(b) = ublas::prod(A,x);
    y.assign(x);
    x.clear();
    test(A,x,b,y);
  }

}

BOOST_AUTO_TEST_SUITE_END();
