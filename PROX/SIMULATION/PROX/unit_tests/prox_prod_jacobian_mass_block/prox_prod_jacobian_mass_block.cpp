#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/test_tools.hpp>

#include <eigenhelperfunctions.h>

BOOST_AUTO_TEST_SUITE(prod_mass_block);

BOOST_AUTO_TEST_CASE(mass_4x6_block_prod_test_case)
{
    //NOTE: THIS TEST WILL FAIL AS IT IS TRANSLATE FROM A UNUSED TYPE! Just here for code completeyion
    // Use Eigen matrices with row-major storage to match the original indexing
    typedef Eigen::Matrix<float, 4, 6, Eigen::RowMajor> jacobian_block_type;
    typedef Eigen::Matrix<float, 6, 6> mass_block_type;

    jacobian_block_type jb;
    jacobian_block_type jb_result = jacobian_block_type::Zero();

    // Fill jb with values starting from 1.0f, incrementing by 1
    float value = 1.0f;
    for (int i = 0; i < jb.rows(); i++)
    {
        for (int j = 0; j < jb.cols(); j++) { jb(i, j) = value++; }
    }

    mass_block_type mb = mass_block_type::Zero();
    mb(0, 0) = -2.0f;

    // Perform the product operation - assuming it's a matrix multiplication
    jb_result = jb * mb;

    // Check the results - note Eigen uses column-major storage by default
    // but we used RowMajor for jb_result, so we can access elements in row-major order
    BOOST_CHECK_CLOSE(jb_result(0, 0), -2.0f, 1e-5f);
    BOOST_CHECK_CLOSE(jb_result(0, 1), -4.0f, 1e-5f);
    BOOST_CHECK_CLOSE(jb_result(0, 2), -6.0f, 1e-5f);
    BOOST_CHECK_CLOSE(jb_result(0, 3), -8.0f, 1e-5f);
    BOOST_CHECK_CLOSE(jb_result(0, 4), -10.0f, 1e-5f);
    BOOST_CHECK_CLOSE(jb_result(0, 5), -12.0f, 1e-5f);

    BOOST_CHECK_CLOSE(jb_result(1, 0), -14.0f, 1e-5f);
    BOOST_CHECK_CLOSE(jb_result(1, 1), -16.0f, 1e-5f);
    BOOST_CHECK_CLOSE(jb_result(1, 2), -18.0f, 1e-5f);
    BOOST_CHECK_CLOSE(jb_result(1, 3), -20.0f, 1e-5f);
    BOOST_CHECK_CLOSE(jb_result(1, 4), -22.0f, 1e-5f);
    BOOST_CHECK_CLOSE(jb_result(1, 5), -24.0f, 1e-5f);

    BOOST_CHECK_CLOSE(jb_result(2, 0), -26.0f, 1e-5f);
    BOOST_CHECK_CLOSE(jb_result(2, 1), -28.0f, 1e-5f);
    BOOST_CHECK_CLOSE(jb_result(2, 2), -30.0f, 1e-5f);
    BOOST_CHECK_CLOSE(jb_result(2, 3), -32.0f, 1e-5f);
    BOOST_CHECK_CLOSE(jb_result(2, 4), -34.0f, 1e-5f);
    BOOST_CHECK_CLOSE(jb_result(2, 5), -36.0f, 1e-5f);

    BOOST_CHECK_CLOSE(jb_result(3, 0), -38.0f, 1e-5f);
    BOOST_CHECK_CLOSE(jb_result(3, 1), -40.0f, 1e-5f);
    BOOST_CHECK_CLOSE(jb_result(3, 2), -42.0f, 1e-5f);
    BOOST_CHECK_CLOSE(jb_result(3, 3), -44.0f, 1e-5f);
    BOOST_CHECK_CLOSE(jb_result(3, 4), -46.0f, 1e-5f);
    BOOST_CHECK_CLOSE(jb_result(3, 5), -48.0f, 1e-5f);
}

BOOST_AUTO_TEST_SUITE_END();
