#include <eigenhelperall.h>
#include <convex.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(convex_simplex);

BOOST_AUTO_TEST_CASE(simplex_testing)
{
    using T = double;
    typedef convex::Simplex<T> S;

    const EigenVector3<T> not_in_simplex = EigenVector3<T>(1.5, 5.0, 1.0);

    S simplex;

    // Vertify that simplex have been initialized correctly as being empty and ``zeroed''
    BOOST_CHECK(simplex.m_bitmask == 0u);
    for (size_t i = 0u; i < 4; ++i)
    {
        BOOST_CHECK(simplex.m_w[0] == 0.0);
        for (size_t j = 0u; j < 3; ++j)
        {
            BOOST_CHECK(simplex.m_v[i](j) == 0.0);
            BOOST_CHECK(simplex.m_a[i](j) == 0.0);
            BOOST_CHECK(simplex.m_b[i](j) == 0.0);
        }
    }

    // Vertify how different query methods behave on an empty Simplex

    BOOST_CHECK(!convex::is_point_in_simplex(not_in_simplex, simplex));

    BOOST_CHECK(!convex::is_full_simplex(simplex));

    BOOST_CHECK(convex::dimension(simplex) == 0u);

    int bit_A = 0;
    size_t idx_A = 0;
    BOOST_CHECK_THROW(convex::get_used_indices(simplex.m_bitmask, idx_A, bit_A), std::invalid_argument);

    int bit_B = 0;
    size_t idx_B = 0;
    BOOST_CHECK_THROW(convex::get_used_indices(simplex.m_bitmask, idx_A, bit_A, idx_B, bit_B), std::invalid_argument);

    int bit_C = 0;
    size_t idx_C = 0;
    BOOST_CHECK_THROW(convex::get_used_indices(simplex.m_bitmask, idx_A, bit_A, idx_B, bit_B, idx_C, bit_C),
                      std::invalid_argument);

    // Next try to insert one simplex vertex into the simplex

    const EigenVector3<T> p1 = EigenVector3<T>(1.0, 0.0, 0.0);
    const EigenVector3<T> a1 = EigenVector3<T>(1.0, 1.0, 0.0);
    const EigenVector3<T> b1 = EigenVector3<T>(1.0, 0.0, 1.0);

    BOOST_CHECK_NO_THROW(convex::add_point_to_simplex(p1, a1, b1, simplex));

    // Verify how differnt query method works on a 1-simplex
    BOOST_CHECK(!convex::is_full_simplex(simplex));
    BOOST_CHECK(convex::dimension(simplex) == 1u);

    bit_A = 0xFFFF;
    idx_A = 0xFFFF;
    BOOST_CHECK_NO_THROW(convex::get_used_indices(simplex.m_bitmask, idx_A, bit_A));
    BOOST_CHECK(bit_A == 1);
    BOOST_CHECK(idx_A == 0);
    BOOST_CHECK(simplex.m_v[idx_A] == p1);
    BOOST_CHECK(simplex.m_a[idx_A] == a1);
    BOOST_CHECK(simplex.m_b[idx_A] == b1);

    bit_A = 0xFFFF;
    idx_A = 0xFFFF;
    bit_B = 0xFFFF;
    idx_B = 0xFFFF;
    BOOST_CHECK_THROW(convex::get_used_indices(simplex.m_bitmask, idx_A, bit_A, idx_B, bit_B), std::logic_error);

    bit_A = 0xFFFF;
    idx_A = 0xFFFF;
    bit_B = 0xFFFF;
    idx_B = 0xFFFF;
    bit_C = 0xFFFF;
    idx_C = 0xFFFF;
    BOOST_CHECK_THROW(convex::get_used_indices(simplex.m_bitmask, idx_A, bit_A, idx_B, bit_B, idx_C, bit_C),
                      std::logic_error);

    BOOST_CHECK(!convex::is_point_in_simplex(not_in_simplex, simplex));
    BOOST_CHECK(convex::is_point_in_simplex(p1, simplex));

    // Next try to insert one more simplex vertex into the simplex

    const EigenVector3<T> p2 = EigenVector3<T>(2.0, 0.5, 1.0);
    const EigenVector3<T> a2 = EigenVector3<T>(2.0, 1.0, 7.0);
    const EigenVector3<T> b2 = EigenVector3<T>(2.0, 0.5, 1.0);

    BOOST_CHECK_NO_THROW(convex::add_point_to_simplex(p2, a2, b2, simplex));

    BOOST_CHECK(!convex::is_full_simplex(simplex));
    BOOST_CHECK(convex::dimension(simplex) == 2u);

    bit_A = 0xFFFF;
    idx_A = 0xFFFF;
    BOOST_CHECK_NO_THROW(convex::get_used_indices(simplex.m_bitmask, idx_A, bit_A));
    BOOST_CHECK(bit_A == 1);
    BOOST_CHECK(idx_A == 0);
    BOOST_CHECK(simplex.m_v[idx_A] == p1);
    BOOST_CHECK(simplex.m_a[idx_A] == a1);
    BOOST_CHECK(simplex.m_b[idx_A] == b1);

    bit_A = 0xFFFF;
    idx_A = 0xFFFF;
    bit_B = 0xFFFF;
    idx_B = 0xFFFF;
    BOOST_CHECK_NO_THROW(convex::get_used_indices(simplex.m_bitmask, idx_A, bit_A, idx_B, bit_B));

    BOOST_CHECK(bit_A == 1);
    BOOST_CHECK(idx_A == 0);
    BOOST_CHECK(simplex.m_v[idx_A] == p1);
    BOOST_CHECK(simplex.m_a[idx_A] == a1);
    BOOST_CHECK(simplex.m_b[idx_A] == b1);

    BOOST_CHECK(bit_B == 2);
    BOOST_CHECK(idx_B == 1);
    BOOST_CHECK(simplex.m_v[idx_B] == p2);
    BOOST_CHECK(simplex.m_a[idx_B] == a2);
    BOOST_CHECK(simplex.m_b[idx_B] == b2);

    bit_A = 0xFFFF;
    idx_A = 0xFFFF;
    bit_B = 0xFFFF;
    idx_B = 0xFFFF;
    bit_C = 0xFFFF;
    idx_C = 0xFFFF;
    BOOST_CHECK_THROW(convex::get_used_indices(simplex.m_bitmask, idx_A, bit_A, idx_B, bit_B, idx_C, bit_C),
                      std::logic_error);

    BOOST_CHECK(!convex::is_point_in_simplex(not_in_simplex, simplex));
    BOOST_CHECK(convex::is_point_in_simplex(p1, simplex));
    BOOST_CHECK(convex::is_point_in_simplex(p2, simplex));

    // Insert one more simplex vertex

    const EigenVector3<T> p3 = EigenVector3<T>(2.3, 7.5, 1.2);
    const EigenVector3<T> a3 = EigenVector3<T>(2.1, 1.1, 2.3);
    const EigenVector3<T> b3 = EigenVector3<T>(2.2, 2.5, 0.1);

    BOOST_CHECK_NO_THROW(convex::add_point_to_simplex(p3, a3, b3, simplex));

    BOOST_CHECK(!convex::is_full_simplex(simplex));
    BOOST_CHECK(convex::dimension(simplex) == 3u);

    bit_A = 0xFFFF;
    idx_A = 0xFFFF;
    BOOST_CHECK_NO_THROW(convex::get_used_indices(simplex.m_bitmask, idx_A, bit_A));
    BOOST_CHECK(bit_A == 1);
    BOOST_CHECK(idx_A == 0);
    BOOST_CHECK(simplex.m_v[idx_A] == p1);
    BOOST_CHECK(simplex.m_a[idx_A] == a1);
    BOOST_CHECK(simplex.m_b[idx_A] == b1);

    bit_A = 0xFFFF;
    idx_A = 0xFFFF;
    bit_B = 0xFFFF;
    idx_B = 0xFFFF;
    BOOST_CHECK_NO_THROW(convex::get_used_indices(simplex.m_bitmask, idx_A, bit_A, idx_B, bit_B));

    BOOST_CHECK(bit_A == 1);
    BOOST_CHECK(idx_A == 0);
    BOOST_CHECK(simplex.m_v[idx_A] == p1);
    BOOST_CHECK(simplex.m_a[idx_A] == a1);
    BOOST_CHECK(simplex.m_b[idx_A] == b1);

    BOOST_CHECK(bit_B == 2);
    BOOST_CHECK(idx_B == 1);
    BOOST_CHECK(simplex.m_v[idx_B] == p2);
    BOOST_CHECK(simplex.m_a[idx_B] == a2);
    BOOST_CHECK(simplex.m_b[idx_B] == b2);

    bit_A = 0xFFFF;
    idx_A = 0xFFFF;
    bit_B = 0xFFFF;
    idx_B = 0xFFFF;
    bit_C = 0xFFFF;
    idx_C = 0xFFFF;
    BOOST_CHECK_NO_THROW(convex::get_used_indices(simplex.m_bitmask, idx_A, bit_A, idx_B, bit_B, idx_C, bit_C));

    BOOST_CHECK(bit_A == 1);
    BOOST_CHECK(idx_A == 0);
    BOOST_CHECK(simplex.m_v[idx_A] == p1);
    BOOST_CHECK(simplex.m_a[idx_A] == a1);
    BOOST_CHECK(simplex.m_b[idx_A] == b1);

    BOOST_CHECK(bit_B == 2);
    BOOST_CHECK(idx_B == 1);
    BOOST_CHECK(simplex.m_v[idx_B] == p2);
    BOOST_CHECK(simplex.m_a[idx_B] == a2);
    BOOST_CHECK(simplex.m_b[idx_B] == b2);

    BOOST_CHECK(bit_C == 4);
    BOOST_CHECK(idx_C == 2);
    BOOST_CHECK(simplex.m_v[idx_C] == p3);
    BOOST_CHECK(simplex.m_a[idx_C] == a3);
    BOOST_CHECK(simplex.m_b[idx_C] == b3);

    BOOST_CHECK(!convex::is_point_in_simplex(not_in_simplex, simplex));
    BOOST_CHECK(convex::is_point_in_simplex(p1, simplex));
    BOOST_CHECK(convex::is_point_in_simplex(p2, simplex));
    BOOST_CHECK(convex::is_point_in_simplex(p3, simplex));

    // Insert one more vertex then we have a full simplex

    const EigenVector3<T> p4 = EigenVector3<T>(1.3, 1.5, 1.2);
    const EigenVector3<T> a4 = EigenVector3<T>(1.1, 1.1, 1.3);
    const EigenVector3<T> b4 = EigenVector3<T>(1.2, 1.5, 1.1);

    BOOST_CHECK_NO_THROW(convex::add_point_to_simplex(p4, a4, b4, simplex));

    BOOST_CHECK(convex::is_full_simplex(simplex));
    BOOST_CHECK(convex::dimension(simplex) == 4u);

    bit_A = 0xFFFF;
    idx_A = 0xFFFF;
    BOOST_CHECK_NO_THROW(convex::get_used_indices(simplex.m_bitmask, idx_A, bit_A));
    BOOST_CHECK(bit_A == 1);
    BOOST_CHECK(idx_A == 0);
    BOOST_CHECK(simplex.m_v[idx_A] == p1);
    BOOST_CHECK(simplex.m_a[idx_A] == a1);
    BOOST_CHECK(simplex.m_b[idx_A] == b1);

    bit_A = 0xFFFF;
    idx_A = 0xFFFF;
    bit_B = 0xFFFF;
    idx_B = 0xFFFF;
    BOOST_CHECK_NO_THROW(convex::get_used_indices(simplex.m_bitmask, idx_A, bit_A, idx_B, bit_B));

    BOOST_CHECK(bit_A == 1);
    BOOST_CHECK(idx_A == 0);
    BOOST_CHECK(simplex.m_v[idx_A] == p1);
    BOOST_CHECK(simplex.m_a[idx_A] == a1);
    BOOST_CHECK(simplex.m_b[idx_A] == b1);

    BOOST_CHECK(bit_B == 2);
    BOOST_CHECK(idx_B == 1);
    BOOST_CHECK(simplex.m_v[idx_B] == p2);
    BOOST_CHECK(simplex.m_a[idx_B] == a2);
    BOOST_CHECK(simplex.m_b[idx_B] == b2);

    bit_A = 0xFFFF;
    idx_A = 0xFFFF;
    bit_B = 0xFFFF;
    idx_B = 0xFFFF;
    bit_C = 0xFFFF;
    idx_C = 0xFFFF;
    BOOST_CHECK_NO_THROW(convex::get_used_indices(simplex.m_bitmask, idx_A, bit_A, idx_B, bit_B, idx_C, bit_C));

    BOOST_CHECK(bit_A == 1);
    BOOST_CHECK(idx_A == 0);
    BOOST_CHECK(simplex.m_v[idx_A] == p1);
    BOOST_CHECK(simplex.m_a[idx_A] == a1);
    BOOST_CHECK(simplex.m_b[idx_A] == b1);

    BOOST_CHECK(bit_B == 2);
    BOOST_CHECK(idx_B == 1);
    BOOST_CHECK(simplex.m_v[idx_B] == p2);
    BOOST_CHECK(simplex.m_a[idx_B] == a2);
    BOOST_CHECK(simplex.m_b[idx_B] == b2);

    BOOST_CHECK(bit_C == 4);
    BOOST_CHECK(idx_C == 2);
    BOOST_CHECK(simplex.m_v[idx_C] == p3);
    BOOST_CHECK(simplex.m_a[idx_C] == a3);
    BOOST_CHECK(simplex.m_b[idx_C] == b3);

    BOOST_CHECK(simplex.m_v[3] == p4);
    BOOST_CHECK(simplex.m_a[3] == a4);
    BOOST_CHECK(simplex.m_b[3] == b4);

    BOOST_CHECK(!convex::is_point_in_simplex(not_in_simplex, simplex));
    BOOST_CHECK(convex::is_point_in_simplex(p1, simplex));
    BOOST_CHECK(convex::is_point_in_simplex(p2, simplex));
    BOOST_CHECK(convex::is_point_in_simplex(p3, simplex));
    BOOST_CHECK(convex::is_point_in_simplex(p4, simplex));

    // Test what happens if we try to insert five vertices into the simplex

    const EigenVector3<T> p5 = EigenVector3<T>(2.3, 2.5, 2.2);
    const EigenVector3<T> a5 = EigenVector3<T>(2.1, 2.1, 2.3);
    const EigenVector3<T> b5 = EigenVector3<T>(2.2, 2.5, 2.1);

    BOOST_CHECK_THROW(convex::add_point_to_simplex(p5, a5, b5, simplex), std::logic_error);

    // Now let us erase one of the simplex vertices

    simplex.m_bitmask = simplex.m_bitmask & ~bit_B;

    BOOST_CHECK(!convex::is_full_simplex(simplex));
    BOOST_CHECK(convex::dimension(simplex) == 3u);

    bit_A = 0xFFFF;
    idx_A = 0xFFFF;
    BOOST_CHECK_NO_THROW(convex::get_used_indices(simplex.m_bitmask, idx_A, bit_A));
    BOOST_CHECK(bit_A == 1);
    BOOST_CHECK(idx_A == 0);
    BOOST_CHECK(simplex.m_v[idx_A] == p1);
    BOOST_CHECK(simplex.m_a[idx_A] == a1);
    BOOST_CHECK(simplex.m_b[idx_A] == b1);

    bit_A = 0xFFFF;
    idx_A = 0xFFFF;
    bit_B = 0xFFFF;
    idx_B = 0xFFFF;
    BOOST_CHECK_NO_THROW(convex::get_used_indices(simplex.m_bitmask, idx_A, bit_A, idx_B, bit_B));

    BOOST_CHECK(bit_A == 1);
    BOOST_CHECK(idx_A == 0);
    BOOST_CHECK(simplex.m_v[idx_A] == p1);
    BOOST_CHECK(simplex.m_a[idx_A] == a1);
    BOOST_CHECK(simplex.m_b[idx_A] == b1);

    BOOST_CHECK(bit_B == 4);
    BOOST_CHECK(idx_B == 2);
    BOOST_CHECK(simplex.m_v[idx_B] == p3);
    BOOST_CHECK(simplex.m_a[idx_B] == a3);
    BOOST_CHECK(simplex.m_b[idx_B] == b3);

    bit_A = 0xFFFF;
    idx_A = 0xFFFF;
    bit_B = 0xFFFF;
    idx_B = 0xFFFF;
    bit_C = 0xFFFF;
    idx_C = 0xFFFF;
    BOOST_CHECK_NO_THROW(convex::get_used_indices(simplex.m_bitmask, idx_A, bit_A, idx_B, bit_B, idx_C, bit_C));

    BOOST_CHECK(bit_A == 1);
    BOOST_CHECK(idx_A == 0);
    BOOST_CHECK(simplex.m_v[idx_A] == p1);
    BOOST_CHECK(simplex.m_a[idx_A] == a1);
    BOOST_CHECK(simplex.m_b[idx_A] == b1);

    BOOST_CHECK(bit_B == 4);
    BOOST_CHECK(idx_B == 2);
    BOOST_CHECK(simplex.m_v[idx_B] == p3);
    BOOST_CHECK(simplex.m_a[idx_B] == a3);
    BOOST_CHECK(simplex.m_b[idx_B] == b3);

    BOOST_CHECK(bit_C == 8);
    BOOST_CHECK(idx_C == 3);
    BOOST_CHECK(simplex.m_v[idx_C] == p4);
    BOOST_CHECK(simplex.m_a[idx_C] == a4);
    BOOST_CHECK(simplex.m_b[idx_C] == b4);

    BOOST_CHECK(!convex::is_point_in_simplex(not_in_simplex, simplex));
    BOOST_CHECK(convex::is_point_in_simplex(p1, simplex));
    BOOST_CHECK(!convex::is_point_in_simplex(p2, simplex));
    BOOST_CHECK(convex::is_point_in_simplex(p3, simplex));
    BOOST_CHECK(convex::is_point_in_simplex(p4, simplex));

    // Let us erase one more simplex

    simplex.m_bitmask = simplex.m_bitmask & ~bit_B;

    BOOST_CHECK(!convex::is_full_simplex(simplex));
    BOOST_CHECK(convex::dimension(simplex) == 2u);

    bit_A = 0xFFFF;
    idx_A = 0xFFFF;
    BOOST_CHECK_NO_THROW(convex::get_used_indices(simplex.m_bitmask, idx_A, bit_A));
    BOOST_CHECK(bit_A == 1);
    BOOST_CHECK(idx_A == 0);
    BOOST_CHECK(simplex.m_v[idx_A] == p1);
    BOOST_CHECK(simplex.m_a[idx_A] == a1);
    BOOST_CHECK(simplex.m_b[idx_A] == b1);

    bit_A = 0xFFFF;
    idx_A = 0xFFFF;
    bit_B = 0xFFFF;
    idx_B = 0xFFFF;
    BOOST_CHECK_NO_THROW(convex::get_used_indices(simplex.m_bitmask, idx_A, bit_A, idx_B, bit_B));

    BOOST_CHECK(bit_A == 1);
    BOOST_CHECK(idx_A == 0);
    BOOST_CHECK(simplex.m_v[idx_A] == p1);
    BOOST_CHECK(simplex.m_a[idx_A] == a1);
    BOOST_CHECK(simplex.m_b[idx_A] == b1);

    BOOST_CHECK(bit_B == 8);
    BOOST_CHECK(idx_B == 3);
    BOOST_CHECK(simplex.m_v[idx_B] == p4);
    BOOST_CHECK(simplex.m_a[idx_B] == a4);
    BOOST_CHECK(simplex.m_b[idx_B] == b4);

    bit_A = 0xFFFF;
    idx_A = 0xFFFF;
    bit_B = 0xFFFF;
    idx_B = 0xFFFF;
    bit_C = 0xFFFF;
    idx_C = 0xFFFF;
    BOOST_CHECK_THROW(convex::get_used_indices(simplex.m_bitmask, idx_A, bit_A, idx_B, bit_B, idx_C, bit_C),
                      std::logic_error);

    BOOST_CHECK(!convex::is_point_in_simplex(not_in_simplex, simplex));
    BOOST_CHECK(convex::is_point_in_simplex(p1, simplex));
    BOOST_CHECK(!convex::is_point_in_simplex(p2, simplex));
    BOOST_CHECK(!convex::is_point_in_simplex(p3, simplex));
    BOOST_CHECK(convex::is_point_in_simplex(p4, simplex));

    // Insert a new simplex vertex

    BOOST_CHECK_NO_THROW(convex::add_point_to_simplex(p5, a5, b5, simplex));

    BOOST_CHECK(!convex::is_full_simplex(simplex));
    BOOST_CHECK(convex::dimension(simplex) == 3u);

    bit_A = 0xFFFF;
    idx_A = 0xFFFF;
    BOOST_CHECK_NO_THROW(convex::get_used_indices(simplex.m_bitmask, idx_A, bit_A));
    BOOST_CHECK(bit_A == 1);
    BOOST_CHECK(idx_A == 0);
    BOOST_CHECK(simplex.m_v[idx_A] == p1);
    BOOST_CHECK(simplex.m_a[idx_A] == a1);
    BOOST_CHECK(simplex.m_b[idx_A] == b1);

    bit_A = 0xFFFF;
    idx_A = 0xFFFF;
    bit_B = 0xFFFF;
    idx_B = 0xFFFF;
    BOOST_CHECK_NO_THROW(convex::get_used_indices(simplex.m_bitmask, idx_A, bit_A, idx_B, bit_B));

    BOOST_CHECK(bit_A == 1);
    BOOST_CHECK(idx_A == 0);
    BOOST_CHECK(simplex.m_v[idx_A] == p1);
    BOOST_CHECK(simplex.m_a[idx_A] == a1);
    BOOST_CHECK(simplex.m_b[idx_A] == b1);

    BOOST_CHECK(bit_B == 2);
    BOOST_CHECK(idx_B == 1);
    BOOST_CHECK(simplex.m_v[idx_B] == p5);
    BOOST_CHECK(simplex.m_a[idx_B] == a5);
    BOOST_CHECK(simplex.m_b[idx_B] == b5);

    bit_A = 0xFFFF;
    idx_A = 0xFFFF;
    bit_B = 0xFFFF;
    idx_B = 0xFFFF;
    bit_C = 0xFFFF;
    idx_C = 0xFFFF;
    BOOST_CHECK_NO_THROW(convex::get_used_indices(simplex.m_bitmask, idx_A, bit_A, idx_B, bit_B, idx_C, bit_C));

    BOOST_CHECK(bit_A == 1);
    BOOST_CHECK(idx_A == 0);
    BOOST_CHECK(simplex.m_v[idx_A] == p1);
    BOOST_CHECK(simplex.m_a[idx_A] == a1);
    BOOST_CHECK(simplex.m_b[idx_A] == b1);

    BOOST_CHECK(bit_B == 2);
    BOOST_CHECK(idx_B == 1);
    BOOST_CHECK(simplex.m_v[idx_B] == p5);
    BOOST_CHECK(simplex.m_a[idx_B] == a5);
    BOOST_CHECK(simplex.m_b[idx_B] == b5);

    BOOST_CHECK(bit_C == 8);
    BOOST_CHECK(idx_C == 3);
    BOOST_CHECK(simplex.m_v[idx_C] == p4);
    BOOST_CHECK(simplex.m_a[idx_C] == a4);
    BOOST_CHECK(simplex.m_b[idx_C] == b4);

    BOOST_CHECK(!convex::is_point_in_simplex(not_in_simplex, simplex));
    BOOST_CHECK(convex::is_point_in_simplex(p1, simplex));
    BOOST_CHECK(!convex::is_point_in_simplex(p2, simplex));
    BOOST_CHECK(!convex::is_point_in_simplex(p3, simplex));
    BOOST_CHECK(convex::is_point_in_simplex(p4, simplex));
    BOOST_CHECK(convex::is_point_in_simplex(p5, simplex));
}

// 2011-11-12: Kenny: The degenerate tests always test if a tetrahedron will be degenerate, what about triangle and line cases?
// 2011-11-12: Kenny: It seems to me that you assume some order on the simplex vertices. In principle one could have bitmask = 11 = 13 = 14 = 14, corresponding to an used mask 0111, 1011, 1101, 1110. Does this order have any influence on your test? The same goes for triangle and line cases.
// 2011-11-12: Kenny: There is no tests that confirm non-degenerate point is correctly classified.

BOOST_AUTO_TEST_CASE(simplex_degenerate_point)
{
    using T = double;
    typedef convex::Simplex<T> S;

    S simplex;

    // Vertify that simplex have been initialized correctly as being empty and ``zeroed''
    BOOST_CHECK(simplex.m_bitmask == 0u);

    for (size_t i = 0u; i < 100u; ++i)
    {
        const EigenVector3<T> rand_in0 = randomEigen<T>();
        const EigenVector3<T> rand_in1 = randomEigen<T>();
        const EigenVector3<T> rand_in2 = randomEigen<T>();

        simplex.m_v[0] = rand_in0;
        simplex.m_v[1] = rand_in1;
        simplex.m_v[2] = rand_in2;
        simplex.m_bitmask = 7;

        bool const is_degenerate = convex::is_degenerate_point(rand_in1, simplex);
        BOOST_CHECK(is_degenerate == true);
    }

    for (size_t i = 0u; i < 100u; ++i)
    {
        const EigenVector3<T> rand_in0 = randomEigen<T>();
        const EigenVector3<T> rand_in1 = randomEigen<T>();
        const EigenVector3<T> rand_in2 = randomEigen<T>();
        simplex.m_v[0] = rand_in0;
        simplex.m_v[1] = rand_in1;
        simplex.m_v[2] = rand_in2;
        simplex.m_bitmask = 7;

        // 2011-11-12: Kenny Why not use tiny::Random value; w0 = value(); ????
        T w0 = std::rand() / RAND_MAX;
        T w1 = std::rand() / RAND_MAX;

        // 2011-11-12: Kenny Why use a while loop? You could just do tiny::Random value = tiny::Random( 0, 0.5) or someting like that? Or just divide by some positive number larger than or equal 2?
        // 2011-11-12: Kenny Why not use value_traits or hardwire 1 to the floating point precision used in the test?
        while (w0 + w1 >= 1)
        {
            w0 = std::rand() / RAND_MAX;
            w1 = std::rand() / RAND_MAX;
        }

        // 2011-11-12: Kenny Why not use value_traits or hardwire 1 to the floating point precision used in the test?
        // 2011-11-12: Kenny: If by design you initialize w0 and w1 to sum to one instead of using the while loop then all w's could be const declared.
        T w2 = 1 - w0 - w1;

        const EigenVector3<T> new_in_plane
            = w0 * rand_in0 + w1 * rand_in1 + w2 * rand_in2;

        bool const is_degenerate = convex::is_degenerate_point(new_in_plane, simplex);
        BOOST_CHECK(is_degenerate == true);
    }
}

BOOST_AUTO_TEST_SUITE_END();
