#include <geometry.h>
#include <eigenhelperall.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

template <typename T, size_t N>
void test_direction_table(geometry::DirectionTable<T, N> const& DT)
{
    BOOST_CHECK_EQUAL(DT.size(), N);

    for (size_t i = 0u; i < N; ++i)
    {
        EigenVector3<T> const di = DT(i);

        BOOST_CHECK_CLOSE(dot(di, di), 1.0f, 0.01f);

        for (size_t j = i + 1u; j < N; ++j)
        {
            EigenVector3<T> const dj = DT(j);

            EigenVector3<T> diff = di - dj;

            BOOST_CHECK_GT(dot(diff, diff), 0.01f);
        }
    }
}

BOOST_AUTO_TEST_SUITE(geometry);

BOOST_AUTO_TEST_CASE(make_direction_tables)
{
    using T = float;

    {
        geometry::DirectionTable<T, 3> A = geometry::make3<T>();
        test_direction_table(A);
    }
    {
        geometry::DirectionTable<T, 4> A = geometry::make4<T>();
        test_direction_table(A);
    }
    {
        geometry::DirectionTable<T, 6> A = geometry::make6<T>();
        test_direction_table(A);
    }
    {
        geometry::DirectionTable<T, 7> A = geometry::make7<T>();
        test_direction_table(A);
    }
    {
        geometry::DirectionTable<T, 9> A = geometry::make9<T>();
        test_direction_table(A);
    }
    {
        geometry::DirectionTable<T, 10> A = geometry::make10<T>();
        test_direction_table(A);
    }
    {
        geometry::DirectionTable<T, 13> A = geometry::make13<T>();
        test_direction_table(A);
    }
}

BOOST_AUTO_TEST_SUITE_END();
