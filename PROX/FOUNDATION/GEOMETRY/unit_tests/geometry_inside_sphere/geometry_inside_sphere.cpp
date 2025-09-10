#include <geometry.h>
#include <eigenhelperall.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(geometry);

BOOST_AUTO_TEST_CASE(inside_sphere_test)
{
    using T = double;

    EigenVector3<T> const center = EigenVector3<T>(0, 0, 0);
    T const radius = 1;

    geometry::Sphere<T> const sphere = geometry::make_sphere(center, radius);

    BOOST_CHECK(geometry::is_valid(sphere));

    {
        EigenVector3<T> const p = EigenVector3<T>(0.0, 0.0, 0.0);
        bool const test = geometry::inside_sphere(p, sphere);
        BOOST_CHECK(test);
    }
    {
        EigenVector3<T> const p = EigenVector3<T>(1.0, 0.0, 0.0);
        bool const test = geometry::inside_sphere(p, sphere);
        BOOST_CHECK(test);
    }
    {
        EigenVector3<T> const p = EigenVector3<T>(2.0, 0.0, 0.0);
        bool const test = geometry::inside_sphere(p, sphere);
        BOOST_CHECK(!test);
    }
}

BOOST_AUTO_TEST_SUITE_END();
