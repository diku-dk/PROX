#include <geometry.h>
#include <eigenhelperall.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(geometry);

BOOST_AUTO_TEST_CASE(inside_cylinder_test)
{
    using T = double;

    EigenVector3<T> const center = EigenVector3<T>(1.0, 2.0, 4.0);
    EigenVector3<T> const axis = EigenVector3<T>(0.0, 0.0, 1.0);
    T const radius = 2.0;
    T const height = 4.0;

    geometry::CylinderEigen<T> const cylinder
        = geometry::make_cylinder(radius, height, axis, center);

    BOOST_CHECK(geometry::is_valid(cylinder));

  // on bottom surface
    {
        EigenVector3<T> const p = EigenVector3<T>(1.0, 2.0, 2.0);
        bool const test = geometry::inside_cylinder(p, cylinder);
        BOOST_CHECK(test);
    }

  // on top surface
    {
        EigenVector3<T> const p = EigenVector3<T>(1.0, 2.0, 6.0);
        bool const test = geometry::inside_cylinder(p, cylinder);
        BOOST_CHECK(test);
    }

  // on cylinder surface
    {
        EigenVector3<T> const p = EigenVector3<T>(3.0, 2.0, 4.0);
        bool const test = geometry::inside_cylinder(p, cylinder);
        BOOST_CHECK(test);
    }

  // inside
    {
        EigenVector3<T> const p = EigenVector3<T>(1.0, 2.0, 2.1);
        bool const test = geometry::inside_cylinder(p, cylinder);
        BOOST_CHECK(test);
    }

  // inside
    {
        EigenVector3<T> const p = EigenVector3<T>(1.0, 2.0, 5.9);
        bool const test = geometry::inside_cylinder(p, cylinder);
        BOOST_CHECK(test);
    }

  // inside
    {
        EigenVector3<T> const p = EigenVector3<T>(2.9, 2.0, 4.0);
        bool const test = geometry::inside_cylinder(p, cylinder);
        BOOST_CHECK(test);
    }

  // outside
    {
        EigenVector3<T> const p = EigenVector3<T>(1.0, 2.0, 1.9);
        bool const test = geometry::inside_cylinder(p, cylinder);
        BOOST_CHECK(!test);
    }

  // outside
    {
        EigenVector3<T> const p = EigenVector3<T>(1.0, 2.0, 6.1);
        bool const test = geometry::inside_cylinder(p, cylinder);
        BOOST_CHECK(!test);
    }

  // outside
    {
        EigenVector3<T> const p = EigenVector3<T>(3.1, 2.0, 4.0);
        bool const test = geometry::inside_cylinder(p, cylinder);
        BOOST_CHECK(!test);
    }
}

BOOST_AUTO_TEST_SUITE_END();
