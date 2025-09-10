#include <geometry.h>
#include <tiny.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

using T = float;

BOOST_AUTO_TEST_SUITE(geometry);

BOOST_AUTO_TEST_CASE(inside_obb_test)
{

    {
        EigenVector3<T> const center = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_ext = EigenVector3<T>(1.0, 1.0, 1.0);
        const EigenQuaternion<T> q = EigenQuaternion<T>::Identity();

        geometry::OBBEigen<T> A = geometry::make_obb<T>(center, q, half_ext);

        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenVector3<T> const p2 = EigenVector3<T>(2.0, 2.0, 2.0);

        BOOST_CHECK(geometry::inside_obb(p0, A) == true);
        BOOST_CHECK(geometry::inside_obb(p1, A) == true);
        BOOST_CHECK(geometry::inside_obb(p2, A) == false);
    }

    {
        EigenVector3<T> const center = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_ext = EigenVector3<T>(1.0, 1.0, 1.0);
        const EigenQuaternion<T> q = Rotatex<T>(std::numbers::pi_v<T> * 0.5);

        geometry::OBBEigen<T> A = geometry::make_obb<T>(center, q, half_ext);

        EigenVector3<T> const p0
            = geometry::transform_from_obb(EigenVector3<T>(0.0, 0.0, 0.0), A);
        EigenVector3<T> const p1
            = geometry::transform_from_obb(EigenVector3<T>(1.0, 1.0, 1.0), A);
        EigenVector3<T> const p2
            = geometry::transform_from_obb(EigenVector3<T>(2.0, 2.0, 2.0), A);

        BOOST_CHECK(geometry::inside_obb(p0, A) == true);
        BOOST_CHECK(geometry::inside_obb(p1, A) == true);
        BOOST_CHECK(geometry::inside_obb(p2, A) == false);
    }
}

BOOST_AUTO_TEST_SUITE_END();
