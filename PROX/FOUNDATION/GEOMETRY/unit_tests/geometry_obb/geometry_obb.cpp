#include <geometry.h>
#include <tiny.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

using T = float;

BOOST_AUTO_TEST_SUITE(geometry);

BOOST_AUTO_TEST_CASE(obb_test)
{

    {
        geometry::OBBEigen<T> A;
        BOOST_CHECK(geometry::is_valid(A) == true);

        BOOST_CHECK_CLOSE(A.center()(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(A.center()(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(A.center()(2), 0.0, 0.01);

        BOOST_CHECK_CLOSE(A.half_extent()(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(A.half_extent()(1), 1.0, 0.01);
        BOOST_CHECK_CLOSE(A.half_extent()(2), 1.0, 0.01);

        BOOST_CHECK_CLOSE(A.orientation().w(), 1.0, 0.01);
        BOOST_CHECK_CLOSE(A.orientation().x(), 0.0, 0.01);
        BOOST_CHECK_CLOSE(A.orientation().y(), 0.0, 0.01);
        BOOST_CHECK_CLOSE(A.orientation().z(), 0.0, 0.01);

        EigenVector3<T> const p000 = geometry::get_local_corner(0, A);
        EigenVector3<T> const p001 = geometry::get_local_corner(1, A);
        EigenVector3<T> const p010 = geometry::get_local_corner(2, A);
        EigenVector3<T> const p011 = geometry::get_local_corner(3, A);
        EigenVector3<T> const p100 = geometry::get_local_corner(4, A);
        EigenVector3<T> const p101 = geometry::get_local_corner(5, A);
        EigenVector3<T> const p110 = geometry::get_local_corner(6, A);
        EigenVector3<T> const p111 = geometry::get_local_corner(7, A);

        BOOST_CHECK_CLOSE(p000(0), -1.0, 0.01);
        BOOST_CHECK_CLOSE(p000(1), -1.0, 0.01);
        BOOST_CHECK_CLOSE(p000(2), -1.0, 0.01);

        BOOST_CHECK_CLOSE(p001(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(p001(1), -1.0, 0.01);
        BOOST_CHECK_CLOSE(p001(2), -1.0, 0.01);

        BOOST_CHECK_CLOSE(p010(0), -1.0, 0.01);
        BOOST_CHECK_CLOSE(p010(1), 1.0, 0.01);
        BOOST_CHECK_CLOSE(p010(2), -1.0, 0.01);

        BOOST_CHECK_CLOSE(p011(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(p011(1), 1.0, 0.01);
        BOOST_CHECK_CLOSE(p011(2), -1.0, 0.01);

        BOOST_CHECK_CLOSE(p100(0), -1.0, 0.01);
        BOOST_CHECK_CLOSE(p100(1), -1.0, 0.01);
        BOOST_CHECK_CLOSE(p100(2), 1.0, 0.01);

        BOOST_CHECK_CLOSE(p101(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(p101(1), -1.0, 0.01);
        BOOST_CHECK_CLOSE(p101(2), 1.0, 0.01);

        BOOST_CHECK_CLOSE(p110(0), -1.0, 0.01);
        BOOST_CHECK_CLOSE(p110(1), 1.0, 0.01);
        BOOST_CHECK_CLOSE(p110(2), 1.0, 0.01);

        BOOST_CHECK_CLOSE(p111(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(p111(1), 1.0, 0.01);
        BOOST_CHECK_CLOSE(p111(2), 1.0, 0.01);
    }

    {
        EigenVector3<T> const center = EigenVector3<T>(1.0, 2.0, 3.0);
        EigenVector3<T> const half_ext = EigenVector3<T>(4.0, 5.0, 6.0);
        const EigenQuaternion<T> q = Rotatex(std::numbers::pi_v<T> * 0.5f);

        geometry::OBBEigen<T> A = geometry::make_obb<T>(center, q, half_ext);

        BOOST_CHECK(geometry::is_valid(A) == true);

        BOOST_CHECK_CLOSE(A.center()(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(A.center()(1), 2.0, 0.01);
        BOOST_CHECK_CLOSE(A.center()(2), 3.0, 0.01);

        BOOST_CHECK_CLOSE(A.half_extent()(0), 4.0, 0.01);
        BOOST_CHECK_CLOSE(A.half_extent()(1), 5.0, 0.01);
        BOOST_CHECK_CLOSE(A.half_extent()(2), 6.0, 0.01);

        BOOST_CHECK_CLOSE(A.orientation().w(), 0.707106769, 0.01);
        BOOST_CHECK_CLOSE(A.orientation().x(), 0.707106769, 0.01);
        BOOST_CHECK_CLOSE(A.orientation().y(), 0.0, 0.01);
        BOOST_CHECK_CLOSE(A.orientation().z(), 0.0, 0.01);
    }

    {
        EigenVector3<T> const center = EigenVector3<T>(1.0, 2.0, 3.0);
        EigenVector3<T> const half_ext = EigenVector3<T>(4.0, 5.0, 6.0);
        const EigenQuaternion<T> q = Rotatex(std::numbers::pi_v<T> * 0.5f);

        geometry::OBBEigen<T> const A
            = geometry::make_obb<T>(center, q, half_ext);
        geometry::OBBEigen<T> const B(A);
        geometry::OBBEigen<T> const C = B;

        BOOST_CHECK(geometry::is_valid(B) == true);
        BOOST_CHECK(geometry::is_valid(C) == true);

        BOOST_CHECK_CLOSE(B.center()(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(B.center()(1), 2.0, 0.01);
        BOOST_CHECK_CLOSE(B.center()(2), 3.0, 0.01);
        BOOST_CHECK_CLOSE(B.half_extent()(0), 4.0, 0.01);
        BOOST_CHECK_CLOSE(B.half_extent()(1), 5.0, 0.01);
        BOOST_CHECK_CLOSE(B.half_extent()(2), 6.0, 0.01);
        BOOST_CHECK_CLOSE(B.orientation().w(), 0.707106769, 0.01);
        BOOST_CHECK_CLOSE(B.orientation().x(), 0.707106769, 0.01);
        BOOST_CHECK_CLOSE(B.orientation().y(), 0.0, 0.01);
        BOOST_CHECK_CLOSE(B.orientation().z(), 0.0, 0.01);

        BOOST_CHECK_CLOSE(C.center()(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(C.center()(1), 2.0, 0.01);
        BOOST_CHECK_CLOSE(C.center()(2), 3.0, 0.01);
        BOOST_CHECK_CLOSE(C.half_extent()(0), 4.0, 0.01);
        BOOST_CHECK_CLOSE(C.half_extent()(1), 5.0, 0.01);
        BOOST_CHECK_CLOSE(C.half_extent()(2), 6.0, 0.01);
        BOOST_CHECK_CLOSE(C.orientation().w(), 0.707106769, 0.01);
        BOOST_CHECK_CLOSE(C.orientation().x(), 0.707106769, 0.01);
        BOOST_CHECK_CLOSE(C.orientation().y(), 0.0, 0.01);
        BOOST_CHECK_CLOSE(C.orientation().z(), 0.0, 0.01);
    }
}

BOOST_AUTO_TEST_SUITE_END();
