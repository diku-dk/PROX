#include <geometry.h>
#include <eigenhelperall.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

using T = float;

BOOST_AUTO_TEST_SUITE(geometry);

BOOST_AUTO_TEST_CASE(tetrahedron_test)
{

    {
        geometry::TetrahedronEigen<T> A;
        BOOST_CHECK(geometry::is_valid(A) == true);

        BOOST_CHECK_CLOSE(A.point(0)(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(A.point(0)(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(A.point(0)(2), 0.0, 0.01);

        BOOST_CHECK_CLOSE(A.point(1)(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(A.point(1)(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(A.point(1)(2), 0.0, 0.01);

        BOOST_CHECK_CLOSE(A.point(2)(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(A.point(2)(1), 1.0, 0.01);
        BOOST_CHECK_CLOSE(A.point(2)(2), 0.0, 0.01);

        BOOST_CHECK_CLOSE(A.point(3)(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(A.point(3)(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(A.point(3)(2), 1.0, 0.01);

        T const vol = geometry::get_signed_volume(A);

        BOOST_CHECK_CLOSE(vol, 1.0 / 6.0, 0.01);

        geometry::Triangle<T> const face_i = geometry::get_opposite_face(0, A);
        geometry::Triangle<T> const face_j = geometry::get_opposite_face(1, A);
        geometry::Triangle<T> const face_k = geometry::get_opposite_face(2, A);
        geometry::Triangle<T> const face_m = geometry::get_opposite_face(3, A);

    //  face opposite vertex 0:  123

        BOOST_CHECK_CLOSE(face_i.point(0)(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(face_i.point(0)(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(face_i.point(0)(2), 0.0, 0.01);

        BOOST_CHECK_CLOSE(face_i.point(1)(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(face_i.point(1)(1), 1.0, 0.01);
        BOOST_CHECK_CLOSE(face_i.point(1)(2), 0.0, 0.01);

        BOOST_CHECK_CLOSE(face_i.point(2)(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(face_i.point(2)(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(face_i.point(2)(2), 1.0, 0.01);

    //  face opposite vertex 1:  032
        BOOST_CHECK_CLOSE(face_j.point(0)(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(face_j.point(0)(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(face_j.point(0)(2), 0.0, 0.01);

        BOOST_CHECK_CLOSE(face_j.point(1)(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(face_j.point(1)(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(face_j.point(1)(2), 1.0, 0.01);

        BOOST_CHECK_CLOSE(face_j.point(2)(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(face_j.point(2)(1), 1.0, 0.01);
        BOOST_CHECK_CLOSE(face_j.point(2)(2), 0.0, 0.01);

    //  face opposite vertex 2:  013
        BOOST_CHECK_CLOSE(face_k.point(0)(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(face_k.point(0)(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(face_k.point(0)(2), 0.0, 0.01);

        BOOST_CHECK_CLOSE(face_k.point(1)(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(face_k.point(1)(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(face_k.point(1)(2), 0.0, 0.01);

        BOOST_CHECK_CLOSE(face_k.point(2)(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(face_k.point(2)(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(face_k.point(2)(2), 1.0, 0.01);

    //  face opposite vertex 3:  021
        BOOST_CHECK_CLOSE(face_m.point(0)(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(face_m.point(0)(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(face_m.point(0)(2), 0.0, 0.01);

        BOOST_CHECK_CLOSE(face_m.point(1)(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(face_m.point(1)(1), 1.0, 0.01);
        BOOST_CHECK_CLOSE(face_m.point(1)(2), 0.0, 0.01);

        BOOST_CHECK_CLOSE(face_m.point(2)(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(face_m.point(2)(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(face_m.point(2)(2), 0.0, 0.01);
    }

    {
        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::TetrahedronEigen<T> A
            = geometry::make_tetrahedron<T>(p0, p1, p2, p3);
        BOOST_CHECK(geometry::is_valid(A) == true);

        BOOST_CHECK_CLOSE(A.point(0)(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(A.point(0)(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(A.point(0)(2), 0.0, 0.01);

        BOOST_CHECK_CLOSE(A.point(1)(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(A.point(1)(1), 1.0, 0.01);
        BOOST_CHECK_CLOSE(A.point(1)(2), 0.0, 0.01);

        BOOST_CHECK_CLOSE(A.point(2)(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(A.point(2)(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(A.point(2)(2), 0.0, 0.01);

        BOOST_CHECK_CLOSE(A.point(3)(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(A.point(3)(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(A.point(3)(2), 1.0, 0.01);

        T const vol = geometry::get_signed_volume(A);

        BOOST_CHECK_CLOSE(vol, -1.0 / 6.0, 0.01);
    }

    {
        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::TetrahedronEigen<T> A
            = geometry::make_tetrahedron<T>(p0, p1, p2, p3);

        EigenVector3<T> const v = EigenVector3<T>(-1.0, -1.0, -1.0);
        EigenVector3<T> const p = A.get_support_point(v);

        BOOST_CHECK_CLOSE(p(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(p(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(p(2), 0.0, 0.01);
    }
    {
        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::TetrahedronEigen<T> A
            = geometry::make_tetrahedron(p0, p1, p2, p3);

        EigenVector3<T> const v = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p = A.get_support_point(v);

        BOOST_CHECK_CLOSE(p(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(p(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(p(2), 0.0, 0.01);
    }
    {
        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::TetrahedronEigen<T> A
            = geometry::make_tetrahedron(p0, p1, p2, p3);

        EigenVector3<T> const v = EigenVector3<T>(0.0, 0.0, 1.0);
        EigenVector3<T> const p = A.get_support_point(v);

        BOOST_CHECK_CLOSE(p(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(p(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(p(2), 1.0, 0.01);
    }
}

BOOST_AUTO_TEST_SUITE_END();
