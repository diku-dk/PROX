#include <tiny.h>
#include <geometry.h>   // needed for geometry::Box, geometry::Sphere
#include <convex_shapes.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(convex_shapes);

BOOST_AUTO_TEST_CASE(case_by_case_testing)
{
    typedef tiny::MathTypes<double> M;
    typedef M::vector3_type V;

    geometry::Box<V> const box;
    convex::Capsule<M> const capsule;
    convex::Cone<M> const cone;
    convex::ConvexHull<M> hull;
    convex::Cylinder<M> const cylinder;
    convex::Ellipsoid<M> const ellipsoid;
    geometry::Sphere<V> const sphere;

    EigenVector3<double> zero{0.0, 0.0, 0.0};
    EigenVector3<double> up{0.0, 0.0, 1.0};
    EigenVector3<double> down{0.0, 0.0, -1.0};
    EigenVector3<double> right{1.0, 0.0, 0.0};
    EigenVector3<double> left{-1.0, 0.0, 0.0};
    EigenVector3<double> forward{0.0, 1.0, 0.0};
    EigenVector3<double> backward{0.0, -1.0, 0.0};

  // Test what happens if zero vector is given as search direction
    {
        auto p = box.get_support_point(zero);
        BOOST_CHECK((p(0) != 0.0) || (p(1) != 0) || (p(2) != 0));

        p = capsule.get_support_point(zero);
        BOOST_CHECK((p(0) != 0.0) || (p(1) != 0) || (p(2) != 0));

        p = cone.get_support_point(zero);
        BOOST_CHECK((p(0) != 0.0) || (p(1) != 0) || (p(2) != 0));

        p = cylinder.get_support_point(zero);
        BOOST_CHECK((p(0) != 0.0) || (p(1) != 0) || (p(2) != 0));

        p = ellipsoid.get_support_point(zero);
        BOOST_CHECK((p(0) != 0.0) || (p(1) != 0) || (p(2) != 0));

        p = sphere.get_support_point(zero);
        BOOST_CHECK((p(0) != 0.0) || (p(1) != 0) || (p(2) != 0));
    }

  // Test what happens if up vector is given
    {
        auto p = box.get_support_point(up);

        BOOST_CHECK(p(0) <= box.half_extent()(0));
        BOOST_CHECK(p(1) <= box.half_extent()(1));
        BOOST_CHECK(p(0) >= -box.half_extent()(0));
        BOOST_CHECK(p(1) >= -box.half_extent()(1));
        BOOST_CHECK(p(2) == box.half_extent()(2));
    }
    {
        auto p = capsule.get_support_point(up);
        BOOST_CHECK(p(0) == 0.0);
        BOOST_CHECK(p(1) == 0.0);
        BOOST_CHECK(p(2) == 2.0);
    }
    {
        auto p = cone.get_support_point(up);
        BOOST_CHECK(p(0) == 0.0);
        BOOST_CHECK(p(1) == 0.0);
        BOOST_CHECK(p(2) == 1.0);
    }
    {
        auto p = cylinder.get_support_point(up);
        BOOST_CHECK(p(0) == 0.0);
        BOOST_CHECK(p(1) == 0.0);
        BOOST_CHECK(p(2) == 1.0);
    }
    {
        auto p = ellipsoid.get_support_point(up);
        BOOST_CHECK(p(0) == 0.0);
        BOOST_CHECK(p(1) == 0.0);
        BOOST_CHECK(p(2) == 1.0);
    }
    {
        auto p = sphere.get_support_point(up);
        BOOST_CHECK(p(0) == 0.0);
        BOOST_CHECK(p(1) == 0.0);
        BOOST_CHECK(p(2) == 1.0);
    }

  // Test what happen if right vector is given
    {
        auto p = box.get_support_point(right);

        BOOST_CHECK(p(0) == box.half_extent()(0));
        BOOST_CHECK(p(1) <= box.half_extent()(1));
        BOOST_CHECK(p(2) <= box.half_extent()(2));
        BOOST_CHECK(p(1) >= -box.half_extent()(1));
        BOOST_CHECK(p(2) >= -box.half_extent()(2));
    }
    {
        auto p = capsule.get_support_point(right);
        BOOST_CHECK(p(0) == 1.0);
        BOOST_CHECK(p(1) == 0.0);
        BOOST_CHECK(p(2) == 0.0);
    }
    {
        auto p = cone.get_support_point(right);
        BOOST_CHECK(p(0) == 1.0);
        BOOST_CHECK(p(1) == 0.0);
        BOOST_CHECK(p(2) == -1.0);
    }
    {
        auto p = cylinder.get_support_point(right);
        BOOST_CHECK(p(0) == cylinder.radius());
        BOOST_CHECK(p(1) == 0.0);
        BOOST_CHECK(p(2) <= cylinder.half_height());
        BOOST_CHECK(p(2) >= -cylinder.half_height());
    }
    {
        auto p = ellipsoid.get_support_point(right);
        BOOST_CHECK(p(0) == 1.0);
        BOOST_CHECK(p(1) == 0.0);
        BOOST_CHECK(p(2) == 0.0);
    }
    {
        auto p = sphere.get_support_point(right);
        BOOST_CHECK(p(0) == 1.0);
        BOOST_CHECK(p(1) == 0.0);
        BOOST_CHECK(p(2) == 0.0);
    }

  // Test what happen if down vector is given
    {
        auto p = box.get_support_point(down);

        BOOST_CHECK(p(2) == -box.half_extent()(2));
        BOOST_CHECK(p(0) <= box.half_extent()(0));
        BOOST_CHECK(p(1) <= box.half_extent()(1));
        BOOST_CHECK(p(0) >= -box.half_extent()(0));
        BOOST_CHECK(p(1) >= -box.half_extent()(1));
    }
    {
        auto p = capsule.get_support_point(down);
        BOOST_CHECK(p(0) == 0.0);
        BOOST_CHECK(p(1) == 0.0);
        BOOST_CHECK(p(2) == -2.0);
    }
    {
        auto p = cone.get_support_point(down);
        BOOST_CHECK(p(0) == 0.0);
        BOOST_CHECK(p(1) == 0.0);
        BOOST_CHECK(p(2) == -1.0);
    }
    {
        auto p = cylinder.get_support_point(down);
        BOOST_CHECK(p(0) == 0.0);
        BOOST_CHECK(p(1) == 0.0);
        BOOST_CHECK(p(2) == -cylinder.half_height());
    }
    {
        auto p = ellipsoid.get_support_point(down);
        BOOST_CHECK(p(0) == 0.0);
        BOOST_CHECK(p(1) == 0.0);
        BOOST_CHECK(p(2) == -1.0);
    }
    {
        auto p = sphere.get_support_point(down);
        BOOST_CHECK(p(0) == 0.0);
        BOOST_CHECK(p(1) == 0.0);
        BOOST_CHECK(p(2) == -1.0);
    }
  // Test the convex hull
    {
        hull.add_point(left);
        hull.add_point(right);
        hull.add_point(backward);
        hull.add_point(forward);
        hull.add_point(down);
        hull.add_point(up);

        auto p = hull.get_support_point(right);
        BOOST_CHECK(p(0) == 1.0);
        BOOST_CHECK(p(1) == 0.0);
        BOOST_CHECK(p(2) == 0.0);

        p = hull.get_support_point(left);
        BOOST_CHECK(p(0) == -1.0);
        BOOST_CHECK(p(1) == 0.0);
        BOOST_CHECK(p(2) == 0.0);

        p = hull.get_support_point(forward);
        BOOST_CHECK(p(0) == 0.0);
        BOOST_CHECK(p(1) == 1.0);
        BOOST_CHECK(p(2) == 0.0);

        p = hull.get_support_point(backward);
        BOOST_CHECK(p(0) == 0.0);
        BOOST_CHECK(p(1) == -1.0);
        BOOST_CHECK(p(2) == 0.0);

        p = hull.get_support_point(up);
        BOOST_CHECK(p(0) == 0.0);
        BOOST_CHECK(p(1) == 0.0);
        BOOST_CHECK(p(2) == 1.0);

        p = hull.get_support_point(down);
        BOOST_CHECK(p(0) == 0.0);
        BOOST_CHECK(p(1) == 0.0);
        BOOST_CHECK(p(2) == -1.0);
    }
}

BOOST_AUTO_TEST_SUITE_END();
