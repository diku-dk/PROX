#include <tiny.h>
#include <geometry.h>  // Needed for geometry::Sphere
#include <convex.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

#include <cmath>   // for std::sqrt

BOOST_AUTO_TEST_SUITE(conservative_advancement);

BOOST_AUTO_TEST_CASE(case_by_case_test)
{
  typedef tiny::MathTypes<double>                 M;
  typedef M::value_traits                         VT;
  typedef M::vector3_type                         V;
  typedef M::real_type                            T;
  typedef M::quaternion_type                      Q;
  typedef M::coordsys_type                        X;

  // Setup a central impact, two sphere hitting in perfect symmetry, totally independent of their rotationnal motion!
  {
    geometry::Sphere<V> const A;
    geometry::Sphere<V> const B;

    X A_from, A_to;
    X B_from, B_to;

    A_from.T() = V::make(-2.0, 0.0, 0.0);
    A_from.Q() = Q::identity();
    A_to.T()   = V::make(0.0, 0.0, 0.0);
    A_to.Q()   = Q::Ry( 5.0 );
    T r_max_a = A.radius();

    B_from.T() = V::make( 2.0, 0.0, 0.0);
    B_from.Q() = Q::identity();
    B_to.T()   = V::make( 0.0, 0.0, 0.0);
    B_to.Q()   = Q::Rz( 5.0 );
    T r_max_b = B.radius();

    size_t const max_iterations = 100u;
    T      const epsilon        = 0.01;
    T            toi            = 0.0;
    size_t iterations;
    V p_a;
    V p_b;

    bool impact = convex::motion_interpolation<M>(
                                                  A_from
                                                  , A_to
                                                  , &A
                                                  , r_max_a
                                                  , B_from
                                                  , B_to
                                                  , &B
                                                  , r_max_b
                                                  , p_a
                                                  , p_b
                                                  , toi
                                                  , iterations
                                                  , epsilon
                                                  , max_iterations
                                                  );

    BOOST_CHECK( impact );
    BOOST_CHECK_CLOSE( toi, 0.5, 1.0 );

    BOOST_CHECK( fabs( p_a(0) ) < epsilon );
    BOOST_CHECK( fabs( p_a(1) ) < epsilon );
    BOOST_CHECK( fabs( p_a(2) ) < epsilon );
    BOOST_CHECK( fabs( p_b(0) ) < epsilon );
    BOOST_CHECK( fabs( p_b(1) ) < epsilon );
    BOOST_CHECK( fabs( p_b(2) ) < epsilon );
  }

  // Setup separating motion, two spheres moving away from each other
  {
    geometry::Sphere<V> const A;
    geometry::Sphere<V> const B;

    X A_from, A_to;
    X B_from, B_to;

    A_from.T() = V::make(-2.0, 0.0, 0.0);
    A_from.Q() = Q::identity();
    A_to.T()   = V::make(-4.0, 0.0, 0.0);
    A_to.Q()   = Q::Ry( 5.0 );
    T r_max_a = A.radius();

    B_from.T() = V::make( 2.0, 0.0, 0.0);
    B_from.Q() = Q::identity();
    B_to.T()   = V::make( 4.0, 0.0, 0.0);
    B_to.Q()   = Q::Rz( 5.0 );
    T r_max_b = B.radius();

    size_t const max_iterations = 100u;
    T      const epsilon        = 0.01;
    T            toi            = 0.0;
    size_t iterations;
    V p_a;
    V p_b;

    bool impact = convex::motion_interpolation<M>(
                                                  A_from
                                                  , A_to
                                                  , &A
                                                  , r_max_a
                                                  , B_from
                                                  , B_to
                                                  , &B
                                                  , r_max_b
                                                  , p_a
                                                  , p_b
                                                  , toi
                                                  , iterations
                                                  , epsilon
                                                  , max_iterations
                                                  );

    BOOST_CHECK( !impact );
  }

  // Setup two sphere moving close by each other but never impacting
  {
    geometry::Sphere<V> const A;
    geometry::Sphere<V> const B;

    X A_from, A_to;
    X B_from, B_to;

    A_from.T() = V::make(-2.0, 1.01, 0.0);
    A_from.Q() = Q::identity();
    A_to.T()   = V::make(0.0, 1.01, 0.0);
    A_to.Q()   = Q::Ry( 5.0 );
    T r_max_a  = A.radius();

    B_from.T() = V::make( 2.0, -1.01, 0.0);
    B_from.Q() = Q::identity();
    B_to.T()   = V::make( 0.0, -1.01, 0.0);
    B_to.Q()   = Q::Rz( 5.0 );
    T r_max_b  = B.radius();

    size_t const max_iterations = 100u;
    T      const epsilon        = 0.01;
    T            toi            = 0.0;
    size_t iterations;
    V p_a;
    V p_b;

    bool impact = convex::motion_interpolation<M>(
                                                  A_from
                                                  , A_to
                                                  , &A
                                                  , r_max_a
                                                  , B_from
                                                  , B_to
                                                  , &B
                                                  , r_max_b
                                                  , p_a
                                                  , p_b
                                                  , toi
                                                  , iterations
                                                  , epsilon
                                                  , max_iterations
                                                  );

    BOOST_CHECK( !impact );
  }

  // Setup test case where only rotational motion is of importance!
  {
    using std::sqrt;

    convex::Cylinder<M> A;
    geometry::Sphere<V> B;

    A.half_height() = 10.0;
    A.radius()      =  1.0;
    B.radius()      =  1.0;

    X A_from, A_to;
    X B_from, B_to;

    A_from.T() = V::make(0.0, 0.0, 2.0);
    A_from.Q() = Q::identity();
    A_to.T()   = V::make(0.0, 0.0, 2.0);
    A_to.Q()   = Q::Rx( -VT::pi() );
    T r_max_a = sqrt( A.half_height()*A.half_height() + A.radius()*A.radius() );

    B_from.T() = V::make( 0.0, 10.0, 0.0);
    B_from.Q() = Q::identity();
    B_to.T()   = V::make( 0.0, 10.0, 0.0);
    B_to.Q()   = Q::identity();
    T r_max_b  = B.radius();

    size_t const max_iterations = 100u;
    T      const epsilon        = 0.01;
    T            toi            = 0.0;
    size_t iterations;
    V p_a;
    V p_b;

    bool impact = convex::motion_interpolation<M>(
                                                  A_from
                                                  , A_to
                                                  , &A
                                                  , r_max_a
                                                  , B_from
                                                  , B_to
                                                  , &B
                                                  , r_max_b
                                                  , p_a
                                                  , p_b
                                                  , toi
                                                  , iterations
                                                  , epsilon
                                                  , max_iterations
                                                  );

    BOOST_CHECK( impact );
    BOOST_CHECK_CLOSE( toi, 0.5, 1.0 );
    BOOST_CHECK( norm(p_a-p_b) < epsilon );
  }

}

BOOST_AUTO_TEST_SUITE_END();
