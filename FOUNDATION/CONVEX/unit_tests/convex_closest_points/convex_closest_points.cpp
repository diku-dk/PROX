#include <tiny.h>
#include <geometry.h>  // needed for geometry::Box and geometry::Sphere
#include <convex.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

#include <cmath>

typedef tiny::MathTypes<double>  M;
typedef M::quaternion_type       Q;
typedef M::vector3_type          V;
typedef M::real_type             T;
typedef M::coordsys_type         C;
typedef M::value_traits          VT;


BOOST_AUTO_TEST_SUITE(convex_compute_closest_points);

BOOST_AUTO_TEST_CASE(sphere_box_compile_test)
{
  {
    geometry::Sphere<V>  B;
    geometry::Box<V>     A;

    B.radius()      = 0.5;
    A.half_extent() = V::make(50, .5, 50);

    C X_A;
    C X_B;

    V p_A = V::zero();
    V p_B = V::zero();

    X_A = C::identity();

    X_B.T() = V::make(0, 1.2, 0);
    X_B.Q() = Q::identity();


    convex::compute_closest_points<M>(
                                      X_A
                                      , &A
                                      , X_B
                                      , &B
                                      , p_A
                                      , p_B
                                      );
  }
  {

    geometry::Box<V>       A;
    convex::Capsule<M>     B;

    A.half_extent() = V::make(50, 1, 50);
    B.radius()      = 0.5;

    C X_A;
    C X_B;

    V p_A = V::zero();
    V p_B = V::zero();

    X_A = C::identity();

    X_B.T() = V::make(0, 1.5, 0);
    X_B.Q() = Q::identity();


    convex::compute_closest_points<M>(
                                      X_A
                                      , &A
                                      , X_B
                                      , &B
                                      , p_A
                                      , p_B
                                      );

  }
}

BOOST_AUTO_TEST_CASE(sphere_sphere_testing)
{
  geometry::Sphere<V> const A;
  geometry::Sphere<V> const B;

  size_t const max_iterations       = 100u;
  T      const absolute_tolerance   = VT::numeric_cast(10e-6);
  T      const relative_tolerance   = VT::numeric_cast(10e-6);
  T      const stagnation_tolerance = VT::numeric_cast(10e-15);

  // Two unit-spheres placed ontop of each other
  {
    C X_A;
    C X_B;

    V p_A;
    V p_B;
    size_t iterations  = 0u;
    size_t status      = 0u;
    T      distance    = VT::infinity();


    X_A = C::identity();
    X_B = C::identity();
    X_B.T() = V::make(VT::zero(), VT::two(), VT::zero());

    convex::compute_closest_points<M>(
                                      X_A
                                      , &A
                                      , X_B
                                      , &B
                                      , p_A
                                      , p_B
                                      , distance
                                      , iterations
                                      , status
                                      , absolute_tolerance
                                      , relative_tolerance
                                      , stagnation_tolerance
                                      , max_iterations
                                      );

    T true_distance = tiny::norm( X_A.T() - X_B.T() ) - VT::two();

    BOOST_CHECK( fabs(distance - true_distance) < 10e-6);

    BOOST_CHECK( p_A(0) == VT::zero() );
    BOOST_CHECK( p_A(1) == VT::one() );
    BOOST_CHECK( p_A(2) == VT::zero() );

    BOOST_CHECK( p_B(0) == VT::zero() );
    BOOST_CHECK( p_B(1) == VT::one()  );
    BOOST_CHECK( p_B(2) == VT::zero() );
  }
  // Two unit-spheres overlapping but both placed on the x-axis
  {
    C X_A;
    C X_B;

    V p_A;
    V p_B;
    size_t iterations  = 0u;
    size_t status      = 0u;
    T      distance    = VT::infinity();

    X_A = C::identity();
    X_B = C::identity();
    X_B.T()(0) = VT::numeric_cast(2.5);

    convex::compute_closest_points<M>(
                                      X_A
                                      , &A
                                      , X_B
                                      , &B
                                      , p_A
                                      , p_B
                                      , distance
                                      , iterations
                                      , status
                                      , absolute_tolerance
                                      , relative_tolerance
                                      , stagnation_tolerance
                                      , max_iterations
                                      );

    T true_distance = tiny::norm( X_A.T() - X_B.T() ) - VT::two();

    BOOST_CHECK( fabs(distance - true_distance) < 10e-6);

    BOOST_CHECK( p_A(0) == VT::one()  );
    BOOST_CHECK( p_A(1) == VT::zero() );
    BOOST_CHECK( p_A(2) == VT::zero() );

    BOOST_CHECK( p_B(0) == VT::numeric_cast(1.5)  );
    BOOST_CHECK( p_B(1) == VT::zero() );
    BOOST_CHECK( p_B(2) == VT::zero() );
  }
  // Two unit-spheres exactly touching in one point (= one intersection point) and but both placed on the x-axis
  {
    C X_A;
    C X_B;

    V p_A;
    V p_B;
    size_t iterations  = 0u;
    size_t status      = 0u;
    T      distance    = VT::infinity();

    X_A.T().clear();
    X_A.T()(0) = -2.0;
    X_A.Q() = Q::identity();
    X_B.T().clear();
    X_B.Q() = Q::identity();

    convex::compute_closest_points<M>(
                                      X_A
                                      , &A
                                      , X_B
                                      , &B
                                      , p_A
                                      , p_B
                                      , distance
                                      , iterations
                                      , status
                                      , absolute_tolerance
                                      , relative_tolerance
                                      , stagnation_tolerance
                                      , max_iterations
                                      );

  }
  // Two unit-spheres non-overlapping but both placed on the x-axis
  {
    C X_A;
    C X_B;

    V p_A;
    V p_B;
    size_t iterations  = 0u;
    size_t status      = 0u;
    T      distance    = VT::infinity();

    X_A.T().clear();
    X_A.T()(0) = -2.5;
    X_A.Q() = Q::identity();
    X_B.T().clear();
    X_B.Q() = Q::identity();

    convex::compute_closest_points<M>(
                                      X_A
                                      , &A
                                      , X_B
                                      , &B
                                      , p_A
                                      , p_B
                                      , distance
                                      , iterations
                                      , status
                                      , absolute_tolerance
                                      , relative_tolerance
                                      , stagnation_tolerance
                                      , max_iterations
                                      );

  }
}

BOOST_AUTO_TEST_CASE(random_test)
{
  geometry::Sphere<V> const A;
  geometry::Sphere<V> const B;

  size_t const max_iterations       = 100u;
  T      const absolute_tolerance   = VT::numeric_cast(10e-6);
  T      const relative_tolerance   = VT::numeric_cast(10e-10);
  T      const stagnation_tolerance = VT::numeric_cast(0.0);

  for(size_t i=0;i<100u;++i)
  {
    C X_A;
    C X_B;

    V p_A;
    V p_B;
    size_t iterations  = 0u;
    size_t status      = 0u;
    T      distance    = VT::infinity();

    X_A.T() = V::random( -2.0, 2.0 );
    X_A.Q() = tiny::unit( Q::random() );
    X_B.T() = V::random( -2.0, 2.0 );
    X_B.Q() = tiny::unit( Q::random() );

    convex::compute_closest_points<M>(
                                      X_A
                                      , &A
                                      , X_B
                                      , &B
                                      , p_A
                                      , p_B
                                      , distance
                                      , iterations
                                      , status
                                      , absolute_tolerance
                                      , relative_tolerance
                                      , stagnation_tolerance
                                      , max_iterations
                                      );

    T const true_distance = tiny::norm( X_A.T() - X_B.T() ) - VT::two();

    if(  true_distance > absolute_tolerance )
    {
      BOOST_CHECK_CLOSE( true_distance, distance, 0.05 );

      BOOST_CHECK( status != convex::ABSOLUTE_CONVERGENCE );            // Would indicate penetration
      BOOST_CHECK( status != convex::INTERSECTION );                    // Would indicate penetration
      BOOST_CHECK( status != convex::ITERATING );                       // Would indicate internal error in CONVEX
      BOOST_CHECK( status != convex::EXCEEDED_MAX_ITERATIONS_LIMIT );   // Would indicate internal error in CONVEX
      BOOST_CHECK( status != convex::NON_DESCEND_DIRECTION );           // Would indicate internal error in CONVEX

      //BOOST_CHECK( status != convex::RELATIVE_CONVERGENCE );          // Would indicate convergence to positive distance
      //BOOST_CHECK( status != convex::SIMPLEX_EXPANSION_FAILED );      // Would indicate convergence to positive distance
      //BOOST_CHECK( status != convex::STAGNATION );                    // Would indicate convergence to positive distance
      //BOOST_CHECK( status != convex::LOWER_ERROR_BOUND_CONVERGENCE ); // Would indicate convergence to positive distance
    }
    else
    {
      BOOST_CHECK( 0.0 <= distance );
      BOOST_CHECK( distance <= absolute_tolerance );
      
      BOOST_CHECK( status != convex::STAGNATION );                     // Can only occur in case of positive distance
      BOOST_CHECK( status != convex::LOWER_ERROR_BOUND_CONVERGENCE );  // Can only occur in case of positive distance
      BOOST_CHECK( status != convex::RELATIVE_CONVERGENCE );           // Can only occur in case of positive distance
      BOOST_CHECK( status != convex::SIMPLEX_EXPANSION_FAILED );       // Can only occur in case of positive distance
      BOOST_CHECK( status != convex::ITERATING );                      // Would indicate internal error in CONVEX
      BOOST_CHECK( status != convex::EXCEEDED_MAX_ITERATIONS_LIMIT );  // Would indicate internal error in CONVEX
      BOOST_CHECK( status != convex::NON_DESCEND_DIRECTION );          // Would indicate internal error in CONVEX

      //BOOST_CHECK( status != convex::ABSOLUTE_CONVERGENCE );         // Indicates penetration
      //BOOST_CHECK( status != convex::INTERSECTION );                 // Indicates penetration
    }
    
  }
}

BOOST_AUTO_TEST_SUITE_END();
