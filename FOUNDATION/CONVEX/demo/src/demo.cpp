#include <geometry.h>  // needed for geometry::Sphere
#include <convex.h>

#include <tiny.h>  // for math types

#include <cmath>  // for std::sqrt
#include <iostream>

void gjk_demo()
{
  std::cout << "GJK Demo" << std::endl;

  // Create some math types
  typedef tiny::MathTypes<double> M;
  
  typedef M::quaternion_type                      Q;
  typedef M::vector3_type                         V;
  typedef M::real_type                            T;
  typedef M::coordsys_type                        C;
  typedef M::value_traits                         VT;
  
  
  // Parameters controlling the behaviour of the GJK 
  size_t const max_iterations       = 100u;
  T      const absolute_tolerance   = VT::numeric_cast(10e-6);
  T      const relative_tolerance   = VT::numeric_cast(10e-6);
  T      const stagnation_tolerance = VT::numeric_cast(10e-15);
  
  // Define the geometries 
  geometry::Sphere<V> const A;
  geometry::Sphere<V> const B;

  // Define the coordinate transformations for the geometries
  C X_A;
  C X_B;

  X_A.T() = V::make( -2.1,12.0,1.0);
  X_A.Q() = Q::identity();
  X_B.T() = V::make( 4.1,11.3,1.0);
  X_B.Q() = Q::identity();

  // Create variables to hold results from the computation
  V p_A;
  V p_B;
  size_t iterations  = 0u;
  size_t status      = 0u;
  T      distance    = VT::infinity();
  
  // Ask GJK for closest points between the two geometries
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
  
  std::cout << "\tstatus = " << convex::get_status_message(status) << std::endl;
  std::cout << "\tdistance = " << distance << std::endl;
  std::cout << "\titerations = " << iterations << std::endl;
  
  std::cout << "Goodbye..." << std::endl;
}

void continuous_demo()
{
  using std::sqrt;
  
  std::cout << "CONTINUOUS Demo" << std::endl;
  
  // Create some math types
  typedef tiny::MathTypes<double> M;
  
  typedef M::quaternion_type                      Q;
  typedef M::vector3_type                         V;
  typedef M::real_type                            T;
  typedef M::coordsys_type                        X;
  typedef M::value_traits                         VT;
  
  // Create some geometries
  convex::Cylinder<M> A;
  geometry::Sphere<V> B;
  
  A.half_height() = 10.0;
  A.radius()      =  1.0;
  B.radius()      =  1.0;
  
  // Create from and to positions of the geometries
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
  
  // Set up parameters and create variables for results
  size_t const max_iterations = 100u;
  T      const epsilon        = 0.0001;
  T            toi            = 0.0;
  size_t iterations;
  V p_a;
  V p_b;
  
  // Call the motion interpolate function to find out if there is an impact and get the estimate of time-of-impact (toi)
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
  
  std::cout << "\timpact = " << impact << std::endl;
  std::cout << "\ttoi = " << toi << std::endl;
  std::cout << "\tdistance at toi = " << norm(p_a-p_b) << std::endl;
  
  std::cout << "Goodbye..." << std::endl;
}

int main()
{
  std::cout << "CONVEX Library Demo" << std::endl;
  
  gjk_demo();
  continuous_demo();
  
	return 0;
}
