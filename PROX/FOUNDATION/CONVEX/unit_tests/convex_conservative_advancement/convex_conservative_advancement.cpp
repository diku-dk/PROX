#include <tiny.h>
#include <geometry.h>  // needed for geometry::Sphere
#include <convex.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

#include <cmath>   // for std::fabs and std::sqrt

BOOST_AUTO_TEST_SUITE(conservative_advancement);

BOOST_AUTO_TEST_CASE(case_by_case_test)
{
    using std::fabs;
    
    typedef tiny::MathTypes<double>                 M;
    typedef M::value_traits                         VT;
    typedef M::vector3_type                         V;
    typedef M::real_type                            T;
    typedef M::quaternion_type                      Q;
    typedef M::coordsys_type                        C;
    
    // Setup a central impact, two sphere hitting in perfect symmetry, totally independent of their rotationnal motion!
    {
        geometry::Sphere<V> const A;
        geometry::Sphere<V> const B;
        
        C T_a;
        C T_b;
        
        T_a.T() = V::make(-2.0, 0.0, 0.0);
        T_a.Q() = Q::identity();
        V const v_a   = V::make( 2.0, 0.0, 0.0);    
        V const w_a   = V::make( 0.0, 5.0, 0.0);    
        T const r_max_a = A.radius();               
        
        T_b.T() = V::make( 2.0, 0.0, 0.0);
        T_b.Q() = Q::identity();
        V const v_b   = V::make(-2.0, 0.0, 0.0);    
        V const w_b   = V::make( 0.0, 0.0, 5.0);    
        T const r_max_b = B.radius();               
        
        size_t const max_iterations = 100u;

        T      const epsilon = 0.01;
        T            toi     = 0.0;
        size_t       iterations;
        V p_a;
        V p_b;
        
        bool const impact = convex::conservative_advancement<M>(
                                                                T_a
                                                                , v_a
                                                                , w_a
                                                                , &A
                                                                , r_max_a
                                                                , T_b
                                                                , v_b
                                                                , w_b
                                                                , &B
                                                                , r_max_b
                                                                , p_a
                                                                , p_b
                                                                , toi
                                                                , iterations
                                                                , epsilon
                                                                , VT::one()
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

        C T_a;
        C T_b;
        
        T_a.T()   = V::make(-2.0, 0.0, 0.0);
        T_a.Q()   = Q::identity();
        V const v_a     = V::make( -2.0, 0.0, 0.0);  
        V const w_a     = V::make( 0.0, 0.0, 0.0);   
        T const r_max_a = A.radius();                
        
        T_b.T()   = V::make( 2.0, 0.0, 0.0);
        T_b.Q()   = Q::identity();
        V const v_b     = V::make( 2.0, 0.0, 0.0);   
        V const w_b     = V::make( 0.0, 0.0, 0.0);   
        T const r_max_b = B.radius();                
        
        size_t const max_iterations = 100u;
        T      const epsilon        = 0.01;
        T            toi            = 0.0;
        size_t iterations;
        V p_a;
        V p_b;
        
        bool const  impact = convex::conservative_advancement<M>(
                                                                 T_a
                                                                 , v_a
                                                                 , w_a
                                                                 , &A
                                                                 , r_max_a
                                                                 , T_b
                                                                 , v_b
                                                                 , w_b
                                                                 , &B
                                                                 , r_max_b
                                                                 , p_a
                                                                 , p_b
                                                                 , toi
                                                                 , iterations
                                                                 , epsilon
                                                                 , VT::one()
                                                                 , max_iterations
                                                                 );
        
        BOOST_CHECK( !impact );
    }
    
    // Setup two sphere moving close by each other but never impacting
    {
      geometry::Sphere<V> const A;
      geometry::Sphere<V> const B;

        C T_a;
        C T_b;
        
        T_a.T()   = V::make(-2.0, 1.01, 0.0);
        T_a.Q()   = Q::identity();
        V const v_a     = V::make( -2.0, 0.0, 0.0);   
        V const w_a     = V::make( 0.0, 5.0, 0.0);    
        T const r_max_a = A.radius();                 
        
        T_b.T()   = V::make( 2.0, -1.01, 0.0);
        T_b.Q()   = Q::identity();
        V const v_b     = V::make( 2.0, 0.0, 0.0);    
        V const w_b     = V::make( 0.0, 0.0, 5.0);    
        T const r_max_b = B.radius();                 
        
        size_t const max_iterations = 100u;
        T      const epsilon        = 0.01;
        T            toi            = 0.0;
        size_t iterations;
        V p_a;
        V p_b;
        
        bool const impact = convex::conservative_advancement<M>(
                                                                T_a
                                                                , v_a
                                                                , w_a
                                                                , &A
                                                                , r_max_a
                                                                , T_b
                                                                , v_b
                                                                , w_b
                                                                , &B
                                                                , r_max_b
                                                                , p_a
                                                                , p_b
                                                                , toi
                                                                , iterations
                                                                , epsilon
                                                                , VT::one()
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
        
        C T_a;
        C T_b;
        
        T_a.T()   = V::make(0.0, 0.0, 2.0);
        T_a.Q()   = Q::identity();
        V const v_a     = V::make( 0.0, 0.0, 0.0);                                              
        V const w_a     = V::make( -VT::pi(), 0.0, 0.0 );                                       
        T const r_max_a = sqrt( A.half_height()*A.half_height() + A.radius()*A.radius() );      
        
        T_b.T()   = V::make( 0.0, 10.0, 0.0);
        T_b.Q()   = Q::identity();
        V const v_b     = V::make( 0.0, 0.0, 0.0);               
        V const w_b     = V::make( 0.0, 0.0, 0.0);               
        T const r_max_b = B.radius();                            
        
        size_t const max_iterations = 100u;
        T      const epsilon        = 0.01;
        T            toi            = 0.0;
        size_t iterations;
        V p_a;
        V p_b;
        
        bool const impact = convex::conservative_advancement<M>(
                                                                T_a
                                                                , v_a
                                                                , w_a
                                                                , &A
                                                                , r_max_a
                                                                , T_b
                                                                , v_b
                                                                , w_b
                                                                , &B
                                                                , r_max_b
                                                                , p_a
                                                                , p_b
                                                                , toi
                                                                , iterations
                                                                , epsilon
                                                                , VT::one()
                                                                , max_iterations
                                                                );
        
        BOOST_CHECK( impact );
        BOOST_CHECK_CLOSE( toi, 0.5, 1.0 );
        BOOST_CHECK( tiny::norm(p_a-p_b) < epsilon );
    }
}

BOOST_AUTO_TEST_SUITE_END();
