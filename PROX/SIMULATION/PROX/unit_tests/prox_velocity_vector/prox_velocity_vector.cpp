#include <narrow.h>

#include <prox_rigid_body.h>
#include <prox_get_velocity_vector.h>
#include <prox_set_velocity_vector.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(velocity_vector);

BOOST_AUTO_TEST_CASE(get_velocity_vector_test)
{

    using T = float;
  // Initialize geometry
    narrow::System<T> narrow;

    size_t gid = narrow.create_geometry();

    narrow::Geometry<T>::sphere_type sphere;
    sphere.radius() = 1.0f;

    narrow::Geometry<T>& moving_geometry = narrow.get_geometry(gid);
    moving_geometry.add_shape(sphere);

  // Initialize objects
    std::vector<prox::RigidBody<T>> bodies;
    bodies.resize(4u);

    bodies[0].set_geometry_idx(gid);
    bodies[0].set_velocity(EigenVector3<T>(5.0, 7.0, 0.0));

    bodies[1].set_geometry_idx(gid);
    bodies[1].set_velocity(EigenVector3<T>(0.0, 1.0, 0.0));

    bodies[2].set_geometry_idx(gid);
    bodies[2].set_velocity(EigenVector3<T>(0.0, 3.0, 0.0));

    bodies[3].set_geometry_idx(gid);
    bodies[3].set_velocity(EigenVector3<T>(2.0, 5.0, 0.0));

    Eigen::VectorX<T> u;

    prox::get_velocity_vector_eigen(bodies.begin(), bodies.end(), u);

    BOOST_CHECK_EQUAL(u(0), 5);
    BOOST_CHECK_EQUAL(u(1), 7);
    BOOST_CHECK_EQUAL(u(2), 0);

    BOOST_CHECK_EQUAL(u(6), 0);
    BOOST_CHECK_EQUAL(u(7), 1);
    BOOST_CHECK_EQUAL(u(8), 0);

    BOOST_CHECK_EQUAL(u(12), 0);
    BOOST_CHECK_EQUAL(u(13), 3);
    BOOST_CHECK_EQUAL(u(14), 0);

    BOOST_CHECK_EQUAL(u(18), 2);
    BOOST_CHECK_EQUAL(u(19), 5);
    BOOST_CHECK_EQUAL(u(20), 0);
    bodies[2].set_velocity(EigenVector3<T>(8.0, 3.0, 0.5));

    prox::get_velocity_vector_eigen(bodies.begin(), bodies.end(), u);

    BOOST_CHECK_EQUAL(u(0), 5);
    BOOST_CHECK_EQUAL(u(1), 7);
    BOOST_CHECK_EQUAL(u(2), 0);

    BOOST_CHECK_EQUAL(u(6), 0);
    BOOST_CHECK_EQUAL(u(7), 1);
    BOOST_CHECK_EQUAL(u(8), 0);

    BOOST_CHECK_EQUAL(u(12), 8);
    BOOST_CHECK_EQUAL(u(13), 3);
    BOOST_CHECK_EQUAL(u(14), 0.5);

    BOOST_CHECK_EQUAL(u(18), 2);
    BOOST_CHECK_EQUAL(u(19), 5);
    BOOST_CHECK_EQUAL(u(20), 0);
}

BOOST_AUTO_TEST_CASE(set_velocity_vector_test)
{
    using T = float;
  // Initialize geometry
    narrow::System<T> narrow;

    size_t gid = narrow.create_geometry();

    narrow::Geometry<T>::sphere_type sphere;
    sphere.radius() = 1.0f;

    narrow::Geometry<T>& moving_geometry = narrow.get_geometry(gid);
    moving_geometry.add_shape(sphere);

  // Initialize objects
    std::vector<prox::RigidBody<T>> bodies;
    bodies.resize(4u);

    bodies[0].set_geometry_idx(gid);
    bodies[1].set_geometry_idx(gid);
    bodies[2].set_geometry_idx(gid);
    bodies[3].set_geometry_idx(gid);

    Eigen::VectorX<T> u;
    u.resize(4 * 6);
    for (uint32_t i = 0; i < 4 * 6; ++i) { u[i] = i % 6; }

    /*    sparse::fill(u(0));// 2009-09-20 Kenny:  error:`fill' is not a member of 'sparse'
    sparse::fill(u(1));// 2009-09-20 Kenny:  error:`fill' is not a member of 'sparse'
    sparse::fill(u(2));// 2009-09-20 Kenny:  error:`fill' is not a member of 'sparse'
    sparse::fill(u(3));// 2009-09-20 Kenny:  error:`fill' is not a member of 'sparse'*/

    prox::set_velocity_vector_eigen(bodies.begin(), bodies.end(), u);

    BOOST_CHECK_EQUAL(bodies[0].get_velocity()(0), 0);
    BOOST_CHECK_EQUAL(bodies[0].get_velocity()(1), 1);
    BOOST_CHECK_EQUAL(bodies[0].get_velocity()(2), 2);

    BOOST_CHECK_EQUAL(bodies[1].get_velocity()(0), 0);
    BOOST_CHECK_EQUAL(bodies[1].get_velocity()(1), 1);
    BOOST_CHECK_EQUAL(bodies[1].get_velocity()(2), 2);

    BOOST_CHECK_EQUAL(bodies[2].get_velocity()(0), 0);
    BOOST_CHECK_EQUAL(bodies[2].get_velocity()(1), 1);
    BOOST_CHECK_EQUAL(bodies[2].get_velocity()(2), 2);

    BOOST_CHECK_EQUAL(bodies[3].get_velocity()(0), 0);
    BOOST_CHECK_EQUAL(bodies[3].get_velocity()(1), 1);
    BOOST_CHECK_EQUAL(bodies[3].get_velocity()(2), 2);

    /*sparse::fill(u(2), 4);// 2009-09-20 Kenny:  error:`fill' is not a member of 'sparse'
    sparse::fill(u(3), 7);// 2009-09-20 Kenny:  error:`fill' is not a member of 'sparse'*/
    for (uint32_t i = 0; i < 4 * 6; ++i)
    {
        u[i] = i % 6;
        if (i >= 11 && i < 18) { u[i] = (i % 6) + 4; }
        if (i > 17) { u[i] = (i % 6) + 7; }
    }

    prox::set_velocity_vector_eigen(bodies.begin(), bodies.end(), u);

    BOOST_CHECK_EQUAL(bodies[0].get_velocity()(0), 0);
    BOOST_CHECK_EQUAL(bodies[0].get_velocity()(1), 1);
    BOOST_CHECK_EQUAL(bodies[0].get_velocity()(2), 2);

    BOOST_CHECK_EQUAL(bodies[1].get_velocity()(0), 0);
    BOOST_CHECK_EQUAL(bodies[1].get_velocity()(1), 1);
    BOOST_CHECK_EQUAL(bodies[1].get_velocity()(2), 2);

    BOOST_CHECK_EQUAL(bodies[2].get_velocity()(0), 4);
    BOOST_CHECK_EQUAL(bodies[2].get_velocity()(1), 5);
    BOOST_CHECK_EQUAL(bodies[2].get_velocity()(2), 6);

    BOOST_CHECK_EQUAL(bodies[3].get_velocity()(0), 7);
    BOOST_CHECK_EQUAL(bodies[3].get_velocity()(1), 8);
    BOOST_CHECK_EQUAL(bodies[3].get_velocity()(2), 9);
}

BOOST_AUTO_TEST_SUITE_END();
