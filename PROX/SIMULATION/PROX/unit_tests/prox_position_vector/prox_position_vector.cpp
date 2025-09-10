#include <narrow.h>

#include <prox_rigid_body.h>
#include <prox_get_position_vector.h>
#include <prox_set_position_vector.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(position_vector);

BOOST_AUTO_TEST_CASE(get_position_vector_test)
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
    bodies[0].set_position(EigenVector3<T>(5.0, 7.0, 0.0));

    bodies[1].set_geometry_idx(gid);
    bodies[1].set_position(EigenVector3<T>(0.0, 1.0, 0.0));

    bodies[2].set_geometry_idx(gid);
    bodies[2].set_position(EigenVector3<T>(0.0, 3.0, 0.0));

    bodies[3].set_geometry_idx(gid);
    bodies[3].set_position(EigenVector3<T>(2.0, 5.0, 0.0));

    Eigen::VectorX<T> q;
    prox::get_position_vector_eigen(bodies.begin(), bodies.end(), q);

    BOOST_CHECK_EQUAL(q(0), 5);
    BOOST_CHECK_EQUAL(q(1), 7);
    BOOST_CHECK_EQUAL(q(2), 0);

    BOOST_CHECK_EQUAL(q(7), 0);
    BOOST_CHECK_EQUAL(q(8), 1);
    BOOST_CHECK_EQUAL(q(9), 0);

    BOOST_CHECK_EQUAL(q(14), 0);
    BOOST_CHECK_EQUAL(q(15), 3);
    BOOST_CHECK_EQUAL(q(16), 0);

    BOOST_CHECK_EQUAL(q(21), 2);
    BOOST_CHECK_EQUAL(q(22), 5);
    BOOST_CHECK_EQUAL(q(23), 0);

    bodies[2].set_position(Eigen::Vector3<T>(8.0, 3.0, 0.5));

    prox::get_position_vector_eigen(bodies.begin(), bodies.end(), q);

    BOOST_CHECK_EQUAL(q(0), 5);
    BOOST_CHECK_EQUAL(q(1), 7);
    BOOST_CHECK_EQUAL(q(2), 0);

    BOOST_CHECK_EQUAL(q(7), 0);
    BOOST_CHECK_EQUAL(q(8), 1);
    BOOST_CHECK_EQUAL(q(9), 0);

    BOOST_CHECK_EQUAL(q(14), 8);
    BOOST_CHECK_EQUAL(q(15), 3);
    BOOST_CHECK_EQUAL(q(16), 0.5);

    BOOST_CHECK_EQUAL(q(21), 2);
    BOOST_CHECK_EQUAL(q(22), 5);
    BOOST_CHECK_EQUAL(q(23), 0);
}

BOOST_AUTO_TEST_CASE(set_position_vector_test)
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

    //vector7_type q(4);
    Eigen::VectorX<T> q;
    q.resize(4 * 7);
    for (uint32_t i = 0; i < 4 * 7; ++i) { q[i] = i % 7; }

    /*    sparse::fill(q(0));// 2009-09-20 Kenny:  error:`fill' is not a member of 'sparse'
    sparse::fill(q(1));// 2009-09-20 Kenny:  error:`fill' is not a member of 'sparse'
    sparse::fill(q(2));// 2009-09-20 Kenny:  error:`fill' is not a member of 'sparse'
    sparse::fill(q(3));// 2009-09-20 Kenny:  error:`fill' is not a member of 'sparse'*/

    prox::set_position_vector_eigen(bodies.begin(), bodies.end(), q);

    BOOST_CHECK_EQUAL(bodies[0].get_position()(0), 0);
    BOOST_CHECK_EQUAL(bodies[0].get_position()(1), 1);
    BOOST_CHECK_EQUAL(bodies[0].get_position()(2), 2);

    BOOST_CHECK_EQUAL(bodies[1].get_position()(0), 0);
    BOOST_CHECK_EQUAL(bodies[1].get_position()(1), 1);
    BOOST_CHECK_EQUAL(bodies[1].get_position()(2), 2);

    BOOST_CHECK_EQUAL(bodies[2].get_position()(0), 0);
    BOOST_CHECK_EQUAL(bodies[2].get_position()(1), 1);
    BOOST_CHECK_EQUAL(bodies[2].get_position()(2), 2);

    BOOST_CHECK_EQUAL(bodies[3].get_position()(0), 0);
    BOOST_CHECK_EQUAL(bodies[3].get_position()(1), 1);
    BOOST_CHECK_EQUAL(bodies[3].get_position()(2), 2);

    /*sparse::fill(q(2), 4);// 2009-09-20 Kenny:  error:`fill' is not a member of 'sparse'
    sparse::fill(q(3), 7);// 2009-09-20 Kenny:  error:`fill' is not a member of 'sparse'*/
    for (uint32_t i = 0; i < 4 * 7; ++i)
    {
        q[i] = i % 7;
        if (i >= 12 && i < 18) { q[i] = (i % 7) + 4; }
        if (i > 18) { q[i] = (i % 7) + 7; }
    }

    prox::set_position_vector_eigen(bodies.begin(), bodies.end(), q);

    BOOST_CHECK_EQUAL(bodies[0].get_position()(0), 0);
    BOOST_CHECK_EQUAL(bodies[0].get_position()(1), 1);
    BOOST_CHECK_EQUAL(bodies[0].get_position()(2), 2);

    BOOST_CHECK_EQUAL(bodies[1].get_position()(0), 0);
    BOOST_CHECK_EQUAL(bodies[1].get_position()(1), 1);
    BOOST_CHECK_EQUAL(bodies[1].get_position()(2), 2);

    BOOST_CHECK_EQUAL(bodies[2].get_position()(0), 4);
    BOOST_CHECK_EQUAL(bodies[2].get_position()(1), 5);
    BOOST_CHECK_EQUAL(bodies[2].get_position()(2), 6);

    BOOST_CHECK_EQUAL(bodies[3].get_position()(0), 7);
    BOOST_CHECK_EQUAL(bodies[3].get_position()(1), 8);
    BOOST_CHECK_EQUAL(bodies[3].get_position()(2), 9);
}

BOOST_AUTO_TEST_SUITE_END();
