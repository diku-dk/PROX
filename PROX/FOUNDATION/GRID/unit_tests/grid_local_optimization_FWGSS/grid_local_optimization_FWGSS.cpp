#include "grid_enclosing_indices.h"
#include "grid_grid.h"
#include "grid_iterators.h"
#include "grid_node_position.h"
#include "grid_value_at.h"

#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>
#include <fstream>
#include <filesystem>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

#include <grid_ccd_golden_section_search.hpp>
#include <grid_local_optimization.hpp>
#include "igl/signed_distance.h"
#include "igl/marching_cubes.h"
#include "igl/writeOBJ.h"
//#include <libigl/signed_distance.h>
//#include <libigl/sphere.h>
//#include <libigl/marching_cubes.h>
//#include <libigl/writeOBJ.h>

#include <igl/read_triangle_mesh.h>

BOOST_AUTO_TEST_SUITE(grid);

template <typename T> bool isClose(EigenVector3<T> a, EigenVector3<T> b)
{
    T eps = 1e-4;
    if (std::abs<T>(a.x() - b.x()) < eps && std::abs<T>(a.y() - b.y()) < eps
        && std::abs<T>(a.z() - b.z()) < eps)
    {
        return true;
    }
    return false;
}

BOOST_AUTO_TEST_CASE(grid_local_strategy)
{
    {
        using D = float;
        using T = float;

        namespace fs = std::filesystem;
        fs::path source_dir = fs::path(__FILE__).parent_path();
        fs::path bunnyRelative = "../../../../../bin/resources/objs/box.obj";
        fs::path bunnyFull = source_dir / bunnyRelative;
        fs::path bunnyNormalized = bunnyFull.lexically_normal();

        // Output
        std::cout << "Bunny path: " << bunnyNormalized << std::endl;
        std::string meshFile = bunnyNormalized.string();

        int res = 32;

        //Read mesh
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        if (!igl::read_triangle_mesh(meshFile, V, F))
        {
            std::cerr << "ERROR: Failed to read mesh: " << meshFile << "\n";
        }
        std::cout << "Read mesh: #V = " << V.rows() << "  #F = " << F.rows()
                  << "\n";

        Eigen::RowVector3d minv = V.colwise().minCoeff();
        Eigen::RowVector3d maxv = V.colwise().maxCoeff();
        Eigen::RowVector3d diag = maxv - minv;
        T longest = diag.maxCoeff();
        //10% padding to our bounding box!
        T pad = 0.10 * longest;
        Eigen::Matrix<T, 3, 1> gmin((T)(minv.x() - pad), (T)(minv.y() - pad),
                                    (T)(minv.z() - pad));
        Eigen::Matrix<T, 3, 1> gmax((T)(maxv.x() + pad), (T)(maxv.y() + pad),
                                    (T)(maxv.z() + pad));

        //create grid (res^3 nodes)
        Eigen::Matrix<size_t, 3, 1> nodes((size_t)res, (size_t)res,
                                          (size_t)res);
        grid::Grid<D, T> G;
        G.create(gmin, gmax, nodes);
        const size_t total = G.m_nodes.x() * G.m_nodes.y() * G.m_nodes.z();
        std::cout << "Created grid: " << G.I() << " x " << G.J() << " x "
                  << G.K() << "  (total nodes = " << total << ")\n";

        //Build the query points matrix P (total x 3) in the same linear order used by grid
        Eigen::MatrixXd P((Eigen::Index)total, 3);
        size_t idx_lin = 0;
        for (size_t k = 0; k < G.K(); ++k)
        {
            for (size_t j = 0; j < G.J(); ++j)
            {
                for (size_t i = 0; i < G.I(); ++i)
                {
                    Eigen::Matrix<size_t, 3, 1> idx(i, j, k);
                    Eigen::Matrix<T, 3, 1> p;
                    grid::node_position(G, idx, p);
                    P((Eigen::Index)idx_lin, 0) = p.x();
                    P((Eigen::Index)idx_lin, 1) = p.y();
                    P((Eigen::Index)idx_lin, 2) = p.z();
                    ++idx_lin;
                }
            }
        }

        std::cout
            << "Computing signed distances (libigl::signed_distance)...\n";
        Eigen::VectorXd S; // signed distances
        Eigen::VectorXi I; // indices to closest triangles (unused here)
        Eigen::MatrixXd C; // closest points on triangles
        Eigen::MatrixXd N; // normals used for signing (pseudonormal)

        igl::SignedDistanceType signType
            = igl::SIGNED_DISTANCE_TYPE_WINDING_NUMBER;
        igl::signed_distance(P, V, F, signType, S, I, C, N);
        if ((size_t)S.size() != total)
        {
            std::cerr << "ERROR: libigl returned unexpected S size: "
                      << S.size() << " expected " << total << "\n";
        }

        //Copy into grid storage
        for (size_t s = 0; s < total; ++s)
        {
            G.data()[s] = static_cast<D>(S((Eigen::Index)s));
        }

        //Write CSV (x,y,z,sdf)
        fs::path bunnySDFRelative = "../../../../../bin/output/sdftest/sdf.csv";
        fs::path bunnySDFFull = source_dir / bunnySDFRelative;
        fs::path bunnySDFNormalized = bunnySDFFull.lexically_normal();
        const std::string outCsv = bunnySDFNormalized.string();

        {
            EigenVector3<T> angularVelA(0, 0, 0);
            EigenVector3<T> angularVelB(0, 0, 0);
            EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
            EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
            EigenVector3<T> transA(2.0, 0.0667, 0.0667);
            EigenVector3<T> transB(0, 0, 0);
            EigenVector3<T> linearA(-100.0, 0, 0.0);
            EigenVector3<T> linearB(0, 0, 0);

            RigidBodyInfo<T> rInfo;
            rInfo.A_angularVel = &angularVelA;
            rInfo.B_angularVel = &angularVelB;
            rInfo.A_centerRotation = &rotA;
            rInfo.B_centerRotation = &rotB;
            rInfo.A_centerTranslation = &transA;
            rInfo.B_centerTranslation = &transB;
            rInfo.A_linearVel = &linearA;
            rInfo.B_linearVel = &linearB;

            rInfo.B_sdf = &G;

            rInfo.A_p0 = EigenVector3<T>(1.0, 0.0, 0.0);
            rInfo.A_p1 = EigenVector3<T>(1.0, 0.2, 0.0);
            rInfo.A_p2 = EigenVector3<T>(1.0, 0.0, 0.2);
            T val = grid::FrankWolfeGSS<T>(T(0.0), T(0.01), rInfo);
            BOOST_TEST(std::abs<T>(val - T(0.005)) < 0.0001f);
        }

        {
            EigenVector3<T> angularVelA(0, 0, 0);
            EigenVector3<T> angularVelB(0, 0, 0);
            EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
            EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
            EigenVector3<T> transA(0.8333 + 1.0, 0.8333 + 1.0, 0.8 + 1.0);
            EigenVector3<T> transB(0, 0, 0);
            EigenVector3<T> linearA(-80, -80, -80);
            EigenVector3<T> linearB(0, 0, 0);

            RigidBodyInfo<T> rInfo;
            rInfo.A_angularVel = &angularVelA;
            rInfo.B_angularVel = &angularVelB;
            rInfo.A_centerRotation = &rotA;
            rInfo.B_centerRotation = &rotB;
            rInfo.A_centerTranslation = &transA;
            rInfo.B_centerTranslation = &transB;
            rInfo.A_linearVel = &linearA;
            rInfo.B_linearVel = &linearB;

            rInfo.B_sdf = &G;

            rInfo.A_p0 = EigenVector3<T>(0.8, 0.8, 0.8);
            rInfo.A_p1 = EigenVector3<T>(0.9, 0.8, 0.8);
            rInfo.A_p2 = EigenVector3<T>(0.8, 0.9, 0.8);
            T val = grid::FrankWolfeGSS<T>(T(0.0), T(0.01), rInfo);
            BOOST_TEST(std::abs<T>(val - T(0.00375)) < 0.0001f);
        }

        {
            EigenVector3<T> angularVelA(0, 0, 0);
            EigenVector3<T> angularVelB(0, 0, 0);
            EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
            EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
            EigenVector3<T> transA(0.766666667 + 1.0, 0.066666667 + 1.0,
                                   0.7 + 1.0);
            EigenVector3<T> transB(0, 0, 0);
            EigenVector3<T> linearA(-60, 0, -60);
            EigenVector3<T> linearB(0, 0, 0);

            RigidBodyInfo<T> rInfo;
            rInfo.A_angularVel = &angularVelA;
            rInfo.B_angularVel = &angularVelB;
            rInfo.A_centerRotation = &rotA;
            rInfo.B_centerRotation = &rotB;
            rInfo.A_centerTranslation = &transA;
            rInfo.B_centerTranslation = &transB;
            rInfo.A_linearVel = &linearA;
            rInfo.B_linearVel = &linearB;

            rInfo.B_sdf = &G;

            rInfo.A_p0 = EigenVector3<T>(0.7, 0.0, 0.7);
            rInfo.A_p1 = EigenVector3<T>(0.7, 0.2, 0.7);
            rInfo.A_p2 = EigenVector3<T>(0.9, 0.0, 0.7);
            T val = grid::FrankWolfeGSS<T>(T(0.0), T(0.01), rInfo);
            BOOST_TEST(std::abs<T>(val - T(0.00346)) < 0.0001f);
        }

        {
            EigenVector3<T> angularVelA(0, 0, 0);
            EigenVector3<T> angularVelB(0, 0, 0);
            EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
            EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
            EigenVector3<T> transA(0.6 + 1.0, -0.033333333 + 1.0,
                                   0.066666667 + 1.0);
            EigenVector3<T> transB(0, 0, 0);
            EigenVector3<T> linearA(-50, 10, 0);
            EigenVector3<T> linearB(0, 0, 0);

            RigidBodyInfo<T> rInfo;
            rInfo.A_angularVel = &angularVelA;
            rInfo.B_angularVel = &angularVelB;
            rInfo.A_centerRotation = &rotA;
            rInfo.B_centerRotation = &rotB;
            rInfo.A_centerTranslation = &transA;
            rInfo.B_centerTranslation = &transB;
            rInfo.A_linearVel = &linearA;
            rInfo.B_linearVel = &linearB;

            rInfo.B_sdf = &G;

            rInfo.A_p0 = EigenVector3<T>(0.6, -0.1, 0.0);
            rInfo.A_p1 = EigenVector3<T>(0.6, 0.1, 0.0);
            rInfo.A_p2 = EigenVector3<T>(0.6, -0.1, 0.2);
            T val = grid::FrankWolfeGSS<T>(T(0.0), T(0.01), rInfo);
            BOOST_TEST(std::abs<T>(val - T(0.002)) < 0.0001f);
        }

        {
            EigenVector3<T> angularVelA(0, 0, 0);
            EigenVector3<T> angularVelB(0, 0, 0);
            EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
            EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
            EigenVector3<T> transA(0.8333 + 1.0, 0.8333 + 1.0, 0.8 + 1.0);
            EigenVector3<T> transB(0, 0, 0);
            EigenVector3<T> linearA(0.0, -100.0, 0.0);
            EigenVector3<T> linearB(0, 0, 0);

            RigidBodyInfo<T> rInfo;
            rInfo.A_angularVel = &angularVelA;
            rInfo.B_angularVel = &angularVelB;
            rInfo.A_centerRotation = &rotA;
            rInfo.B_centerRotation = &rotB;
            rInfo.A_centerTranslation = &transA;
            rInfo.B_centerTranslation = &transB;
            rInfo.A_linearVel = &linearA;
            rInfo.B_linearVel = &linearB;

            rInfo.B_sdf = &G;

            rInfo.A_p0 = EigenVector3<T>(0.8, 0.8, 0.8);
            rInfo.A_p1 = EigenVector3<T>(0.9, 0.8, 0.8);
            rInfo.A_p2 = EigenVector3<T>(0.8, 0.9, 0.8);
            T val = grid::FrankWolfeGSS<T>(T(0.0), T(0.01), rInfo);
            BOOST_TEST(std::abs<T>(val - T(0.01)) < 0.0001f);
        }

        {
            EigenVector3<T> angularVelA(0, 0, 0);
            EigenVector3<T> angularVelB(0, 0, 0);
            EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
            EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
            EigenVector3<T> transA(0.733333333 + 1.0, 0.333333333 + 1.0,
                                   0.4 + 1.0);
            EigenVector3<T> transB(0, 0, 0);
            EigenVector3<T> linearA(-40.0, -20.0, -30.0);
            EigenVector3<T> linearB(0, 0, 0);

            RigidBodyInfo<T> rInfo;
            rInfo.A_angularVel = &angularVelA;
            rInfo.B_angularVel = &angularVelB;
            rInfo.A_centerRotation = &rotA;
            rInfo.B_centerRotation = &rotB;
            rInfo.A_centerTranslation = &transA;
            rInfo.B_centerTranslation = &transB;
            rInfo.A_linearVel = &linearA;
            rInfo.B_linearVel = &linearB;

            rInfo.B_sdf = &G;

            rInfo.A_p0 = EigenVector3<T>(0.7, 0.3, 0.4);
            rInfo.A_p1 = EigenVector3<T>(0.8, 0.3, 0.4);
            rInfo.A_p2 = EigenVector3<T>(0.7, 0.4, 0.4);
            T val = grid::FrankWolfeGSS<T>(T(0.0), T(0.01), rInfo);
            BOOST_TEST(std::abs<T>(val - T(0.005)) < 0.0001f);
        }

        {
            EigenVector3<T> angularVelA(0, 0, 0);
            EigenVector3<T> angularVelB(0, 0, 0);
            EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
            EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
            EigenVector3<T> transA(0.066666667 + 1.0, 0.066666667 + 1.0,
                                   1.0 + 1.0);
            EigenVector3<T> transB(0, 0, 0);
            EigenVector3<T> linearA(0.0, 0.0, -100.0);
            EigenVector3<T> linearB(0, 0, 0);

            RigidBodyInfo<T> rInfo;
            rInfo.A_angularVel = &angularVelA;
            rInfo.B_angularVel = &angularVelB;
            rInfo.A_centerRotation = &rotA;
            rInfo.B_centerRotation = &rotB;
            rInfo.A_centerTranslation = &transA;
            rInfo.B_centerTranslation = &transB;
            rInfo.A_linearVel = &linearA;
            rInfo.B_linearVel = &linearB;

            rInfo.B_sdf = &G;

            rInfo.A_p0 = EigenVector3<T>(0.0, 0.0, 1.0);
            rInfo.A_p1 = EigenVector3<T>(0.2, 0.0, 1.0);
            rInfo.A_p2 = EigenVector3<T>(0.0, 0.2, 1.0);
            T val = grid::FrankWolfeGSS<T>(T(0.0), T(0.01), rInfo);
            BOOST_TEST(std::abs<T>(val - T(0.005)) < 0.0001f);
        }

        {
            EigenVector3<T> angularVelA(0, 0, 0);
            EigenVector3<T> angularVelB(0, 0, 0);
            EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
            EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
            EigenVector3<T> transA(0.8333 + 1.0, 0.8333 + 1.0, 0.8 + 1.0);
            EigenVector3<T> transB(0, 0, 0);
            EigenVector3<T> linearA(-50.0, -50.0, -50.0);
            EigenVector3<T> linearB(0, 0, 0);

            RigidBodyInfo<T> rInfo;
            rInfo.A_angularVel = &angularVelA;
            rInfo.B_angularVel = &angularVelB;
            rInfo.A_centerRotation = &rotA;
            rInfo.B_centerRotation = &rotB;
            rInfo.A_centerTranslation = &transA;
            rInfo.B_centerTranslation = &transB;
            rInfo.A_linearVel = &linearA;
            rInfo.B_linearVel = &linearB;

            rInfo.B_sdf = &G;

            rInfo.A_p0 = EigenVector3<T>(0.8, 0.8, 0.8);
            rInfo.A_p1 = EigenVector3<T>(0.9, 0.8, 0.8);
            rInfo.A_p2 = EigenVector3<T>(0.8, 0.9, 0.8);
            T val = grid::FrankWolfeGSS<T>(T(0.0), T(0.01), rInfo);
            BOOST_TEST(std::abs<T>(val - T(0.00618)) < 0.0001f);
        }

        {
            EigenVector3<T> angularVelA(0, 0, 0);
            EigenVector3<T> angularVelB(0, 0, 0);
            EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
            EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
            EigenVector3<T> transA(0.583333333 + 1.0, 0.033333333 + 1.0,
                                   0.55 + 1.0);
            EigenVector3<T> transB(0, 0, 0);
            EigenVector3<T> linearA(-5.0, 0.0, -5.0);
            EigenVector3<T> linearB(0, 0, 0);

            RigidBodyInfo<T> rInfo;
            rInfo.A_angularVel = &angularVelA;
            rInfo.B_angularVel = &angularVelB;
            rInfo.A_centerRotation = &rotA;
            rInfo.B_centerRotation = &rotB;
            rInfo.A_centerTranslation = &transA;
            rInfo.B_centerTranslation = &transB;
            rInfo.A_linearVel = &linearA;
            rInfo.B_linearVel = &linearB;

            rInfo.B_sdf = &G;

            rInfo.A_p0 = EigenVector3<T>(0.55, 0.0, 0.55);
            rInfo.A_p1 = EigenVector3<T>(0.55, 0.1, 0.55);
            rInfo.A_p2 = EigenVector3<T>(0.65, 0.0, 0.55);
            T val = grid::FrankWolfeGSS<T>(T(0.0), T(0.01), rInfo);
            BOOST_TEST(std::abs<T>(val - T(0.01)) < 0.0001f);
        }

        {
            //No Collision
            EigenVector3<T> angularVelA(0, 0, 0);
            EigenVector3<T> angularVelB(0, 0, 0);
            EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
            EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
            EigenVector3<T> transA(0.633333333 + 1.0, -0.033333333 + 1.0,
                                   0.6 + 1.0);
            EigenVector3<T> transB(0, 0, 0);
            EigenVector3<T> linearA(0, -1, 0);
            EigenVector3<T> linearB(0, 0, 0);

            RigidBodyInfo<T> rInfo;
            rInfo.A_angularVel = &angularVelA;
            rInfo.B_angularVel = &angularVelB;
            rInfo.A_centerRotation = &rotA;
            rInfo.B_centerRotation = &rotB;
            rInfo.A_centerTranslation = &transA;
            rInfo.B_centerTranslation = &transB;
            rInfo.A_linearVel = &linearA;
            rInfo.B_linearVel = &linearB;

            rInfo.B_sdf = &G;

            rInfo.A_p0 = EigenVector3<T>(0.6, 0.9, 0.6);
            rInfo.A_p1 = EigenVector3<T>(0.6, 0.8, 0.6);
            rInfo.A_p2 = EigenVector3<T>(0.7, 0.9, 0.6);
            T val = grid::FrankWolfeGSS<T>(T(0.0), T(0.01), rInfo);
            BOOST_TEST(std::abs<T>(val - T(0.01)) < 0.0001f);
        }

        {
            //Collision with Y-Face
            EigenVector3<T> angularVelA(0, 0, 0);
            EigenVector3<T> angularVelB(0, 0, 0);
            EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
            EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
            EigenVector3<T> transA(0.066666667 + 1.0, 1.0 + 1.0,
                                   0.066666667 + 1.0);
            EigenVector3<T> transB(0, 0, 0);
            EigenVector3<T> linearA(0, -100, 0);
            EigenVector3<T> linearB(0, 0, 0);

            RigidBodyInfo<T> rInfo;
            rInfo.A_angularVel = &angularVelA;
            rInfo.B_angularVel = &angularVelB;
            rInfo.A_centerRotation = &rotA;
            rInfo.B_centerRotation = &rotB;
            rInfo.A_centerTranslation = &transA;
            rInfo.B_centerTranslation = &transB;
            rInfo.A_linearVel = &linearA;
            rInfo.B_linearVel = &linearB;

            rInfo.B_sdf = &G;

            rInfo.A_p0 = EigenVector3<T>(0.0, 1.0, 0.0);
            rInfo.A_p1 = EigenVector3<T>(0.2, 1.0, 0.0);
            rInfo.A_p2 = EigenVector3<T>(0.0, 1.0, 0.2);
            T val = grid::FrankWolfeGSS<T>(T(0.0), T(0.01), rInfo);
            BOOST_TEST(std::abs<T>(val - T(0.005)) < 0.0001f);
        }

        {
            //Large Triangle Covering Multiple Faces
            EigenVector3<T> angularVelA(0, 0, 0);
            EigenVector3<T> angularVelB(0, 0, 0);
            EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
            EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
            EigenVector3<T> transA(1.0 + 1.0, 0 + 1.0, 0.066666667 + 1.0);
            EigenVector3<T> transB(0, 0, 0);
            EigenVector3<T> linearA(-100, 0, 0);
            EigenVector3<T> linearB(0, 0, 0);

            RigidBodyInfo<T> rInfo;
            rInfo.A_angularVel = &angularVelA;
            rInfo.B_angularVel = &angularVelB;
            rInfo.A_centerRotation = &rotA;
            rInfo.B_centerRotation = &rotB;
            rInfo.A_centerTranslation = &transA;
            rInfo.B_centerTranslation = &transB;
            rInfo.A_linearVel = &linearA;
            rInfo.B_linearVel = &linearB;

            rInfo.B_sdf = &G;

            rInfo.A_p0 = EigenVector3<T>(1.0, -0.2, 0.0);
            rInfo.A_p1 = EigenVector3<T>(1.0, 0.2, 0.0);
            rInfo.A_p2 = EigenVector3<T>(1.0, 0.0, 0.2);
            T val = grid::FrankWolfeGSS<T>(T(0.0), T(0.01), rInfo);
            BOOST_TEST(std::abs<T>(val - T(0.005)) < 0.0001f);
        }

        {
            //Initial Penetration -- (almost)
            EigenVector3<T> angularVelA(0, 0, 0);
            EigenVector3<T> angularVelB(0, 0, 0);
            EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
            EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
            EigenVector3<T> transA(0.5 + 1.0, 0.066666667 + 1.0, 0.0 + 1.0);
            EigenVector3<T> transB(0, 0, 0);
            EigenVector3<T> linearA(-10, 0, 0);
            EigenVector3<T> linearB(0, 0, 0);

            RigidBodyInfo<T> rInfo;
            rInfo.A_angularVel = &angularVelA;
            rInfo.B_angularVel = &angularVelB;
            rInfo.A_centerRotation = &rotA;
            rInfo.B_centerRotation = &rotB;
            rInfo.A_centerTranslation = &transA;
            rInfo.B_centerTranslation = &transB;
            rInfo.A_linearVel = &linearA;
            rInfo.B_linearVel = &linearB;

            rInfo.B_sdf = &G;

            rInfo.A_p0 = EigenVector3<T>(0.4, 0.0, 0.0);
            rInfo.A_p1 = EigenVector3<T>(0.6, 0.0, 0.0);
            rInfo.A_p2 = EigenVector3<T>(0.5, 0.2, 0.0);
            T val = grid::FrankWolfeGSS<T>(T(0.0), T(0.01), rInfo);
            BOOST_TEST(std::abs<T>(val - T(0.000)) < 0.0001f);
        }

        {
            //Moving Away from Box (I think undefined, should be removed in narrow phase?)
            //The reason I set it to 0 is because our gradient direction di would
            // see us moving in the negative direction, thus it thinks we have to go
            // backwards in time, but we are already at t=0, so backwards in time
            // is not possible so we stop. This is incorrect behaviour, but I am honestly
            // not sure if this is an issue!
            EigenVector3<T> angularVelA(0, 0, 0);
            EigenVector3<T> angularVelB(0, 0, 0);
            EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
            EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
            EigenVector3<T> transA(0.633333333 + 1.0, 0.033333333 + 1.0,
                                   0.6 + 1.0);
            EigenVector3<T> transB(0, 0, 0);
            EigenVector3<T> linearA(10, 0, 10);
            EigenVector3<T> linearB(0, 0, 0);

            RigidBodyInfo<T> rInfo;
            rInfo.A_angularVel = &angularVelA;
            rInfo.B_angularVel = &angularVelB;
            rInfo.A_centerRotation = &rotA;
            rInfo.B_centerRotation = &rotB;
            rInfo.A_centerTranslation = &transA;
            rInfo.B_centerTranslation = &transB;
            rInfo.A_linearVel = &linearA;
            rInfo.B_linearVel = &linearB;

            rInfo.B_sdf = &G;

            rInfo.A_p0 = EigenVector3<T>(0.6, 0.0, 0.6);
            rInfo.A_p1 = EigenVector3<T>(0.6, 0.1, 0.6);
            rInfo.A_p2 = EigenVector3<T>(0.7, 0.0, 0.6);
            T val = grid::FrankWolfeGSS<T>(T(0.0), T(0.01), rInfo);
            BOOST_TEST(std::abs<T>(val - T(0.010)) < 0.0001f);
        }

        {
            //Close to surface
            EigenVector3<T> angularVelA(0, 0, 0);
            EigenVector3<T> angularVelB(0, 0, 0);
            EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
            EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
            EigenVector3<T> transA(0.533333333 + 1.0, 0.066666667 + 1.0,
                                   0.533333333 + 1.0);
            EigenVector3<T> transB(0, 0, 0);
            EigenVector3<T> linearA(0, -10, 0);
            EigenVector3<T> linearB(0, 0, 0);

            RigidBodyInfo<T> rInfo;
            rInfo.A_angularVel = &angularVelA;
            rInfo.B_angularVel = &angularVelB;
            rInfo.A_centerRotation = &rotA;
            rInfo.B_centerRotation = &rotB;
            rInfo.A_centerTranslation = &transA;
            rInfo.B_centerTranslation = &transB;
            rInfo.A_linearVel = &linearA;
            rInfo.B_linearVel = &linearB;

            rInfo.B_sdf = &G;

            rInfo.A_p0 = EigenVector3<T>(0.5, 0.6, 0.5);
            rInfo.A_p1 = EigenVector3<T>(0.5, 0.6, 0.6);
            rInfo.A_p2 = EigenVector3<T>(0.6, 0.6, 0.5);
            T val = grid::FrankWolfeGSS<T>(T(0.0), T(0.01), rInfo);
            BOOST_TEST(std::abs<T>(val - T(0.01)) < 0.0001f);
        }

        {
            //Very fast-moving "small triangle
            EigenVector3<T> angularVelA(0, 0, 0);
            EigenVector3<T> angularVelB(0, 0, 0);
            EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
            EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
            EigenVector3<T> transA(2.0 + 1.0, 0.003333333 + 1.0,
                                   0.003333333 + 1.0);
            EigenVector3<T> transB(0, 0, 0);
            EigenVector3<T> linearA(-500, 0, 0);
            EigenVector3<T> linearB(0, 0, 0);

            RigidBodyInfo<T> rInfo;
            rInfo.A_angularVel = &angularVelA;
            rInfo.B_angularVel = &angularVelB;
            rInfo.A_centerRotation = &rotA;
            rInfo.B_centerRotation = &rotB;
            rInfo.A_centerTranslation = &transA;
            rInfo.B_centerTranslation = &transB;
            rInfo.A_linearVel = &linearA;
            rInfo.B_linearVel = &linearB;

            rInfo.B_sdf = &G;

            rInfo.A_p0 = EigenVector3<T>(1.0, 0.0, 0.0);
            rInfo.A_p1 = EigenVector3<T>(1.0, 0.1, 0.0);
            rInfo.A_p2 = EigenVector3<T>(1.0, 0.0, 0.1);
            T val = grid::FrankWolfeGSS<T>(T(0.0), T(0.01), rInfo);
            BOOST_TEST(std::abs<T>(val - T(0.001)) < 0.0001f);
        }

        {
            //Very fast-moving "small triangle
            EigenVector3<T> angularVelA(0, 0, 0);
            EigenVector3<T> angularVelB(0, 0, 0);
            EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
            EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
            EigenVector3<T> transA(4.0 + 1.0, 0.003333333 + 1.0,
                                   0.003333333 + 1.0);
            EigenVector3<T> transB(0, 0, 0);
            EigenVector3<T> linearA(-500, 0, 0);
            EigenVector3<T> linearB(0, 0, 0);

            RigidBodyInfo<T> rInfo;
            rInfo.A_angularVel = &angularVelA;
            rInfo.B_angularVel = &angularVelB;
            rInfo.A_centerRotation = &rotA;
            rInfo.B_centerRotation = &rotB;
            rInfo.A_centerTranslation = &transA;
            rInfo.B_centerTranslation = &transB;
            rInfo.A_linearVel = &linearA;
            rInfo.B_linearVel = &linearB;

            rInfo.B_sdf = &G;

            rInfo.A_p0 = EigenVector3<T>(4.0, 0.0, 0.0);
            rInfo.A_p1 = EigenVector3<T>(4.0, 0.01, 0.0);
            rInfo.A_p2 = EigenVector3<T>(4.0, 0.0, 0.01);
            T val = grid::FrankWolfeGSS<T>(T(0.0), T(0.01), rInfo);
            BOOST_TEST(std::abs<T>(val - T(0.007)) < 0.0001f);
        }

        {
            //Very fast-moving "small triangle
            EigenVector3<T> angularVelA(0, 0, 0);
            EigenVector3<T> angularVelB(0, 0, 0);
            EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
            EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
            EigenVector3<T> transA(0.733333333 + 1.0, 0.7 + 1.0,
                                   0.003333333 + 1.0);
            EigenVector3<T> transB(0, 0, 0);
            EigenVector3<T> linearA(-50, -50, 0);
            EigenVector3<T> linearB(0, 0, 0);

            RigidBodyInfo<T> rInfo;
            rInfo.A_angularVel = &angularVelA;
            rInfo.B_angularVel = &angularVelB;
            rInfo.A_centerRotation = &rotA;
            rInfo.B_centerRotation = &rotB;
            rInfo.A_centerTranslation = &transA;
            rInfo.B_centerTranslation = &transB;
            rInfo.A_linearVel = &linearA;
            rInfo.B_linearVel = &linearB;

            rInfo.B_sdf = &G;

            rInfo.A_p0 = EigenVector3<T>(0.7, 0.7, 0.0);
            rInfo.A_p1 = EigenVector3<T>(0.7, 0.7, 0.1);
            rInfo.A_p2 = EigenVector3<T>(0.8, 0.7, 0.0);
            T val = grid::FrankWolfeGSS<T>(T(0.0), T(0.01), rInfo);
            BOOST_TEST(std::abs<T>(val - T(0.00382)) < 0.0001f);
        }
    }
}

BOOST_AUTO_TEST_CASE(grid_local_opt_FWGSS_torus)
{
    {
        using D = float;
        using T = float;

        namespace fs = std::filesystem;
        fs::path source_dir = fs::path(__FILE__).parent_path();
        fs::path bunnyRelative = "../../../../../bin/resources/objs/torus.obj";
        fs::path bunnyFull = source_dir / bunnyRelative;
        fs::path bunnyNormalized = bunnyFull.lexically_normal();

        // Output
        std::cout << "Bunny path: " << bunnyNormalized << std::endl;
        std::string meshFile = bunnyNormalized.string();

        int res = 32;

        //Read mesh
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        if (!igl::read_triangle_mesh(meshFile, V, F))
        {
            std::cerr << "ERROR: Failed to read mesh: " << meshFile << "\n";
        }
        std::cout << "Read mesh: #V = " << V.rows() << "  #F = " << F.rows()
                  << "\n";

        Eigen::RowVector3d minv = V.colwise().minCoeff();
        Eigen::RowVector3d maxv = V.colwise().maxCoeff();
        Eigen::RowVector3d diag = maxv - minv;
        T longest = diag.maxCoeff();
        //10% padding to our bounding box!
        T pad = 0.10 * longest;
        Eigen::Matrix<T, 3, 1> gmin((T)(minv.x() - pad), (T)(minv.y() - pad),
                                    (T)(minv.z() - pad));
        Eigen::Matrix<T, 3, 1> gmax((T)(maxv.x() + pad), (T)(maxv.y() + pad),
                                    (T)(maxv.z() + pad));

        //create grid (res^3 nodes)
        Eigen::Matrix<size_t, 3, 1> nodes((size_t)res, (size_t)res,
                                          (size_t)res);
        grid::Grid<D, T> G;
        G.create(gmin, gmax, nodes);
        const size_t total = G.m_nodes.x() * G.m_nodes.y() * G.m_nodes.z();
        std::cout << "Created grid: " << G.I() << " x " << G.J() << " x "
                  << G.K() << "  (total nodes = " << total << ")\n";

        //Build the query points matrix P (total x 3) in the same linear order used by grid
        Eigen::MatrixXd P((Eigen::Index)total, 3);
        size_t idx_lin = 0;
        for (size_t k = 0; k < G.K(); ++k)
        {
            for (size_t j = 0; j < G.J(); ++j)
            {
                for (size_t i = 0; i < G.I(); ++i)
                {
                    Eigen::Matrix<size_t, 3, 1> idx(i, j, k);
                    Eigen::Matrix<T, 3, 1> p;
                    grid::node_position(G, idx, p);
                    P((Eigen::Index)idx_lin, 0) = p.x();
                    P((Eigen::Index)idx_lin, 1) = p.y();
                    P((Eigen::Index)idx_lin, 2) = p.z();
                    ++idx_lin;
                }
            }
        }

        std::cout
            << "Computing signed distances (libigl::signed_distance)...\n";
        Eigen::VectorXd S; // signed distances
        Eigen::VectorXi I; // indices to closest triangles (unused here)
        Eigen::MatrixXd C; // closest points on triangles
        Eigen::MatrixXd N; // normals used for signing (pseudonormal)

        igl::SignedDistanceType signType
            = igl::SIGNED_DISTANCE_TYPE_WINDING_NUMBER;
        igl::signed_distance(P, V, F, signType, S, I, C, N);
        if ((size_t)S.size() != total)
        {
            std::cerr << "ERROR: libigl returned unexpected S size: "
                      << S.size() << " expected " << total << "\n";
        }

        //Copy into grid storage
        for (size_t s = 0; s < total; ++s)
        {
            G.data()[s] = static_cast<D>(S((Eigen::Index)s));
        }

        //Write CSV (x,y,z,sdf)
        fs::path bunnySDFRelative = "../../../../../bin/output/sdftest/sdf.csv";
        fs::path bunnySDFFull = source_dir / bunnySDFRelative;
        fs::path bunnySDFNormalized = bunnySDFFull.lexically_normal();
        const std::string outCsv = bunnySDFNormalized.string();

        {
            //From top
            EigenVector3<T> angularVelA(0, 0, 0);
            EigenVector3<T> angularVelB(0, 0, 0);
            EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
            EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
            EigenVector3<T> transA(0.566666667 + 1.0, 1.0 + 1.0, 0.0667 + 1.0);
            EigenVector3<T> transB(0, 0, 0);
            EigenVector3<T> linearA(0.0, -100, 0.0);
            EigenVector3<T> linearB(0, 0, 0);

            RigidBodyInfo<T> rInfo;
            rInfo.A_angularVel = &angularVelA;
            rInfo.B_angularVel = &angularVelB;
            rInfo.A_centerRotation = &rotA;
            rInfo.B_centerRotation = &rotB;
            rInfo.A_centerTranslation = &transA;
            rInfo.B_centerTranslation = &transB;
            rInfo.A_linearVel = &linearA;
            rInfo.B_linearVel = &linearB;

            rInfo.B_sdf = &G;

            rInfo.A_p0 = EigenVector3<T>(0.5, 1.0, 0.0);
            rInfo.A_p1 = EigenVector3<T>(0.7, 1.0, 0.0);
            rInfo.A_p2 = EigenVector3<T>(0.5, 1.0, 0.2);

            T val = grid::FrankWolfeGSS<T>(T(0.0), T(0.01), rInfo);

            /*std::cerr << "VALUE: "
                      << grid::valueAtProjection(
                             G, EigenVector3<T>(0.5, 0.0106432, 0),
                             *(rInfo.B_centerTranslation),
                             *(rInfo.B_centerRotation));*/
            BOOST_TEST(std::abs<T>(val - T(0.0075)) < 0.0001f);
        }

        {
            //From side
            EigenVector3<T> angularVelA(0, 0, 0);
            EigenVector3<T> angularVelB(0, 0, 0);
            EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
            EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
            EigenVector3<T> transA(1.066666667 + 1.0, 0.066666667 + 1.0,
                                   0.0 + 1.0);
            EigenVector3<T> transB(0, 0, 0);
            EigenVector3<T> linearA(-100, 0.0, 0.0);
            EigenVector3<T> linearB(0, 0, 0);

            RigidBodyInfo<T> rInfo;
            rInfo.A_angularVel = &angularVelA;
            rInfo.B_angularVel = &angularVelB;
            rInfo.A_centerRotation = &rotA;
            rInfo.B_centerRotation = &rotB;
            rInfo.A_centerTranslation = &transA;
            rInfo.B_centerTranslation = &transB;
            rInfo.A_linearVel = &linearA;
            rInfo.B_linearVel = &linearB;

            rInfo.B_sdf = &G;

            rInfo.A_p0 = EigenVector3<T>(1.0, 0.0, 0.0);
            rInfo.A_p1 = EigenVector3<T>(1.2, 0.0, 0.0);
            rInfo.A_p2 = EigenVector3<T>(1.0, 0.2, 0.0);

            T val = grid::FrankWolfeGSS<T>(T(0.0), T(0.01), rInfo);

            BOOST_TEST(std::abs<T>(val - T(0.0025)) < 0.0001f);
        }

        {
            //Inside torus to outside
            EigenVector3<T> angularVelA(0, 0, 0);
            EigenVector3<T> angularVelB(0, 0, 0);
            EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
            EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
            EigenVector3<T> transA(0.366666667 + 1.0, 0.066666667 + 1.0,
                                   0.0 + 1.0);
            EigenVector3<T> transB(0, 0, 0);
            EigenVector3<T> linearA(-50, 0.0, 0.0);
            EigenVector3<T> linearB(0, 0, 0);

            RigidBodyInfo<T> rInfo;
            rInfo.A_angularVel = &angularVelA;
            rInfo.B_angularVel = &angularVelB;
            rInfo.A_centerRotation = &rotA;
            rInfo.B_centerRotation = &rotB;
            rInfo.A_centerTranslation = &transA;
            rInfo.B_centerTranslation = &transB;
            rInfo.A_linearVel = &linearA;
            rInfo.B_linearVel = &linearB;

            rInfo.B_sdf = &G;

            rInfo.A_p0 = EigenVector3<T>(0.3, 0.0, 0.0);
            rInfo.A_p1 = EigenVector3<T>(0.5, 0.0, 0.0);
            rInfo.A_p2 = EigenVector3<T>(0.3, 0.2, 0.0);

            T val = grid::FrankWolfeGSS<T>(T(0.0), T(0.01), rInfo);

            BOOST_TEST(std::abs<T>(val - T(0.000)) < 0.0001f);
        }

        {
            //Bottom up
            EigenVector3<T> angularVelA(0, 0, 0);
            EigenVector3<T> angularVelB(0, 0, 0);
            EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
            EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
            EigenVector3<T> transA(0.566666667 + 1.0, -1.0 + 1.0,
                                   0.066666667 + 1.0);
            EigenVector3<T> transB(0, 0, 0);
            EigenVector3<T> linearA(0.0, 100.0, 0.0);
            EigenVector3<T> linearB(0, 0, 0);

            RigidBodyInfo<T> rInfo;
            rInfo.A_angularVel = &angularVelA;
            rInfo.B_angularVel = &angularVelB;
            rInfo.A_centerRotation = &rotA;
            rInfo.B_centerRotation = &rotB;
            rInfo.A_centerTranslation = &transA;
            rInfo.B_centerTranslation = &transB;
            rInfo.A_linearVel = &linearA;
            rInfo.B_linearVel = &linearB;

            rInfo.B_sdf = &G;

            rInfo.A_p0 = EigenVector3<T>(0.5, -1.0, 0.0);
            rInfo.A_p1 = EigenVector3<T>(0.7, -1.0, 0.0);
            rInfo.A_p2 = EigenVector3<T>(0.5, -1.0, 0.2);

            T val = grid::FrankWolfeGSS<T>(T(0.0), T(0.01), rInfo);

            BOOST_TEST(std::abs<T>(val - T(0.0075)) < 0.0001f);
        }

        {
            //Surface hit
            EigenVector3<T> angularVelA(0, 0, 0);
            EigenVector3<T> angularVelB(0, 0, 0);
            EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
            EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
            EigenVector3<T> transA(0.866666667 + 1.0, 0.066666667 + 1.0,
                                   0.0 + 1.0);
            EigenVector3<T> transB(0, 0, 0);
            EigenVector3<T> linearA(-50, 0.0, 0.0);
            EigenVector3<T> linearB(0, 0, 0);

            RigidBodyInfo<T> rInfo;
            rInfo.A_angularVel = &angularVelA;
            rInfo.B_angularVel = &angularVelB;
            rInfo.A_centerRotation = &rotA;
            rInfo.B_centerRotation = &rotB;
            rInfo.A_centerTranslation = &transA;
            rInfo.B_centerTranslation = &transB;
            rInfo.A_linearVel = &linearA;
            rInfo.B_linearVel = &linearB;

            rInfo.B_sdf = &G;

            rInfo.A_p0 = EigenVector3<T>(0.8, 0.0, 0.0);
            rInfo.A_p1 = EigenVector3<T>(1.0, 0.0, 0.0);
            rInfo.A_p2 = EigenVector3<T>(0.8, 0.2, 0.0);

            T val = grid::FrankWolfeGSS<T>(T(0.0), T(0.01), rInfo);

            BOOST_TEST(std::abs<T>(val - T(0.001)) < 0.0001f);
        }

        {
            //No collision
            EigenVector3<T> angularVelA(0, 0, 0);
            EigenVector3<T> angularVelB(0, 0, 0);
            EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
            EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
            EigenVector3<T> transA(0.066666667 + 1.0, 1.0 + 1.0,
                                   0.066666667 + 1.0);
            EigenVector3<T> transB(0, 0, 0);
            EigenVector3<T> linearA(100, 0.0, 0.0);
            EigenVector3<T> linearB(0, 0, 0);

            RigidBodyInfo<T> rInfo;
            rInfo.A_angularVel = &angularVelA;
            rInfo.B_angularVel = &angularVelB;
            rInfo.A_centerRotation = &rotA;
            rInfo.B_centerRotation = &rotB;
            rInfo.A_centerTranslation = &transA;
            rInfo.B_centerTranslation = &transB;
            rInfo.A_linearVel = &linearA;
            rInfo.B_linearVel = &linearB;

            rInfo.B_sdf = &G;

            rInfo.A_p0 = EigenVector3<T>(0.0, 1.0, 0.0);
            rInfo.A_p1 = EigenVector3<T>(0.2, 1.0, 0.0);
            rInfo.A_p2 = EigenVector3<T>(0.0, 1.0, 0.2);

            T val = grid::FrankWolfeGSS<T>(T(0.0), T(0.01), rInfo);

            BOOST_TEST(std::abs<T>(val - T(0.01)) < 0.0001f);
        }

        {
            //Front
            EigenVector3<T> angularVelA(0, 0, 0);
            EigenVector3<T> angularVelB(0, 0, 0);
            EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
            EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
            EigenVector3<T> transA(0.066666667 + 1.0, 0.066666667 + 1.0,
                                   1.0 + 1.0);
            EigenVector3<T> transB(0, 0, 0);
            EigenVector3<T> linearA(0.0, 0.0, -100.0);
            EigenVector3<T> linearB(0, 0, 0);

            RigidBodyInfo<T> rInfo;
            rInfo.A_angularVel = &angularVelA;
            rInfo.B_angularVel = &angularVelB;
            rInfo.A_centerRotation = &rotA;
            rInfo.B_centerRotation = &rotB;
            rInfo.A_centerTranslation = &transA;
            rInfo.B_centerTranslation = &transB;
            rInfo.A_linearVel = &linearA;
            rInfo.B_linearVel = &linearB;

            rInfo.B_sdf = &G;

            rInfo.A_p0 = EigenVector3<T>(0.0, 0.0, 1.0);
            rInfo.A_p1 = EigenVector3<T>(0.2, 0.0, 1.0);
            rInfo.A_p2 = EigenVector3<T>(0.0, 0.2, 1.0);

            T val = grid::FrankWolfeGSS<T>(T(0.0), T(0.01), rInfo);

            BOOST_TEST(std::abs<T>(val - T(0.0025)) < 0.0001f);
        }

        {
            //Diagonal approach
            EigenVector3<T> angularVelA(0, 0, 0);
            EigenVector3<T> angularVelB(0, 0, 0);
            EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
            EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
            EigenVector3<T> transA(1.066666667 + 1.0, 0.566666667 + 1.0,
                                   1.0 + 1.0);
            EigenVector3<T> transB(0, 0, 0);
            EigenVector3<T> linearA(-100.0, -100.0, -100.0);
            EigenVector3<T> linearB(0, 0, 0);

            RigidBodyInfo<T> rInfo;
            rInfo.A_angularVel = &angularVelA;
            rInfo.B_angularVel = &angularVelB;
            rInfo.A_centerRotation = &rotA;
            rInfo.B_centerRotation = &rotB;
            rInfo.A_centerTranslation = &transA;
            rInfo.B_centerTranslation = &transB;
            rInfo.A_linearVel = &linearA;
            rInfo.B_linearVel = &linearB;

            rInfo.B_sdf = &G;

            rInfo.A_p0 = EigenVector3<T>(1.0, 0.5, 1.0);
            rInfo.A_p1 = EigenVector3<T>(1.2, 0.5, 1.0);
            rInfo.A_p2 = EigenVector3<T>(1.0, 0.7, 1.0);

            T val = grid::FrankWolfeGSS<T>(T(0.0), T(0.01), rInfo);

            BOOST_TEST(std::abs<T>(val - T(0.00471)) < 0.0001f);
        }

        {
            //Side
            EigenVector3<T> angularVelA(0, 0, 0);
            EigenVector3<T> angularVelB(0, 0, 0);
            EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
            EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
            EigenVector3<T> transA(1.066666667 + 1.0, 0.8 + 1.0, 0.0 + 1.0);
            EigenVector3<T> transB(0, 0, 0);
            EigenVector3<T> linearA(-100.0, 0.0, 0.0);
            EigenVector3<T> linearB(0, 0, 0);

            RigidBodyInfo<T> rInfo;
            rInfo.A_angularVel = &angularVelA;
            rInfo.B_angularVel = &angularVelB;
            rInfo.A_centerRotation = &rotA;
            rInfo.B_centerRotation = &rotB;
            rInfo.A_centerTranslation = &transA;
            rInfo.B_centerTranslation = &transB;
            rInfo.A_linearVel = &linearA;
            rInfo.B_linearVel = &linearB;

            rInfo.B_sdf = &G;

            rInfo.A_p0 = EigenVector3<T>(1.0, 0.2, 0.0);
            rInfo.A_p1 = EigenVector3<T>(1.2, 0.2, 0.0);
            rInfo.A_p2 = EigenVector3<T>(1.0, 0.4, 0.0);

            T val = grid::FrankWolfeGSS<T>(T(0.0), T(0.01), rInfo);

            BOOST_TEST(std::abs<T>(val - T(0.0035)) < 0.0001f);
        }

        {
            //No collision inside hole? Maybe it will orrecctly be incorrect as it should be removed in narrow phase
            EigenVector3<T> angularVelA(0, 0, 0);
            EigenVector3<T> angularVelB(0, 0, 0);
            EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
            EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
            EigenVector3<T> transA(0.266666667 + 1.0, 0.066666667 + 1.0,
                                   0.0 + 1.0);
            EigenVector3<T> transB(0, 0, 0);
            EigenVector3<T> linearA(-100.0, 0.0, 0.0);
            EigenVector3<T> linearB(0, 0, 0);

            RigidBodyInfo<T> rInfo;
            rInfo.A_angularVel = &angularVelA;
            rInfo.B_angularVel = &angularVelB;
            rInfo.A_centerRotation = &rotA;
            rInfo.B_centerRotation = &rotB;
            rInfo.A_centerTranslation = &transA;
            rInfo.B_centerTranslation = &transB;
            rInfo.A_linearVel = &linearA;
            rInfo.B_linearVel = &linearB;

            rInfo.B_sdf = &G;

            rInfo.A_p0 = EigenVector3<T>(0.2, 0.0, 0.0);
            rInfo.A_p1 = EigenVector3<T>(0.4, 0.0, 0.0);
            rInfo.A_p2 = EigenVector3<T>(0.2, 0.2, 0.0);

            T val = grid::FrankWolfeGSS<T>(T(0.0), T(0.01), rInfo);

            BOOST_TEST(std::abs<T>(val - T(0.00)) < 0.0001f);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END();
