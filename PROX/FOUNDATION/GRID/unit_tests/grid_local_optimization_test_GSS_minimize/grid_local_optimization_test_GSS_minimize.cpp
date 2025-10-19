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
        using D = double;
        using T = double;

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
            EigenVector3<T> transA(0, 0, 0);
            EigenVector3<T> transB(0, 0, 0);
            EigenVector3<T> linearA(0, 0, -3.0);
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

            rInfo.A_p0 = EigenVector3<T>(0.0, 0.0, 2.0);
            rInfo.A_p1 = EigenVector3<T>(1.0, 0.0, 2.0);
            rInfo.A_p2 = EigenVector3<T>(0.0, 1.0, 2.0);
            DistanceAtTimeParams<T> params{.u = T(1),
                                           .v = T(0),
                                           .w = T(0),
                                           .grid = &G,
                                           .p0 = rInfo.A_p0,
                                           .p1 = rInfo.A_p1,
                                           .p2 = rInfo.A_p2};
            std::cerr << "TEST1.1: "
                      << grid::GSSMinimize_WHAT_MODIFIED(
                             0.0, 1.0, SignedDistanceAtTime<T>, params, rInfo)
                      << "\n";

            rInfo.A_p0 = EigenVector3<T>(2.0, 0.0, 0.0);
            rInfo.A_p1 = EigenVector3<T>(2.0, 1.0, 0.0);
            rInfo.A_p2 = EigenVector3<T>(2.0, 0.0, 1.0);
            EigenVector3<T> linearA2(1, 0, 0.0);
            rInfo.A_linearVel = &linearA2;
            DistanceAtTimeParams<T> params2{.u = T(0.5),
                                            .v = T(0.5),
                                            .w = T(0),
                                            .grid = &G,
                                            .p0 = rInfo.A_p0,
                                            .p1 = rInfo.A_p1,
                                            .p2 = rInfo.A_p2};
            std::cerr << "TEST1.1: "
                      << grid::GSSMinimize_WHAT_MODIFIED(
                             0.0, 1.0, SignedDistanceAtTime<T>, params2, rInfo)
                      << "\n";

            rInfo.A_p0 = EigenVector3<T>(2.0, 0.0, 0.0);
            rInfo.A_p1 = EigenVector3<T>(2.0, 1.0, 0.0);
            rInfo.A_p2 = EigenVector3<T>(2.0, 0.0, 1.0);
            EigenVector3<T> linearA3(-1, 0, 0.0);
            rInfo.A_linearVel = &linearA3;
            DistanceAtTimeParams<T> params3{.u = T(0.5),
                                            .v = T(0.5),
                                            .w = T(0),
                                            .grid = &G,
                                            .p0 = rInfo.A_p0,
                                            .p1 = rInfo.A_p1,
                                            .p2 = rInfo.A_p2};
            std::cerr << "TEST1.2: "
                      << grid::GSSMinimize_WHAT_MODIFIED(
                             0.0, 1.0, SignedDistanceAtTime<T>, params3, rInfo)
                      << "\n";

            rInfo.A_p0 = EigenVector3<T>(-1.0, 0.0, 0.0);
            rInfo.A_p1 = EigenVector3<T>(-1.0, 1.0, 0.0);
            rInfo.A_p2 = EigenVector3<T>(-1.0, 0.0, 1.0);
            EigenVector3<T> linearA4(2, 0, 0.0);
            rInfo.A_linearVel = &linearA4;
            DistanceAtTimeParams<T> params4{.u = T(0.5),
                                            .v = T(0.5),
                                            .w = T(0),
                                            .grid = &G,
                                            .p0 = rInfo.A_p0,
                                            .p1 = rInfo.A_p1,
                                            .p2 = rInfo.A_p2};
            std::cerr << "TEST1.3: "
                      << grid::GSSMinimize_WHAT_MODIFIED(
                             0.0, 1.0, SignedDistanceAtTime<T>, params4, rInfo)
                      << "\n";

            rInfo.A_p0 = EigenVector3<T>(0.0, 2.0, 0.0);
            rInfo.A_p1 = EigenVector3<T>(0.5, 2.0, 0.0);
            rInfo.A_p2 = EigenVector3<T>(0.0, 2.0, 0.5);
            EigenVector3<T> linearA5(0.0, -3.0, 0.0);
            rInfo.A_linearVel = &linearA5;
            DistanceAtTimeParams<T> params5{.u = T(0.5),
                                            .v = T(0.5),
                                            .w = T(0),
                                            .grid = &G,
                                            .p0 = rInfo.A_p0,
                                            .p1 = rInfo.A_p1,
                                            .p2 = rInfo.A_p2};
            std::cerr << "TEST1.4: "
                      << grid::GSSMinimize_WHAT_MODIFIED(
                             0.0, 1.0, SignedDistanceAtTime<T>, params5, rInfo)
                      << "\n";
            /*            std::cerr << "TEST1.4: "
                      << grid::GSSMinimize_WHAT(
                             0.0, 1.0, SignedDistanceAtTime<T>, params5, rInfo)
                      << "\n";
            std::cerr << "TEST1.4: "
                      << grid::GSSMinimize(0.0, 1.0, SignedDistanceAtTime<T>,
                                           params5, rInfo)
                      << "\n";*/
        }

        {
            EigenVector3<T> angularVelA(0, 0, 0);
            EigenVector3<T> angularVelB(0, 0, 0);
            EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
            EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
            EigenVector3<T> transA(0, 0, 0);
            EigenVector3<T> transB(0, 0, 0);
            EigenVector3<T> linearA(1, 0, 0.0);
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

            rInfo.A_p0 = EigenVector3<T>(0.5, 0.0, 0.0);
            rInfo.A_p1 = EigenVector3<T>(1.0, 0.0, 0.0);
            rInfo.A_p2 = EigenVector3<T>(0.5, 1.0, 0.0);
            DistanceAtTimeParams<T> params{.u = T(1),
                                           .v = T(0),
                                           .w = T(0),
                                           .grid = &G,
                                           .p0 = rInfo.A_p0,
                                           .p1 = rInfo.A_p1,
                                           .p2 = rInfo.A_p2};
            std::cerr << "TEST2.1: "
                      << grid::GSSMinimize_WHAT_MODIFIED(
                             0.0, 1.0, UnsignedDistanceAtTime<T>, params, rInfo)
                      << "\n";

            rInfo.A_p0 = EigenVector3<T>(0.0, 0.0, 0.0);
            rInfo.A_p1 = EigenVector3<T>(1.0, 0.0, 0.0);
            rInfo.A_p2 = EigenVector3<T>(0.0, 1.0, 0.0);
            EigenVector3<T> linearA2(-1, 0, 0.0);
            rInfo.A_linearVel = &linearA2;
            DistanceAtTimeParams<T> params2{.u = T(0.5),
                                            .v = T(0.5),
                                            .w = T(0),
                                            .grid = &G,
                                            .p0 = rInfo.A_p0,
                                            .p1 = rInfo.A_p1,
                                            .p2 = rInfo.A_p2};
            std::cerr << "TEST2.2: "
                      << grid::GSSMinimize_WHAT_MODIFIED(
                             0.0, 1.0, UnsignedDistanceAtTime<T>, params2,
                             rInfo)
                      << "\n";

            rInfo.A_p0 = EigenVector3<T>(0.0, 0.0, 0.0);
            rInfo.A_p1 = EigenVector3<T>(0.0, 1.5, 0.0);
            rInfo.A_p2 = EigenVector3<T>(0.0, 0.0, 1.0);
            EigenVector3<T> linearA3(0.0, 0.0, 0.5);
            rInfo.A_linearVel = &linearA3;
            DistanceAtTimeParams<T> params3{.u = T(0.5),
                                            .v = T(0.5),
                                            .w = T(0),
                                            .grid = &G,
                                            .p0 = rInfo.A_p0,
                                            .p1 = rInfo.A_p1,
                                            .p2 = rInfo.A_p2};
            std::cerr << "TEST2.3: "
                      << grid::GSSMinimize_WHAT_MODIFIED(
                             0.0, 1.0, UnsignedDistanceAtTime<T>, params3,
                             rInfo)
                      << "\n";

            rInfo.A_p0 = EigenVector3<T>(0.0, 0.5, 0.0);
            rInfo.A_p1 = EigenVector3<T>(1.0, 0.5, 0.0);
            rInfo.A_p2 = EigenVector3<T>(0.0, 0.5, 1.0);
            EigenVector3<T> linearA4(1, 0, 0.0);
            rInfo.A_linearVel = &linearA4;
            DistanceAtTimeParams<T> params4{.u = T(0.5),
                                            .v = T(0.5),
                                            .w = T(0),
                                            .grid = &G,
                                            .p0 = rInfo.A_p0,
                                            .p1 = rInfo.A_p1,
                                            .p2 = rInfo.A_p2};
            std::cerr << "TEST2.4: "
                      << grid::GSSMinimize_WHAT_MODIFIED(
                             0.0, 1.0, UnsignedDistanceAtTime<T>, params4,
                             rInfo)
                      << "\n";

            rInfo.A_p0 = EigenVector3<T>(0.0, 0.0, 0.0);
            rInfo.A_p1 = EigenVector3<T>(0.0, 0.0, 1.0);
            rInfo.A_p2 = EigenVector3<T>(0.0, 1.0, 0.0);
            EigenVector3<T> linearA5(1.0, 0.0, 0.0);
            rInfo.A_linearVel = &linearA5;
            DistanceAtTimeParams<T> params5{.u = T(1.0),
                                            .v = T(0.0),
                                            .w = T(0.0),
                                            .grid = &G,
                                            .p0 = rInfo.A_p0,
                                            .p1 = rInfo.A_p1,
                                            .p2 = rInfo.A_p2};
            std::cerr << "TEST2.5: "
                      << grid::GSSMinimize_WHAT_MODIFIED(
                             0.0, 2.0, UnsignedDistanceAtTime<T>, params5,
                             rInfo)
                      << "\n";
            /*            std::cerr << "TEST2.5: "
                      << grid::GSSMinimize_WHAT(0.0, 2.0,
                                                UnsignedDistanceAtTime<T>,
                                                params5, rInfo)
                      << "\n";
            std::cerr << "TEST2.5: "
                      << grid::GSSMinimize(0.0, 2.0, UnsignedDistanceAtTime<T>,
                                           params5, rInfo)
                      << "\n";*/
        }

        {
            EigenVector3<T> angularVelA(0, 0, 0);
            EigenVector3<T> angularVelB(0, 0, 0);
            EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
            EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
            EigenVector3<T> transA(0, 0, 0);
            EigenVector3<T> transB(0, 0, 0);
            EigenVector3<T> linearA(0, 0, -3.0);
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

            rInfo.A_p0 = EigenVector3<T>(0.0, 0.0, 2.0);
            rInfo.A_p1 = EigenVector3<T>(1.0, 0.0, 2.0);
            rInfo.A_p2 = EigenVector3<T>(0.0, 1.0, 2.0);
            DistanceAtTimeParams<T> params{.u = T(1),
                                           .v = T(0),
                                           .w = T(0),
                                           .grid = &G,
                                           .p0 = rInfo.A_p0,
                                           .p1 = rInfo.A_p1,
                                           .p2 = rInfo.A_p2};
            std::cerr << "TEST1.1: "
                      << grid::GSSMinimize_WHAT_MODIFIED(
                             0.0, 1.0, SignedDistanceAtTime<T>, params, rInfo)
                      << "\n";

            rInfo.A_p0 = EigenVector3<T>(2.0, 0.0, 0.0);
            rInfo.A_p1 = EigenVector3<T>(2.0, 1.0, 0.0);
            rInfo.A_p2 = EigenVector3<T>(2.0, 0.0, 1.0);
            EigenVector3<T> linearA2(1, 0, 0.0);
            rInfo.A_linearVel = &linearA2;
            DistanceAtTimeParams<T> params2{.u = T(0.5),
                                            .v = T(0.5),
                                            .w = T(0),
                                            .grid = &G,
                                            .p0 = rInfo.A_p0,
                                            .p1 = rInfo.A_p1,
                                            .p2 = rInfo.A_p2};
            std::cerr << "TEST1.1: "
                      << grid::GSSMinimize_WHAT_MODIFIED(
                             0.0, 1.0, SignedDistanceAtTime<T>, params2, rInfo)
                      << "\n";

            rInfo.A_p0 = EigenVector3<T>(2.0, 0.0, 0.0);
            rInfo.A_p1 = EigenVector3<T>(2.0, 1.0, 0.0);
            rInfo.A_p2 = EigenVector3<T>(2.0, 0.0, 1.0);
            EigenVector3<T> linearA3(-1, 0, 0.0);
            rInfo.A_linearVel = &linearA3;
            DistanceAtTimeParams<T> params3{.u = T(0.5),
                                            .v = T(0.5),
                                            .w = T(0),
                                            .grid = &G,
                                            .p0 = rInfo.A_p0,
                                            .p1 = rInfo.A_p1,
                                            .p2 = rInfo.A_p2};
            std::cerr << "TEST1.2: "
                      << grid::GSSMinimize_WHAT_MODIFIED(
                             0.0, 1.0, SignedDistanceAtTime<T>, params3, rInfo)
                      << "\n";

            rInfo.A_p0 = EigenVector3<T>(-1.0, 0.0, 0.0);
            rInfo.A_p1 = EigenVector3<T>(-1.0, 1.0, 0.0);
            rInfo.A_p2 = EigenVector3<T>(-1.0, 0.0, 1.0);
            EigenVector3<T> linearA4(2, 0, 0.0);
            rInfo.A_linearVel = &linearA4;
            DistanceAtTimeParams<T> params4{.u = T(0.5),
                                            .v = T(0.5),
                                            .w = T(0),
                                            .grid = &G,
                                            .p0 = rInfo.A_p0,
                                            .p1 = rInfo.A_p1,
                                            .p2 = rInfo.A_p2};
            std::cerr << "TEST1.3: "
                      << grid::GSSMinimize_WHAT_MODIFIED(
                             0.0, 1.0, SignedDistanceAtTime<T>, params4, rInfo)
                      << "\n";

            rInfo.A_p0 = EigenVector3<T>(0.0, 2.0, 0.0);
            rInfo.A_p1 = EigenVector3<T>(0.5, 2.0, 0.0);
            rInfo.A_p2 = EigenVector3<T>(0.0, 2.0, 0.5);
            EigenVector3<T> linearA5(0.0, -3.0, 0.0);
            rInfo.A_linearVel = &linearA5;
            DistanceAtTimeParams<T> params5{.u = T(0.5),
                                            .v = T(0.5),
                                            .w = T(0),
                                            .grid = &G,
                                            .p0 = rInfo.A_p0,
                                            .p1 = rInfo.A_p1,
                                            .p2 = rInfo.A_p2};
            std::cerr << "TEST1.4: "
                      << grid::GSSMinimize_WHAT_MODIFIED(
                             0.0, 1.0, SignedDistanceAtTime<T>, params5, rInfo)
                      << "\n";
            /*            std::cerr << "TEST1.4: "
                      << grid::GSSMinimize_WHAT(
                             0.0, 1.0, SignedDistanceAtTime<T>, params5, rInfo)
                      << "\n";
            std::cerr << "TEST1.4: "
                      << grid::GSSMinimize(0.0, 1.0, SignedDistanceAtTime<T>,
                                           params5, rInfo)
                      << "\n";*/
        }
    }
}

BOOST_AUTO_TEST_SUITE_END();
