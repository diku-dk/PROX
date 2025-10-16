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
        double longest = diag.maxCoeff();
        //10% padding to our bounding box!
        double pad = 0.10 * longest;
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

        BOOST_ASSERT(
            isClose(computeGradient_Working(EigenVector3<T>(0, 5, 0), G),
                    EigenVector3<T>(0, 1, 0)));
        BOOST_ASSERT(
            isClose(computeGradient_Working(EigenVector3<T>(0, -5, 0), G),
                    EigenVector3<T>(0, -1, 0)));

        BOOST_ASSERT(
            isClose(computeGradient_Working(EigenVector3<T>(5, 0, 0), G),
                    EigenVector3<T>(1, 0, 0)));
        BOOST_ASSERT(
            isClose(computeGradient_Working(EigenVector3<T>(-5, 0, 0), G),
                    EigenVector3<T>(-1, 0, 0)));

        BOOST_ASSERT(
            isClose(computeGradient_Working(EigenVector3<T>(0, 0, 5), G),
                    EigenVector3<T>(0, 0, 1)));
        BOOST_ASSERT(
            isClose(computeGradient_Working(EigenVector3<T>(0, 0, -5), G),
                    EigenVector3<T>(0, 0, -1)));

        std::cerr << computeGradient_Working(EigenVector3<T>(0, 5.0, 0), G)
                  << "\n";
        std::cerr << computeGradient_Working(EigenVector3<T>(0, 0.99, 0), G)
                  << "\n";
        std::cerr << computeGradient_Working(EigenVector3<T>(0, 0.5, 0), G)
                  << "\n";
        std::cerr << computeGradient_Working(EigenVector3<T>(0, 0.01, 0), G)
                  << "\n";
        std::cerr << computeGradient_Working(EigenVector3<T>(0, 0.00, 0), G)
                  << "\n";
        std::cerr << computeGradient_Working(EigenVector3<T>(0, -0.01, 0), G)
                  << "\n";
        std::cerr << computeGradient_Working(EigenVector3<T>(0, -0.5, 0), G)
                  << "\n";
        std::cerr << computeGradient_Working(EigenVector3<T>(0, -0.99, 0), G)
                  << "\n";
        std::cerr << computeGradient_Working(EigenVector3<T>(0, -5.0, 0), G)
                  << "\n";

        std::cerr << computeGradient_Working(EigenVector3<T>(3.0, 0, 0), G)
                  << "\n";
        // Expected: (1, 0, 0) - points away from +X face

        std::cerr << computeGradient_Working(EigenVector3<T>(0, 0, -4.0), G)
                  << "\n";
        // Expected: (0, 0, -1) - points away from -Z face

        // Test corners/edges
        std::cerr << computeGradient_Working(EigenVector3<T>(0.6, 0.6, 0), G)
                  << "\n";
        // Expected: normalized(1, 1, 0) - points away from edge between +X and +Y faces

        std::cerr << computeGradient_Working(EigenVector3<T>(0.6, 0.6, 0.6), G)
                  << "\n";
        // Expected: normalized(1, 1, 1) - points away from corner

        // Test points exactly at center
        std::cerr << computeGradient_Working(EigenVector3<T>(0, 10.0, 0), G)
                  << "\n";
        // Expected: (0, 1, 0) or similar - implementation dependent choice of direction

        // Test points very close to boundary from inside
        std::cerr << computeGradient_Working(EigenVector3<T>(0.49, 0, 0), G)
                  << "\n";
        // Expected: (1, 0, 0) - points toward +X boundary

        std::cerr << "SDF at (0.6,0.6,0.6): "
                  << value_at_2(G, EigenVector3<T>(0.6, 0.6, 0.6)) << "\n";
        std::cerr << "SDF at (0.6,0.6,0): "
                  << value_at_2(G, EigenVector3<T>(0.6, 0.6, 0)) << "\n";
        std::cerr << "SDF at (0.5,0,0): "
                  << value_at_2(G, EigenVector3<T>(0.5, 0, 0)) << "\n";
        std::cerr << "SDF at (0,0,0): "
                  << value_at_2(G, EigenVector3<T>(0, 0, 0)) << "\n";

        Eigen::Matrix<T, 3, 1> diff = (G.m_max - G.m_min);
        Eigen::Matrix<T, 3, 1> cell_size = Eigen::Matrix<T, 3, 1>(
            diff.x() / (G.m_nodes.x() - 1), diff.y() / (G.m_nodes.y() - 1),
            diff.z() / (G.m_nodes.z() - 1));
        T h = cell_size.x() * 0.1; // Very small step
        T dx = (value_at_2(G, EigenVector3<T>(0.6 + h, 0.6, 0.0))
                - value_at_2(G, EigenVector3<T>(0.6 - h, 0.6, 0.0)))
             / (2 * h);
        std::cerr << "dx: " << dx << "\n";
    }
}

BOOST_AUTO_TEST_SUITE_END();
