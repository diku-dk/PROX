#include "eigenhelperall.h"
#include "grid_enclosing_indices.h"
#include "grid_grid.h"
#include "grid_iterators.h"
#include "grid_node_position.h"
#include "grid_sdf_sdf_voxelize.hpp"
#include "grid_sdf_sdf_ccd_gradient_descent.hpp"
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

#include <grid_ccd_frank_wolfe_include_all.hpp>
#include <grid_local_optimization.hpp>
#include "igl/signed_distance.h"
#include "igl/marching_cubes.h"
#include "igl/writeOBJ.h"
//#include <libigl/signed_distance.h>
//#include <libigl/sphere.h>
//#include <libigl/marching_cubes.h>
//#include <libigl/writeOBJ.h>

#include <igl/read_triangle_mesh.h>
#include <grid_sdf_sdf_ccd_gradient.hpp>
#include <grid_local_optimization.hpp>

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

template <typename T>
void makeSDF(grid::Grid<T, T>& grid, std::string fileName,
             std::vector<SDFSDFContact::SDFVoxel<T>>& finishedVoxels,
             std::vector<SDFSDFContact::SDFVoxel<T>>& voxels, int res = 64,
             int voxelsCnt = 8)
{
    using D = T;
    namespace fs = std::filesystem;
    fs::path source_dir = fs::path(__FILE__).parent_path();
    fs::path bunnyRelative = "../../../../../bin/resources/objs/" + fileName;
    fs::path bunnyFull = source_dir / bunnyRelative;
    fs::path bunnyNormalized = bunnyFull.lexically_normal();

    // Output
    std::cout << "Bunny path: " << bunnyNormalized << std::endl;
    std::string meshFile = bunnyNormalized.string();

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
    Eigen::Matrix<size_t, 3, 1> nodes((size_t)res, (size_t)res, (size_t)res);
    grid::Grid<D, T> G;
    G.create(gmin, gmax, nodes);
    const size_t total = G.m_nodes.x() * G.m_nodes.y() * G.m_nodes.z();
    std::cout << "Created grid: " << G.I() << " x " << G.J() << " x " << G.K()
              << "  (total nodes = " << total << ")\n";

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

    std::cout << "Computing signed distances (libigl::signed_distance)...\n";
    Eigen::VectorXd S; // signed distances
    Eigen::VectorXi I; // indices to closest triangles (unused here)
    Eigen::MatrixXd C; // closest points on triangles
    Eigen::MatrixXd N; // normals used for signing (pseudonormal)

    igl::SignedDistanceType signType = igl::SIGNED_DISTANCE_TYPE_WINDING_NUMBER;
    igl::signed_distance(P, V, F, signType, S, I, C, N);
    if ((size_t)S.size() != total)
    {
        std::cerr << "ERROR: libigl returned unexpected S size: " << S.size()
                  << " expected " << total << "\n";
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

        Eigen::Matrix<T, 3, 1> diff = (G.m_max - G.m_min);
        Eigen::Matrix<T, 3, 1> cell_size = Eigen::Matrix<T, 3, 1>(
            diff.x() / (G.m_nodes.x() - 1), diff.y() / (G.m_nodes.y() - 1),
            diff.z() / (G.m_nodes.z() - 1));
        T cellSpacing = cell_size.x();

        // voxels per axis
        const int Nv = voxelsCnt;
        uint32_t nx, ny, nz;

        voxels = SDFSDFContact::createVoxels(G.min(), G.max(), Nv, cellSpacing,
                                             nx, ny, nz);
        //Filter voxels;
        finishedVoxels
            = SDFSDFContact::voxelFilterAll(voxels, G.min(), cellSpacing, G);

        SDFSDFContact::selectFeaturePointsPerVoxel(finishedVoxels, G.min(),
                                                   cellSpacing, nx, ny, nz, G);

        //FInish everything!
        grid = G;
    }
}

/*
BOOST_AUTO_TEST_CASE(grid_test_bunny_many_points_main)
{
    using T = double;
    //todo
    grid::Grid<T, T> SDFA;
    std::vector<SDFSDFContact::SDFVoxel<T>> finishedVoxelsA;
    std::vector<SDFSDFContact::SDFVoxel<T>> voxelsA;
    makeSDF(SDFA, "bunny.obj", finishedVoxelsA, voxelsA, 64, 8);
    std::cerr << "Finished filtering SDF 1/2!" << "\n";
    grid::Grid<T, T> SDFB;
    std::vector<SDFSDFContact::SDFVoxel<T>> voxelsB;
    std::vector<SDFSDFContact::SDFVoxel<T>> finishedVoxelsB;
    makeSDF(SDFB, "bunny.obj", finishedVoxelsB, voxelsB, 64, 8);
    std::cerr << "Finished filtering SDF 2/2!" << "\n";

    EigenVector3<T> angularVelA(0, 0, 0);
    EigenVector3<T> angularVelB(0, 0, 0);
    EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
    EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
    EigenVector3<T> transA(0.0, 0.0, 0.0);
    EigenVector3<T> transB(0.0, 0.5, 0.0);
    EigenVector3<T> linearA(0.0, 1.0, 0.0);
    EigenVector3<T> linearB(0, 0, 0);
    SDFSDFContact::SingleRigidBodyInfo<T> rInfoA;
    rInfoA.A_angularVel = &angularVelA;
    rInfoA.A_centerRotation = &rotA;
    rInfoA.A_centerTranslation = &transA;
    rInfoA.A_linearVel = &linearA;
    rInfoA.sdf = &SDFA;

    SDFSDFContact::SingleRigidBodyInfo<T> rInfoB;
    rInfoB.A_angularVel = &angularVelB;
    rInfoB.A_centerRotation = &rotB;
    rInfoB.A_centerTranslation = &transB;
    rInfoB.A_linearVel = &linearB;
    rInfoB.sdf = &SDFB;


    T toi = SDFSDFContact::getSDFSDFTOI(finishedVoxelsA, finishedVoxelsB,
                                        rInfoA, rInfoB, T(0.0), T(1.0));
    std::cerr << "GOT TOI = " << toi << "\n";
}*/

BOOST_AUTO_TEST_CASE(grid_test_sphere_correct_placement)
{
    //EXPECTED CONTACT POINT (0, 5, 0)!
    //With starting point (0, 1, 0) at 10 m/s, we expect contact to be at 0.4s
    using T = double;
    //todo
    grid::Grid<T, T> SDFA;
    std::vector<SDFSDFContact::SDFVoxel<T>> finishedVoxelsA;
    std::vector<SDFSDFContact::SDFVoxel<T>> voxelsA;
    makeSDF(SDFA, "blender_sphere.obj", finishedVoxelsA, voxelsA, 64, 8);
    std::cerr << "Finished filtering SDF 1/2!" << "\n";
    grid::Grid<T, T> SDFB;
    std::vector<SDFSDFContact::SDFVoxel<T>> voxelsB;
    std::vector<SDFSDFContact::SDFVoxel<T>> finishedVoxelsB;
    makeSDF(SDFB, "blender_sphere.obj", finishedVoxelsB, voxelsB, 64, 8);
    std::cerr << "Finished filtering SDF 2/2!" << "\n";

    EigenVector3<T> angularVelA(0, 0, 0);
    EigenVector3<T> angularVelB(0, 0, 0);
    EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
    EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
    EigenVector3<T> transA(0.0, 0.0, 0.0);
    EigenVector3<T> transB(0.0, 6.0, 0.0);
    EigenVector3<T> linearA(0.0, 10.0, 0.0);
    EigenVector3<T> linearB(0, 0, 0);
    SDFSDFContact::SingleRigidBodyInfo<T> rInfoA;
    rInfoA.A_angularVel = &angularVelA;
    rInfoA.A_centerRotation = &rotA;
    rInfoA.A_centerTranslation = &transA;
    rInfoA.A_linearVel = &linearA;
    rInfoA.sdf = &SDFA;

    SDFSDFContact::SingleRigidBodyInfo<T> rInfoB;
    rInfoB.A_angularVel = &angularVelB;
    rInfoB.A_centerRotation = &rotB;
    rInfoB.A_centerTranslation = &transB;
    rInfoB.A_linearVel = &linearB;
    rInfoB.sdf = &SDFB;

    for (size_t i = 0; i < finishedVoxelsA.size(); ++i)
    {
        for (size_t j = 0; j < finishedVoxelsA[i].selected.size(); ++j)
        {
            finishedVoxelsA[i].selected.clear();
        }
    }
    for (size_t i = 0; i < finishedVoxelsB.size(); ++i)
    {
        for (size_t j = 0; j < finishedVoxelsB[i].selected.size(); ++j)
        {
            finishedVoxelsB[i].selected.clear();
        }
    }
    SDFSDFContact::SelectedPoint<T> sp;
    sp.pos = EigenVector3<T>(0.0f, 1.0f, 0.0f);
    finishedVoxelsA[0].selected.push_back(sp);

    std::vector<EigenVector3<T>> dummy;
    T toi = SDFSDFContact::getSDFSDFTOI(finishedVoxelsA, finishedVoxelsB,
                                        rInfoA, rInfoB, T(0.0), T(1.0), dummy);
    std::cerr << "GOT TOI = " << toi << "\n";
}

BOOST_AUTO_TEST_CASE(grid_test_sphere_side_placement)
{
    using T = double;
    //todo
    //EXPECTED CONTACT POINT ???!
    //With starting point (0, 1, 0) at 10 m/s, we expect contact to be at 0.427s
    grid::Grid<T, T> SDFA;
    std::vector<SDFSDFContact::SDFVoxel<T>> finishedVoxelsA;
    std::vector<SDFSDFContact::SDFVoxel<T>> voxelsA;
    makeSDF(SDFA, "blender_sphere.obj", finishedVoxelsA, voxelsA, 64, 8);
    std::cerr << "Finished filtering SDF 1/2!" << "\n";
    grid::Grid<T, T> SDFB;
    std::vector<SDFSDFContact::SDFVoxel<T>> voxelsB;
    std::vector<SDFSDFContact::SDFVoxel<T>> finishedVoxelsB;
    makeSDF(SDFB, "blender_sphere.obj", finishedVoxelsB, voxelsB, 64, 8);
    std::cerr << "Finished filtering SDF 2/2!" << "\n";

    EigenVector3<T> angularVelA(0, 0, 0);
    EigenVector3<T> angularVelB(0, 0, 0);
    EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
    EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
    EigenVector3<T> transA(0.0, 0.0, 0.0);
    EigenVector3<T> transB(0.0, 6.0, 1.0);
    EigenVector3<T> linearA(0.0, 10.0, 0.0);
    EigenVector3<T> linearB(0, 0, 0);
    SDFSDFContact::SingleRigidBodyInfo<T> rInfoA;
    rInfoA.A_angularVel = &angularVelA;
    rInfoA.A_centerRotation = &rotA;
    rInfoA.A_centerTranslation = &transA;
    rInfoA.A_linearVel = &linearA;
    rInfoA.sdf = &SDFA;

    SDFSDFContact::SingleRigidBodyInfo<T> rInfoB;
    rInfoB.A_angularVel = &angularVelB;
    rInfoB.A_centerRotation = &rotB;
    rInfoB.A_centerTranslation = &transB;
    rInfoB.A_linearVel = &linearB;
    rInfoB.sdf = &SDFB;

    for (size_t i = 0; i < finishedVoxelsA.size(); ++i)
    {
        for (size_t j = 0; j < finishedVoxelsA[i].selected.size(); ++j)
        {
            finishedVoxelsA[i].selected.clear();
        }
    }
    for (size_t i = 0; i < finishedVoxelsB.size(); ++i)
    {
        for (size_t j = 0; j < finishedVoxelsB[i].selected.size(); ++j)
        {
            finishedVoxelsB[i].selected.clear();
        }
    }
    SDFSDFContact::SelectedPoint<T> sp;
    sp.pos = EigenVector3<T>(0.0f, 1.0f, 0.0);
    finishedVoxelsA[0].selected.push_back(sp);

    std::vector<EigenVector3<T>> dummy;
    T toi = SDFSDFContact::getSDFSDFTOI(finishedVoxelsA, finishedVoxelsB,
                                        rInfoA, rInfoB, T(0.0), T(1.0), dummy);
    std::cerr << "GOT TOI = " << toi << "\n";
}

BOOST_AUTO_TEST_CASE(grid_test_box_top_placement)
{
    //Configuration: Box colliding with box
    //With starting point (0, 0.5, 0) at 10 m/s, we expect contact to be at 0.6s.
    using T = double;
    grid::Grid<T, T> SDFA;
    std::vector<SDFSDFContact::SDFVoxel<T>> finishedVoxelsA;
    std::vector<SDFSDFContact::SDFVoxel<T>> voxelsA;
    makeSDF(SDFA, "box.obj", finishedVoxelsA, voxelsA, 64, 8);
    std::cerr << "Finished filtering SDF 1/2!" << "\n";
    grid::Grid<T, T> SDFB;
    std::vector<SDFSDFContact::SDFVoxel<T>> voxelsB;
    std::vector<SDFSDFContact::SDFVoxel<T>> finishedVoxelsB;
    makeSDF(SDFB, "box.obj", finishedVoxelsB, voxelsB, 64, 8);
    std::cerr << "Finished filtering SDF 2/2!" << "\n";

    EigenVector3<T> angularVelA(0, 0, 0);
    EigenVector3<T> angularVelB(0, 0, 0);
    EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
    EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
    EigenVector3<T> transA(0.0, 0.0, 0.0);
    EigenVector3<T> transB(0.0, 7.0, 0.0);
    EigenVector3<T> linearA(0.0, 10.0, 0.0);
    EigenVector3<T> linearB(0, 0, 0);
    SDFSDFContact::SingleRigidBodyInfo<T> rInfoA;
    rInfoA.A_angularVel = &angularVelA;
    rInfoA.A_centerRotation = &rotA;
    rInfoA.A_centerTranslation = &transA;
    rInfoA.A_linearVel = &linearA;
    rInfoA.sdf = &SDFA;

    SDFSDFContact::SingleRigidBodyInfo<T> rInfoB;
    rInfoB.A_angularVel = &angularVelB;
    rInfoB.A_centerRotation = &rotB;
    rInfoB.A_centerTranslation = &transB;
    rInfoB.A_linearVel = &linearB;
    rInfoB.sdf = &SDFB;

    for (size_t i = 0; i < finishedVoxelsA.size(); ++i)
    {
        for (size_t j = 0; j < finishedVoxelsA[i].selected.size(); ++j)
        {
            finishedVoxelsA[i].selected.clear();
        }
    }
    for (size_t i = 0; i < finishedVoxelsB.size(); ++i)
    {
        for (size_t j = 0; j < finishedVoxelsB[i].selected.size(); ++j)
        {
            finishedVoxelsB[i].selected.clear();
        }
    }
    SDFSDFContact::SelectedPoint<T> sp;
    sp.pos = EigenVector3<T>(0.0f, 0.5f, 0.0);
    finishedVoxelsA[0].selected.push_back(sp);

    std::vector<EigenVector3<T>> dummy;
    T toi = SDFSDFContact::getSDFSDFTOI(finishedVoxelsA, finishedVoxelsB,
                                        rInfoA, rInfoB, T(0.0), T(1.0), dummy);
    BOOST_TEST(std::abs<T>(toi - 0.6) <= 0.0001);
    std::cerr << "GOT TOI = " << toi << "\n";
}

BOOST_AUTO_TEST_CASE(grid_test_box_corner_placement)
{
    //Configuration: Box colliding with box
    //With starting point (0.5, 0.5, 0.5) at 10 m/s, we expect contact to be at 1.81 (some corner is cut off).
    using T = double;
    grid::Grid<T, T> SDFA;
    std::vector<SDFSDFContact::SDFVoxel<T>> finishedVoxelsA;
    std::vector<SDFSDFContact::SDFVoxel<T>> voxelsA;
    makeSDF(SDFA, "box.obj", finishedVoxelsA, voxelsA, 64, 8);
    std::cerr << "Finished filtering SDF 1/2!" << "\n";
    grid::Grid<T, T> SDFB;
    std::vector<SDFSDFContact::SDFVoxel<T>> voxelsB;
    std::vector<SDFSDFContact::SDFVoxel<T>> finishedVoxelsB;
    makeSDF(SDFB, "blender_sphere.obj", finishedVoxelsB, voxelsB, 64, 8);
    std::cerr << "Finished filtering SDF 2/2!" << "\n";

    EigenVector3<T> angularVelA(0, 0, 0);
    EigenVector3<T> angularVelB(0, 0, 0);
    EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
    EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
    EigenVector3<T> transA(0.0, 0.0, 0.0);
    EigenVector3<T> transB(1.0, 3.0, 1.0);
    EigenVector3<T> linearA(0.0, 10.0, 0.0);
    EigenVector3<T> linearB(0, 0, 0);
    SDFSDFContact::SingleRigidBodyInfo<T> rInfoA;
    rInfoA.A_angularVel = &angularVelA;
    rInfoA.A_centerRotation = &rotA;
    rInfoA.A_centerTranslation = &transA;
    rInfoA.A_linearVel = &linearA;
    rInfoA.sdf = &SDFA;

    SDFSDFContact::SingleRigidBodyInfo<T> rInfoB;
    rInfoB.A_angularVel = &angularVelB;
    rInfoB.A_centerRotation = &rotB;
    rInfoB.A_centerTranslation = &transB;
    rInfoB.A_linearVel = &linearB;
    rInfoB.sdf = &SDFB;

    for (size_t i = 0; i < finishedVoxelsA.size(); ++i)
    {
        for (size_t j = 0; j < finishedVoxelsA[i].selected.size(); ++j)
        {
            finishedVoxelsA[i].selected.clear();
        }
    }
    for (size_t i = 0; i < finishedVoxelsB.size(); ++i)
    {
        for (size_t j = 0; j < finishedVoxelsB[i].selected.size(); ++j)
        {
            finishedVoxelsB[i].selected.clear();
        }
    }
    SDFSDFContact::SelectedPoint<T> sp;
    sp.pos = EigenVector3<T>(0.5f, 0.5f, 0.5);
    finishedVoxelsA[0].selected.push_back(sp);

    std::vector<EigenVector3<T>> dummy;
    T toi = SDFSDFContact::getSDFSDFTOI(finishedVoxelsA, finishedVoxelsB,
                                        rInfoA, rInfoB, T(0.0), T(1.0), dummy);
    BOOST_TEST(std::abs<T>(toi - 0.181541) <= 0.0001);
    std::cerr << "GOT TOI = " << toi << "\n";
}

BOOST_AUTO_TEST_CASE(grid_test_box_non_corner_placement)
{
    //Configuration: Box colliding with box
    //With starting point (0.5, 0.5, 0.5) at 10 m/s, we expect contact to be at 1.81 (some corner is cut off).
    using T = double;
    grid::Grid<T, T> SDFA;
    std::vector<SDFSDFContact::SDFVoxel<T>> finishedVoxelsA;
    std::vector<SDFSDFContact::SDFVoxel<T>> voxelsA;
    makeSDF(SDFA, "box.obj", finishedVoxelsA, voxelsA, 64, 8);
    std::cerr << "Finished filtering SDF 1/2!" << "\n";
    grid::Grid<T, T> SDFB;
    std::vector<SDFSDFContact::SDFVoxel<T>> voxelsB;
    std::vector<SDFSDFContact::SDFVoxel<T>> finishedVoxelsB;
    makeSDF(SDFB, "blender_sphere.obj", finishedVoxelsB, voxelsB, 64, 8);
    std::cerr << "Finished filtering SDF 2/2!" << "\n";

    EigenVector3<T> angularVelA(0, 0, 0);
    EigenVector3<T> angularVelB(0, 0, 0);
    EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
    EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
    EigenVector3<T> transA(0.0, 0.0, 0.0);
    EigenVector3<T> transB(1.0, 3.0, 1.0);
    EigenVector3<T> linearA(0.0, 10.0, 0.0);
    EigenVector3<T> linearB(0, 0, 0);
    SDFSDFContact::SingleRigidBodyInfo<T> rInfoA;
    rInfoA.A_angularVel = &angularVelA;
    rInfoA.A_centerRotation = &rotA;
    rInfoA.A_centerTranslation = &transA;
    rInfoA.A_linearVel = &linearA;
    rInfoA.sdf = &SDFA;

    SDFSDFContact::SingleRigidBodyInfo<T> rInfoB;
    rInfoB.A_angularVel = &angularVelB;
    rInfoB.A_centerRotation = &rotB;
    rInfoB.A_centerTranslation = &transB;
    rInfoB.A_linearVel = &linearB;
    rInfoB.sdf = &SDFB;

    for (size_t i = 0; i < finishedVoxelsA.size(); ++i)
    {
        for (size_t j = 0; j < finishedVoxelsA[i].selected.size(); ++j)
        {
            finishedVoxelsA[i].selected.clear();
        }
    }
    for (size_t i = 0; i < finishedVoxelsB.size(); ++i)
    {
        for (size_t j = 0; j < finishedVoxelsB[i].selected.size(); ++j)
        {
            finishedVoxelsB[i].selected.clear();
        }
    }
    SDFSDFContact::SelectedPoint<T> sp;
    sp.pos = EigenVector3<T>(0.0f, 0.5f, 0.0);
    finishedVoxelsA[0].selected.push_back(sp);
    std::vector<EigenVector3<T>> dummy;
    T toi = SDFSDFContact::getSDFSDFTOI(finishedVoxelsA, finishedVoxelsB,
                                        rInfoA, rInfoB, T(0.0), T(1.0), dummy);
    BOOST_TEST(std::abs<T>(toi - 0.181541) <= 0.0001);
    std::cerr << "GOT TOI = " << toi << "\n";
}

BOOST_AUTO_TEST_CASE(grid_test_star_star)
{
    //Configuration: Box colliding with box
    //With starting point (0.5, 0.5, 0.5) at 10 m/s, we expect contact to be at.
    using T = double;
    grid::Grid<T, T> SDFA;
    std::vector<SDFSDFContact::SDFVoxel<T>> finishedVoxelsA;
    std::vector<SDFSDFContact::SDFVoxel<T>> voxelsA;
    makeSDF(SDFA, "blender_star.obj", finishedVoxelsA, voxelsA, 64, 8);
    std::cerr << "Finished filtering SDF 1/2!" << "\n";
    grid::Grid<T, T> SDFB;
    std::vector<SDFSDFContact::SDFVoxel<T>> voxelsB;
    std::vector<SDFSDFContact::SDFVoxel<T>> finishedVoxelsB;
    makeSDF(SDFB, "blender_star.obj", finishedVoxelsB, voxelsB, 64, 8);
    std::cerr << "Finished filtering SDF 2/2!" << "\n";

    EigenVector3<T> angularVelA(0, 0, 0);
    EigenVector3<T> angularVelB(0, 0, 0);
    EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
    EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
    EigenVector3<T> transA(0.0, 0.0, 0.0);
    EigenVector3<T> transB(0.0, 4.0, 0.0);
    EigenVector3<T> linearA(0.0, 10.0, 0.0);
    EigenVector3<T> linearB(0, 0, 0);
    SDFSDFContact::SingleRigidBodyInfo<T> rInfoA;
    rInfoA.A_angularVel = &angularVelA;
    rInfoA.A_centerRotation = &rotA;
    rInfoA.A_centerTranslation = &transA;
    rInfoA.A_linearVel = &linearA;
    rInfoA.sdf = &SDFA;

    SDFSDFContact::SingleRigidBodyInfo<T> rInfoB;
    rInfoB.A_angularVel = &angularVelB;
    rInfoB.A_centerRotation = &rotB;
    rInfoB.A_centerTranslation = &transB;
    rInfoB.A_linearVel = &linearB;
    rInfoB.sdf = &SDFB;

    for (size_t i = 0; i < finishedVoxelsA.size(); ++i)
    {
        for (size_t j = 0; j < finishedVoxelsA[i].selected.size(); ++j)
        {
            finishedVoxelsA[i].selected.clear();
        }
    }
    for (size_t i = 0; i < finishedVoxelsB.size(); ++i)
    {
        for (size_t j = 0; j < finishedVoxelsB[i].selected.size(); ++j)
        {
            finishedVoxelsB[i].selected.clear();
        }
    }
    SDFSDFContact::SelectedPoint<T> sp;
    sp.pos = EigenVector3<T>(0.001112, 0.797714, -0.57828);
    finishedVoxelsA[0].selected.push_back(sp);

    std::vector<EigenVector3<T>> dummy;
    T toi = SDFSDFContact::getSDFSDFTOI(finishedVoxelsA, finishedVoxelsB,
                                        rInfoA, rInfoB, T(0.0), T(1.0), dummy);
    BOOST_TEST(std::abs<T>(toi - 0.289893) <= 0.0001);
    std::cerr << "GOT TOI = " << toi << "\n";
}

BOOST_AUTO_TEST_SUITE_END();
