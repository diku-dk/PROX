#include "eigenhelperall.h"
#include "grid_enclosing_indices.h"
#include "grid_grid.h"
#include "grid_iterators.h"
#include "grid_node_position.h"
#include "grid_sdf_sdf_voxelize.hpp"
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
        fs::path bunnyRelative
            = "../../../../../bin/resources/objs/blender_star.obj";
        fs::path bunnyFull = source_dir / bunnyRelative;
        fs::path bunnyNormalized = bunnyFull.lexically_normal();

        // Output
        std::cout << "Bunny path: " << bunnyNormalized << std::endl;
        std::string meshFile = bunnyNormalized.string();

        int res = 64;

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

            Eigen::Matrix<T, 3, 1> diff = (G.m_max - G.m_min);
            Eigen::Matrix<T, 3, 1> cell_size = Eigen::Matrix<T, 3, 1>(
                diff.x() / (G.m_nodes.x() - 1), diff.y() / (G.m_nodes.y() - 1),
                diff.z() / (G.m_nodes.z() - 1));
            T cellSpacing = cell_size.x();

            uint32_t a, b, c;
            std::vector<SDFSDFContact::SDFVoxel<T>> voxels
                = SDFSDFContact::createVoxels(G.min(), G.max(), 8, cellSpacing,
                                              a, b, c);
            //Filter voxels;
            std::vector<SDFSDFContact::SDFVoxel<T>> filteredVoxels
                = SDFSDFContact::voxelFilterAll(voxels, G.min(), cellSpacing,
                                                G);

            SDFSDFContact::selectFeaturePointsPerVoxel(filteredVoxels, G.min(),
                                                       cellSpacing, a, b, c, G);

            int show = 70;
            int printed = 0;
            for (size_t i = 0; i < filteredVoxels.size() && printed < show; ++i)
            {
                const SDFSDFContact::SDFVoxel<T>& v = filteredVoxels[i];
                if (v.selected.empty()) continue;
                std::cout << "Voxel " << i << " bounds min(" << v.vmin.x()
                          << "," << v.vmin.y() << "," << v.vmin.z() << ") max("
                          << v.vmax.x() << "," << v.vmax.y() << ","
                          << v.vmax.z() << ")\n";
                for (size_t s = 0; s < v.selected.size(); ++s)
                {
                    const SDFSDFContact::SelectedPoint<T>& sp = v.selected[s];
                    std::string t;
                    switch (sp.type)
                    {
                    case SDFSDFContact::SelectedPoint<T>::PEAK:
                        t = "PEAK";
                        break;
                    case SDFSDFContact::SelectedPoint<T>::VALLEY:
                        t = "VALLEY";
                        break;
                    case SDFSDFContact::SelectedPoint<T>::DISCONTINUITY:
                        t = "DISCONTINUITY";
                        break;
                    case SDFSDFContact::SelectedPoint<T>::PROJECTED_ISOSURFACE:
                        t = "PROJECTED_ISOSURFACE";
                        break;
                    case SDFSDFContact::SelectedPoint<T>::CENTER_ISOSURFACE:
                        t = "CENTER_ISOSURFACE";
                        break;
                    case SDFSDFContact::SelectedPoint<T>::NEAREST_GRID:
                        t = "NEAREST_GRID";
                        break;
                    }
                    std::cout << "  [" << s << "] type=" << t << " pos=("
                              << sp.pos.x() << "," << sp.pos.y() << ","
                              << sp.pos.z() << ")" << " sdf=" << sp.sdf;
                    if (sp.gx >= 0)
                        std::cout << " grid=(" << sp.gx << "," << sp.gy << ","
                                  << sp.gz << ")";
                    std::cout << "\n";
                }
                printed++;
            }

            // Done
            std::cout << "Processed " << filteredVoxels.size()
                      << " voxels. Example printed: " << printed << "\n";
            BOOST_TEST(true);
        }

        {

            Eigen::Matrix<T, 3, 1> diff = (G.m_max - G.m_min);
            Eigen::Matrix<T, 3, 1> cell_size = Eigen::Matrix<T, 3, 1>(
                diff.x() / (G.m_nodes.x() - 1), diff.y() / (G.m_nodes.y() - 1),
                diff.z() / (G.m_nodes.z() - 1));
            T cellSpacing = cell_size.x();

            // voxels per axis
            const int Nv = 4;
            uint32_t nx, ny, nz;
            std::vector<SDFSDFContact::SDFVoxel<T>> voxels
                = SDFSDFContact::createVoxels(G.min(), G.max(), Nv, cellSpacing,
                                              nx, ny, nz);
            //Filter voxels;
            std::vector<SDFSDFContact::SDFVoxel<T>> filteredVoxels
                = SDFSDFContact::voxelFilterAll(voxels, G.min(), cellSpacing,
                                                G);

            SDFSDFContact::selectFeaturePointsPerVoxel(
                filteredVoxels, G.min(), cellSpacing, nx, ny, nz, G);

            //For writing!
            //First write SDF grid!
            size_t total = size_t(nx) * size_t(ny) * size_t(nz);
            std::vector<T> sdf_values;
            sdf_values.resize(total, 0.0);
            auto idx_flat = [&](int ix, int iy, int iz) -> size_t
            {
                // ix fastest, then iy, then iz
                return size_t(iz) * size_t(ny) * size_t(nx)
                     + size_t(iy) * size_t(nx) + size_t(ix);
            };

            EigenVector3<T> mmin = G.m_min;
            EigenVector3<T> mmax = G.m_max;
            for (int iz = 0; iz < int(nz); ++iz)
            {
                for (int iy = 0; iy < int(ny); ++iy)
                {
                    for (int ix = 0; ix < int(nx); ++ix)
                    {
                        EigenVector3<T> p = SDFSDFContact::gridIndexToPos<T>(
                            ix, iy, iz, G.m_min, cellSpacing);
                        T v = grid::value_at_2(G, p);
                        sdf_values[idx_flat(ix, iy, iz)] = v;
                    }
                }
            }

            // Write header
            {
                std::ofstream f("sdf_header.txt");
                f << std::setprecision(17);
                f << "mmin " << mmin.x() << " " << mmin.y() << " " << mmin.z()
                  << "\n";
                f << "mmax " << mmax.x() << " " << mmax.y() << " " << mmax.z()
                  << "\n";
                f << "spacing " << cellSpacing << "\n";
                f << "nx " << nx << " ny " << ny << " nz " << nz << "\n";
                f << "Nv " << Nv << "\n";
                f.close();
            }

            //Actual write SDF
            {
                std::ofstream f("sdf_values.bin", std::ios::binary);
                f.write(reinterpret_cast<const char*>(sdf_values.data()),
                        sdf_values.size() * sizeof(double));
                f.close();
            }

            //Write voxels.txt (one voxel per line)
            {
                std::ofstream f("voxels.txt");
                for (size_t i = 0; i < filteredVoxels.size(); ++i)
                {
                    const SDFSDFContact::SDFVoxel<T>& v = filteredVoxels[i];
                    f << v.vmin.x() << " " << v.vmin.y() << " " << v.vmin.z()
                      << " " << v.vmax.x() << " " << v.vmax.y() << " "
                      << v.vmax.z() << "\n";
                }
                f.close();
            }

            //Write points.txt (voxel_index pos.x pos.y pos.z sdf type_int)
            {
                std::ofstream f("points.txt");
                for (size_t i = 0; i < filteredVoxels.size(); ++i)
                {
                    SDFSDFContact::SDFVoxel<T>& v = filteredVoxels[i];
                    for (auto& sp : v.selected)
                    {
                        sp.pos = sp.pos
                               - grid::value_at_2(G, sp.pos)
                                     * grid::computeGradient_Working(sp.pos, G);
                        f << i << " " << sp.pos.x() << " " << sp.pos.y() << " "
                          << sp.pos.z() << " " << sp.sdf << " " << int(sp.type)
                          << "\n";
                    }
                }
                f.close();
            }

            // Write types mapping
            {
                std::ofstream f("types.txt");
                f << "0 PEAK\n1 VALLEY\n2 DISCONTINUITY\n3 "
                     "PROJECTED_ISOSURFACE\n4 CENTER_ISOSURFACE\n5 "
                     "NEAREST_GRID\n";
                f.close();
            }

            std::ofstream csv("sdf.csv");
            if (!csv)
            {
                std::cerr << "ERROR: cannot open " << outCsv
                          << " for writing\n";
            }
            csv << "x,y,z,sdf\n";
            idx_lin = 0;
            for (size_t k = 0; k < G.K(); ++k)
            {
                for (size_t j = 0; j < G.J(); ++j)
                {
                    for (size_t i = 0; i < G.I(); ++i)
                    {
                        Eigen::Matrix<size_t, 3, 1> idx(i, j, k);
                        Eigen::Matrix<T, 3, 1> p;
                        grid::node_position(G, idx, p);
                        csv << std::setprecision(9) << p.x() << "," << p.y()
                            << "," << p.z() << "," << G(idx) << "\n";
                        ++idx_lin;
                    }
                }
            }
            csv.close();
            std::cout << "Wrote: " << outCsv << " (" << total << " rows)\n";
        }
    }
}
BOOST_AUTO_TEST_SUITE_END();
