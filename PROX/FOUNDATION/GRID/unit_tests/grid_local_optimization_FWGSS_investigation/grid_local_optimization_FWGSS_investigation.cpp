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

#include <thread>
#include <vector>
#include <mutex>

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

template <typename T> struct ComputedTimers
{
    std::vector<T> minimizevals;
    T wholeFuncVal;
};

template <typename T> class CCDStatistics
{
public:
    int almostAtSolution = 0;
    int outsideButTOINotEnd = 0;
    int outsideAndTOIEqualsEnd = 0;
    int insideButTOINotStart = 0;
    int insideAndTOIEqualsStart = 0;
    int otherCases = 0;

    void collectStatistics(const std::vector<bool>& configurations)
    {
        // Reset counters
        almostAtSolution = 0;
        outsideButTOINotEnd = 0;
        outsideAndTOIEqualsEnd = 0;
        insideButTOINotStart = 0;
        insideAndTOIEqualsStart = 0;
        otherCases = 0;

        for (bool config : configurations)
        {
            // Simulate different scenarios based on random config
            if (config) { almostAtSolution++; }
            else { outsideButTOINotEnd++; }
        }
    }

    void printStatistics()
    {
        std::cout << "CCD Statistics:\n";
        std::cout << "1. Almost at solution (small interpenetration): "
                  << almostAtSolution << "\n";
        std::cout << "2. Outside solution, but TOI != tEnd: "
                  << outsideButTOINotEnd << "\n";
        std::cout << "3. TOI == tEnd and outside solution: "
                  << outsideAndTOIEqualsEnd << "\n";
        std::cout << "4. Inside solution, but TOI != tStart: "
                  << insideButTOINotStart << "\n";
        std::cout << "5. TOI == tStart and inside solution: "
                  << insideAndTOIEqualsStart << "\n";
        std::cout << "6. Other cases: " << otherCases << "\n";
    }

    void printAverages(const std::vector<ComputedTimers<T>>& timers)
    {
        if (timers.empty())
        {
            std::cerr << "No timers available.\n";
            return;
        }

        long double totalIterations = 0.0;
        long double totalWholeFuncVal = 0.0;
        long double totalMinimizeVal = 0.0;
        std::size_t totalMinimizeCount = 0;

        for (const auto& t : timers)
        {
            totalIterations += static_cast<long double>(t.minimizevals.size())
                             / static_cast<long double>(2.0);
            totalWholeFuncVal += static_cast<long double>(t.wholeFuncVal);

            for (const auto& v : t.minimizevals)
            {
                totalMinimizeVal += static_cast<long double>(v);
                ++totalMinimizeCount;
            }
        }

        long double avgIterations = totalIterations / timers.size();
        long double avgWholeFuncVal = totalWholeFuncVal / timers.size();
        long double avgMinimizeVal = totalMinimizeCount
                                       ? (totalMinimizeVal / totalMinimizeCount)
                                       : 0.0L;

        std::cerr << "Average iterations: " << avgIterations << '\n';
        std::cerr << "Average minimizevals value: " << avgMinimizeVal << '\n';
        std::cerr << "Average wholeFuncVal: " << avgWholeFuncVal << '\n';
    }
};

template <typename T> struct PointPenetration
{
    EigenVector3<T> point;
    T distToSDF;
};

template <typename T>
std::vector<PointPenetration<T>>
getPenetrationForTris(const EigenVector3<T>& A, const EigenVector3<T>& B,
                      const EigenVector3<T>& C, const RigidBodyInfo<T>& info,
                      int N = 446)
{
    if (N <= 0) throw std::invalid_argument("N must be positive");

    size_t M = static_cast<size_t>(N + 1) * static_cast<size_t>(N + 2) / 2;
    std::vector<PointPenetration<T>> pts;
    pts.reserve(M);

    for (int i = 0; i <= N; ++i)
    {
        int jmax = N - i;
        for (int j = 0; j <= jmax; ++j)
        {
            int k = N - i - j;
            T u = static_cast<T>(i) / static_cast<T>(N);
            T v = static_cast<T>(j) / static_cast<T>(N);
            T w = static_cast<T>(k) / static_cast<T>(N);
            PointPenetration<T> pointPenetration;
            pointPenetration.point = u * A + v * B + w * C;
            pointPenetration.distToSDF = valueAtProjection<T>(
                *(info.B_sdf), pointPenetration.point,
                *(info.B_centerTranslation), *(info.B_centerRotation));
            pts.push_back(pointPenetration);
        }
    }
    return pts;
}

template <typename T> class TriangleCCDTester
{
private:
    const T tStart = T(0.0);
    const T tEnd = T(0.01);
    const T sphereRadius = T(10.0);
    const T boxSize = T(2.0);

public:
    void runTests(int numTriangles, int type = 0)
    {
        CCDStatistics<T> stats;
        std::vector<bool> testResults;

        std::random_device rd;
        std::mt19937 gen(rd());

        using D = T;

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

        std::uniform_real_distribution<T> posDist(-sphereRadius, sphereRadius);
        std::uniform_real_distribution<T> smallAngleDist(-0.1, 0.1);
        std::uniform_real_distribution<T> largeAngleDist(-M_PI, M_PI);

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

        std::vector<ComputedTimers<T>> allComputedTimers;
        T minDiff = std::numeric_limits<T>::max();
        T maxDiff = std::numeric_limits<T>::min();
        std::vector<long double> allDifs(numTriangles);

        T minDiffDist = std::numeric_limits<T>::max();
        T maxDiffDist = std::numeric_limits<T>::min();
        std::vector<long double> allDifsDist(numTriangles);

        for (int i = 0; i < numTriangles; i++)
        {
            //Generate random triangle vertices around a sphere
            EigenVector3<T> center = randomPointOnSphere(gen, sphereRadius);

            //Create triangle oriented toward box center
            auto triangle = generateTriangleFacingOrigin(center, gen);

            //Compute velocity toward box center (origin)
            // Normalize so it reaches center at tEnd
            EigenVector3<T> velocity = (-center / tEnd) * T(1.0);

            transA = center;
            linearA = velocity;
            rInfo.A_centerTranslation = &transA;
            rInfo.A_linearVel = &velocity;
            rInfo.A_p0 = std::get<0>(triangle);
            rInfo.A_p1 = std::get<1>(triangle);
            rInfo.A_p2 = std::get<2>(triangle);

            T dt = 0.0;

            EigenVector3<T> lastContactPoint;
            while (dt <= 0.01000)
            {
                TriangleAtTimeInfo<T> tri = getTriangleAtTime(dt, rInfo);
                EigenVector3<T> p0AtTimeTi = tri.A_p0;
                EigenVector3<T> p1AtTimeTi = tri.A_p1;
                EigenVector3<T> p2AtTimeTi = tri.A_p2;

                EigenVector3<T> contactPoint;
                EigenVector3<T> normalDummy;
                T penetration;
                /*optimizeTriangleFW_Working(p0AtTimeTi, p1AtTimeTi, p2AtTimeTi,
                                           *(rInfo.B_sdf), contactPoint,
                                           normalDummy, penetration, 1000);
                T penetrations = valueAtProjection(*(rInfo.B_sdf), contactPoint,
                                                   *(rInfo.B_centerTranslation),
                                                   *(rInfo.B_centerRotation));*/
                std::vector<PointPenetration<T>> pens = getPenetrationForTris(
                    p0AtTimeTi, p1AtTimeTi, p2AtTimeTi, rInfo, 100);
                T minPenetration = std::numeric_limits<T>::max();
                EigenVector3<T> minPoint;
                for (size_t i = 0; i < pens.size(); ++i)
                {
                    PointPenetration<T> currPen = pens[i];
                    if (currPen.distToSDF < minPenetration)
                    {
                        minPoint = currPen.point;
                        minPenetration = currPen.distToSDF;
                    }
                }
                lastContactPoint = minPoint;
                if (minPenetration <= 0.0) { break; }
                dt += 1e-1;
            }
            std::cerr << ".";
            // Run CCD
            EigenVector3<T> firstIntersectPoint;
            std::vector<T> minimizerTimes;
            T totalAlgorithmTime;
            T toi;
            if (type == 0)
            {
                auto gss_start = std::chrono::high_resolution_clock::now();

                toi = FrankWolfeGSS_BENCHMARK_TIME<T>(
                    0.0, 0.01, rInfo, firstIntersectPoint, minimizerTimes);
                auto gss_end = std::chrono::high_resolution_clock::now();
                auto gss_us
                    = std::chrono::duration_cast<std::chrono::nanoseconds>(
                          gss_end - gss_start)
                          .count();
                totalAlgorithmTime = T(gss_us);
            }

            else if (type == 1)
            {
                auto gss_start = std::chrono::high_resolution_clock::now();

                toi = FrankWolfeBacktracking_BENCHMARK_TIME<T>(
                    0.0, 0.01, rInfo, firstIntersectPoint, minimizerTimes);
                auto gss_end = std::chrono::high_resolution_clock::now();
                auto gss_us
                    = std::chrono::duration_cast<std::chrono::nanoseconds>(
                          gss_end - gss_start)
                          .count();
                totalAlgorithmTime = T(gss_us);
            }
            else if (type == 2)
            {
                auto gss_start = std::chrono::high_resolution_clock::now();

                toi = FrankWolfeGSSBisection_BENCHMARK_TIME<T>(
                    0.0, 0.01, rInfo, firstIntersectPoint, minimizerTimes);
                auto gss_end = std::chrono::high_resolution_clock::now();
                auto gss_us
                    = std::chrono::duration_cast<std::chrono::nanoseconds>(
                          gss_end - gss_start)
                          .count();
                totalAlgorithmTime = T(gss_us);
            }
            else if (type == 3)
            {
                auto gss_start = std::chrono::high_resolution_clock::now();

                toi = performProjectedGradientDescent<T>(
                    0.0, 0.01, rInfo, firstIntersectPoint, minimizerTimes);
                auto gss_end = std::chrono::high_resolution_clock::now();
                auto gss_us
                    = std::chrono::duration_cast<std::chrono::nanoseconds>(
                          gss_end - gss_start)
                          .count();
                totalAlgorithmTime = T(gss_us);
            }
            else
            {
                auto gss_start = std::chrono::high_resolution_clock::now();

                toi = FrankWolfeBRENT_BENCHMARK_TIME<T>(
                    0.0, 0.01, rInfo, firstIntersectPoint, minimizerTimes);
                toi = std::min<T>(toi, 0.01);
                auto gss_end = std::chrono::high_resolution_clock::now();
                auto gss_us
                    = std::chrono::duration_cast<std::chrono::nanoseconds>(
                          gss_end - gss_start)
                          .count();
                totalAlgorithmTime = T(gss_us);
            }
            ComputedTimers<T> computedTimers;
            computedTimers.minimizevals = minimizerTimes;
            computedTimers.wholeFuncVal = totalAlgorithmTime;
            allComputedTimers.push_back(computedTimers);

            T dist = grid::valueAtProjection(G, firstIntersectPoint,
                                             *(rInfo.B_centerTranslation),
                                             *(rInfo.B_centerRotation));
            T distGT = grid::valueAtProjection(G, lastContactPoint,
                                               *(rInfo.B_centerTranslation),
                                               *(rInfo.B_centerRotation));

            T diffDistance = dist - distGT;

            if (dist < 0.0 && dist < -T(1e-3))
            {
                std::cerr << "XTI: " << firstIntersectPoint << "\n";
                std::cerr << "TOI: " << toi << "\n";
                std::cerr << "DIST: " << dist << "\n";
            }
            bool hasIntersection = (toi >= 0.01);
            T diff = toi - dt;

            if (diffDistance <= -1e-5)
            {
                std::cerr << "=============================================\n";
                std::cerr
                    << "Triangle configuration of diffDistance <= -1e-5: p0=("
                    << rInfo.A_p0.x() << ", " << rInfo.A_p0.y() << ", "
                    << rInfo.A_p0.z() << "), p1=(" << rInfo.A_p1.x() << ", "
                    << rInfo.A_p1.y() << ", " << rInfo.A_p1.z() << "), p2=("
                    << rInfo.A_p2.x() << ", " << rInfo.A_p2.y() << ", "
                    << rInfo.A_p2.z() << "), velocity=("
                    << (*(rInfo.A_linearVel)).x() << ", "
                    << (*(rInfo.A_linearVel)).y() << ", "
                    << (*(rInfo.A_linearVel)).z()
                    << "), exact TOI vs dt = " << toi << " vs " << dt << "\n";
                std::cerr << "Note diffDistance is" << diffDistance
                          << ", diffGT is " << distGT << " and diff is " << dist
                          << "\n";
                std::cerr << "=============================================\n";
            }
            if (diffDistance >= 1e-6)
            {
                std::cerr << "=============================================\n";
                std::cerr
                    << "Triangle configuration of diffDistance >= 1e-6: p0=("
                    << rInfo.A_p0.x() << ", " << rInfo.A_p0.y() << ", "
                    << rInfo.A_p0.z() << "), p1=(" << rInfo.A_p1.x() << ", "
                    << rInfo.A_p1.y() << ", " << rInfo.A_p1.z() << "), p2=("
                    << rInfo.A_p2.x() << ", " << rInfo.A_p2.y() << ", "
                    << rInfo.A_p2.z() << "), velocity=("
                    << (*(rInfo.A_linearVel)).x() << ", "
                    << (*(rInfo.A_linearVel)).y() << ", "
                    << (*(rInfo.A_linearVel)).z()
                    << "), exact TOI vs dt = " << toi << " vs " << dt << "\n";
                std::cerr << "Note diffDistance is" << diffDistance
                          << ", diffGT is " << distGT << " and diff is " << dist
                          << "\n";
                std::cerr << "=============================================\n";
            }

            if (toi > dt + 1e-5)
            {
                std::cerr << "=============================================\n";
                std::cerr
                    << "WE GOT A TOI that is larger thann DT! (1e5 prec)\n";
                std::cerr << "Triangle configuration of TOI=1e-5: p0=("
                          << rInfo.A_p0.x() << ", " << rInfo.A_p0.y() << ", "
                          << rInfo.A_p0.z() << "), p1=(" << rInfo.A_p1.x()
                          << ", " << rInfo.A_p1.y() << ", " << rInfo.A_p1.z()
                          << "), p2=(" << rInfo.A_p2.x() << ", "
                          << rInfo.A_p2.y() << ", " << rInfo.A_p2.z()
                          << "), velocity=(" << (*(rInfo.A_linearVel)).x()
                          << ", " << (*(rInfo.A_linearVel)).y() << ", "
                          << (*(rInfo.A_linearVel)).z()
                          << "), exact TOI vs dt = " << toi << " vs " << dt
                          << "\n";
                std::cerr << "=============================================\n";
            }

            else if (toi > dt + 1e-8)
            {
                std::cerr
                    << "WE GOT A TOI that is larger thann DT! (1e8 prec)\n";
            }
            else if (toi > dt + 1e-11)
            {
                std::cerr
                    << "WE GOT A TOI that is larger thann DT! (no prec)\n";
            }

            minDiff = std::min<T>(minDiff, diff);
            maxDiff = std::max<T>(maxDiff, diff);
            allDifs.push_back(static_cast<long double>(diff));

            minDiffDist = std::min<T>(minDiffDist, diffDistance);
            maxDiffDist = std::max<T>(maxDiffDist, diffDistance);
            allDifsDist.push_back(static_cast<long double>(diffDistance));
            classifyResult(stats, dist, toi, hasIntersection);
        }

        stats.printStatistics();
        stats.printAverages(allComputedTimers);

        {
            // --- Mean ---
            long double mean
                = std::accumulate(allDifs.begin(), allDifs.end(), 0.0L)
                / allDifs.size();

            // --- Median ---
            std::vector<long double> sorted = allDifs;
            std::sort(sorted.begin(), sorted.end());
            long double median;
            long double two = 2.0;
            long double one = 1.0;
            if (sorted.size() % 2 == 0)
                median = (sorted[sorted.size() / two - one]
                          + sorted[sorted.size() / two])
                       / 2.0L;
            else
                median = sorted[sorted.size() / two];

            // --- Standard Deviation ---
            long double sumSqDiff = 0.0L;
            for (auto x : allDifs) sumSqDiff += (x - mean) * (x - mean);
            long double stddev = std::sqrtl(sumSqDiff / (allDifs.size() - 1));

            // --- Output ---
            std::cout << "Mean diff of toi:   " << mean << "\n";
            std::cout << "Median diff of toi: " << std::setprecision(20)
                      << median << "\n";
            std::cout << "StdDev diff of toi: " << stddev << "\n";
            std::cout << "min diff of toi: " << minDiff << "\n";
            std::cout << "max diff of toi: " << maxDiff << "\n";
        }
        {
            // --- Mean ---
            long double mean
                = std::accumulate(allDifsDist.begin(), allDifsDist.end(), 0.0L)
                / allDifsDist.size();

            // --- Median ---
            std::vector<long double> sorted = allDifsDist;
            std::sort(sorted.begin(), sorted.end());
            long double median;
            long double two = 2.0;
            long double one = 1.0;
            if (sorted.size() % 2 == 0)
                median = (sorted[sorted.size() / two - one]
                          + sorted[sorted.size() / two])
                       / two;
            else
                median = sorted[sorted.size() / two];

            // --- Standard Deviation ---
            long double sumSqDiff = 0.0L;
            for (auto x : allDifsDist) sumSqDiff += (x - mean) * (x - mean);
            long double stddev
                = std::sqrtl(sumSqDiff / (allDifsDist.size() - 1));

            // --- Output ---
            std::cout << "Mean diff (phi (x)) of toi:   " << mean << "\n";
            std::cout << "Median diff (phi (x)) of toi: "
                      << std::setprecision(20) << median << "\n";
            std::cout << "StdDev diff (phi (x)) of toi: " << stddev << "\n";
            std::cout << "min diff (phi (x)) of toi: " << minDiffDist << "\n";
            std::cout << "max diff (phi (x)) of toi: " << maxDiffDist << "\n";
        }
    }

    void runTestsParallel(int numTriangles, int type = 0)
    {
        CCDStatistics<T> stats;
        std::vector<bool> testResults;

        std::random_device rd;
        std::mt19937 gen(rd());

        using D = T;

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

        std::uniform_real_distribution<T> posDist(-sphereRadius, sphereRadius);
        std::uniform_real_distribution<T> smallAngleDist(-0.1, 0.1);
        std::uniform_real_distribution<T> largeAngleDist(-M_PI, M_PI);

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

        std::vector<ComputedTimers<T>> allComputedTimers;
        T minDiff = std::numeric_limits<T>::max();
        T maxDiff = std::numeric_limits<T>::min();
        std::vector<long double> allDifs(numTriangles);

        T minDiffDist = std::numeric_limits<T>::max();
        T maxDiffDist = std::numeric_limits<T>::min();
        std::vector<long double> allDifsDist(numTriangles);

        {
            std::vector<std::thread> threads;
            const int numThreads = 14;
            std::mutex mutex;

            // Variables that need to be protected by mutex
            std::vector<ComputedTimers<T>> thread_allComputedTimers;
            T thread_minDiff = std::numeric_limits<T>::max();
            T thread_maxDiff = std::numeric_limits<T>::min();
            std::vector<long double> thread_allDifs;
            T thread_minDiffDist = std::numeric_limits<T>::max();
            T thread_maxDiffDist = std::numeric_limits<T>::min();
            std::vector<long double> thread_allDifsDist;
            CCDStatistics<T> thread_stats;

            auto worker = [&](int start, int end)
            {
                // Local variables for each thread
                std::vector<ComputedTimers<T>> local_allComputedTimers;
                T local_minDiff = std::numeric_limits<T>::max();
                T local_maxDiff = std::numeric_limits<T>::min();
                std::vector<long double> local_allDifs;
                T local_minDiffDist = std::numeric_limits<T>::max();
                T local_maxDiffDist = std::numeric_limits<T>::min();
                std::vector<long double> local_allDifsDist;
                CCDStatistics<T> local_stats;

                // Thread-local random number generator
                std::random_device rd;
                std::mt19937 gen(rd());

                for (int i = start; i < end; i++)
                {
                    //Generate random triangle vertices around a sphere
                    EigenVector3<T> center
                        = randomPointOnSphere(gen, sphereRadius);

                    //Create triangle oriented toward box center
                    auto triangle = generateTriangleFacingOrigin(center, gen);

                    //Compute velocity toward box center (origin)
                    // Normalize so it reaches center at tEnd
                    EigenVector3<T> velocity = (-center / tEnd) * T(1.0);

                    // Create thread-local copies of these variables since they're modified
                    EigenVector3<T> local_transA = center;
                    EigenVector3<T> local_linearA = velocity;

                    // Create a local copy of rInfo for thread safety
                    RigidBodyInfo<T> local_rInfo = rInfo;
                    local_rInfo.A_centerTranslation = &local_transA;
                    local_rInfo.A_linearVel = &local_linearA;
                    local_rInfo.A_p0 = std::get<0>(triangle);
                    local_rInfo.A_p1 = std::get<1>(triangle);
                    local_rInfo.A_p2 = std::get<2>(triangle);

                    T dt = 0.0;

                    EigenVector3<T> lastContactPoint;
                    while (dt <= 0.01000)
                    {
                        TriangleAtTimeInfo<T> tri
                            = getTriangleAtTime(dt, local_rInfo);
                        EigenVector3<T> p0AtTimeTi = tri.A_p0;
                        EigenVector3<T> p1AtTimeTi = tri.A_p1;
                        EigenVector3<T> p2AtTimeTi = tri.A_p2;

                        EigenVector3<T> contactPoint;
                        EigenVector3<T> normalDummy;
                        T penetration;
                        std::vector<PointPenetration<T>> pens
                            = getPenetrationForTris(p0AtTimeTi, p1AtTimeTi,
                                                    p2AtTimeTi, local_rInfo,
                                                    150);
                        T minPenetration = std::numeric_limits<T>::max();
                        EigenVector3<T> minPoint;
                        for (size_t j = 0; j < pens.size(); ++j)
                        {
                            PointPenetration<T> currPen = pens[j];
                            if (currPen.distToSDF < minPenetration)
                            {
                                minPoint = currPen.point;
                                minPenetration = currPen.distToSDF;
                            }
                        }
                        lastContactPoint = minPoint;
                        if (minPenetration <= 0.0) { break; }
                        dt += 1e-6 * 0.5;
                    }

                    // Protect std::cerr output with mutex
                    {
                        std::lock_guard<std::mutex> lock(mutex);
                        std::cerr << ".";
                    }

                    // Run CCD
                    EigenVector3<T> firstIntersectPoint;
                    std::vector<T> minimizerTimes;
                    T totalAlgorithmTime;
                    T toi;
                    if (type == 0)
                    {
                        auto gss_start
                            = std::chrono::high_resolution_clock::now();
                        toi = FrankWolfeGSS_BENCHMARK_TIME<T>(
                            0.0, 0.01, local_rInfo, firstIntersectPoint,
                            minimizerTimes);
                        auto gss_end
                            = std::chrono::high_resolution_clock::now();
                        auto gss_us
                            = std::chrono::duration_cast<
                                  std::chrono::nanoseconds>(gss_end - gss_start)
                                  .count();
                        totalAlgorithmTime = T(gss_us);
                    }
                    else if (type == 1)
                    {
                        auto gss_start
                            = std::chrono::high_resolution_clock::now();
                        toi = FrankWolfeBacktracking_BENCHMARK_TIME<T>(
                            0.0, 0.01, local_rInfo, firstIntersectPoint,
                            minimizerTimes);
                        auto gss_end
                            = std::chrono::high_resolution_clock::now();
                        auto gss_us
                            = std::chrono::duration_cast<
                                  std::chrono::nanoseconds>(gss_end - gss_start)
                                  .count();
                        totalAlgorithmTime = T(gss_us);
                    }
                    else if (type == 2)
                    {
                        auto gss_start
                            = std::chrono::high_resolution_clock::now();
                        toi = FrankWolfeGSSBisection_BENCHMARK_TIME<T>(
                            0.0, 0.01, local_rInfo, firstIntersectPoint,
                            minimizerTimes);
                        auto gss_end
                            = std::chrono::high_resolution_clock::now();
                        auto gss_us
                            = std::chrono::duration_cast<
                                  std::chrono::nanoseconds>(gss_end - gss_start)
                                  .count();
                        totalAlgorithmTime = T(gss_us);
                    }
                    else if (type == 3)
                    {
                        auto gss_start
                            = std::chrono::high_resolution_clock::now();
                        toi = performProjectedGradientDescent<T>(
                            0.0, 0.01, local_rInfo, firstIntersectPoint,
                            minimizerTimes);
                        auto gss_end
                            = std::chrono::high_resolution_clock::now();
                        auto gss_us
                            = std::chrono::duration_cast<
                                  std::chrono::nanoseconds>(gss_end - gss_start)
                                  .count();
                        totalAlgorithmTime = T(gss_us);
                    }
                    else
                    {
                        auto gss_start
                            = std::chrono::high_resolution_clock::now();
                        toi = FrankWolfeBRENT_BENCHMARK_TIME<T>(
                            0.0, 0.01, local_rInfo, firstIntersectPoint,
                            minimizerTimes);
                        toi = std::min<T>(toi, 0.01);
                        auto gss_end
                            = std::chrono::high_resolution_clock::now();
                        auto gss_us
                            = std::chrono::duration_cast<
                                  std::chrono::nanoseconds>(gss_end - gss_start)
                                  .count();
                        totalAlgorithmTime = T(gss_us);
                    }

                    ComputedTimers<T> computedTimers;
                    computedTimers.minimizevals = minimizerTimes;
                    computedTimers.wholeFuncVal = totalAlgorithmTime;
                    local_allComputedTimers.push_back(computedTimers);

                    T dist = grid::valueAtProjection(
                        G, firstIntersectPoint,
                        *(local_rInfo.B_centerTranslation),
                        *(local_rInfo.B_centerRotation));
                    T distGT = grid::valueAtProjection(
                        G, lastContactPoint, *(local_rInfo.B_centerTranslation),
                        *(local_rInfo.B_centerRotation));

                    T diffDistance = dist - distGT;

                    if (dist < 0.0 && dist < -T(1e-3))
                    {
                        std::lock_guard<std::mutex> lock(mutex);
                        std::cerr << "XTI: " << firstIntersectPoint << "\n";
                        std::cerr << "TOI: " << toi << "\n";
                        std::cerr << "DIST: " << dist << "\n";
                    }
                    bool hasIntersection = (toi >= 0.01);
                    T diff = toi - dt;

                    if (diffDistance <= -1e-5)
                    {
                        std::lock_guard<std::mutex> lock(mutex);
                        std::cerr << "========================================="
                                     "====\n";
                        std::cerr << "Triangle configuration of diffDistance "
                                     "<= -1e-5: p0=("
                                  << local_rInfo.A_p0.x() << ", "
                                  << local_rInfo.A_p0.y() << ", "
                                  << local_rInfo.A_p0.z() << "), p1=("
                                  << local_rInfo.A_p1.x() << ", "
                                  << local_rInfo.A_p1.y() << ", "
                                  << local_rInfo.A_p1.z() << "), p2=("
                                  << local_rInfo.A_p2.x() << ", "
                                  << local_rInfo.A_p2.y() << ", "
                                  << local_rInfo.A_p2.z() << "), velocity=("
                                  << (*(local_rInfo.A_linearVel)).x() << ", "
                                  << (*(local_rInfo.A_linearVel)).y() << ", "
                                  << (*(local_rInfo.A_linearVel)).z()
                                  << "), exact TOI vs dt = " << toi << " vs "
                                  << dt << "\n";
                        std::cerr << "Note diffDistance is" << diffDistance
                                  << ", diffGT is " << distGT << " and diff is "
                                  << dist << "\n";
                        std::cerr << "========================================="
                                     "====\n";
                    }
                    if (diffDistance >= 1e-2)
                    {
                        std::lock_guard<std::mutex> lock(mutex);
                        std::cerr << "========================================="
                                     "====\n";
                        std::cerr << "Triangle configuration of diffDistance "
                                     ">= 1e-2: p0=("
                                  << local_rInfo.A_p0.x() << ", "
                                  << local_rInfo.A_p0.y() << ", "
                                  << local_rInfo.A_p0.z() << "), p1=("
                                  << local_rInfo.A_p1.x() << ", "
                                  << local_rInfo.A_p1.y() << ", "
                                  << local_rInfo.A_p1.z() << "), p2=("
                                  << local_rInfo.A_p2.x() << ", "
                                  << local_rInfo.A_p2.y() << ", "
                                  << local_rInfo.A_p2.z() << "), velocity=("
                                  << (*(local_rInfo.A_linearVel)).x() << ", "
                                  << (*(local_rInfo.A_linearVel)).y() << ", "
                                  << (*(local_rInfo.A_linearVel)).z()
                                  << "), exact TOI vs dt = " << toi << " vs "
                                  << dt << "\n";
                        std::cerr << "Note diffDistance is" << diffDistance
                                  << ", diffGT is " << distGT << " and diff is "
                                  << dist << "\n";
                        std::cerr << "========================================="
                                     "====\n";
                    }

                    if (toi > dt + 1e-5)
                    {
                        std::lock_guard<std::mutex> lock(mutex);
                        std::cerr << "========================================="
                                     "====\n";
                        std::cerr << "WE GOT A TOI that is larger thann DT! "
                                     "(1e5 prec)\n";
                        std::cerr << "Triangle configuration of TOI=1e-5: p0=("
                                  << local_rInfo.A_p0.x() << ", "
                                  << local_rInfo.A_p0.y() << ", "
                                  << local_rInfo.A_p0.z() << "), p1=("
                                  << local_rInfo.A_p1.x() << ", "
                                  << local_rInfo.A_p1.y() << ", "
                                  << local_rInfo.A_p1.z() << "), p2=("
                                  << local_rInfo.A_p2.x() << ", "
                                  << local_rInfo.A_p2.y() << ", "
                                  << local_rInfo.A_p2.z() << "), velocity=("
                                  << (*(local_rInfo.A_linearVel)).x() << ", "
                                  << (*(local_rInfo.A_linearVel)).y() << ", "
                                  << (*(local_rInfo.A_linearVel)).z()
                                  << "), exact TOI vs dt = " << toi << " vs "
                                  << dt << "\n";
                        std::cerr << "========================================="
                                     "====\n";
                    }
                    else if (toi > dt + 1e-8)
                    {
                        std::lock_guard<std::mutex> lock(mutex);
                        std::cerr << "WE GOT A TOI that is larger thann DT! "
                                     "(1e8 prec)\n";
                    }
                    else if (toi > dt + 1e-11)
                    {
                        std::lock_guard<std::mutex> lock(mutex);
                        std::cerr << "WE GOT A TOI that is larger thann DT! "
                                     "(no prec)\n";
                    }

                    local_minDiff = std::min<T>(local_minDiff, diff);
                    local_maxDiff = std::max<T>(local_maxDiff, diff);
                    local_allDifs.push_back(static_cast<long double>(diff));

                    local_minDiffDist
                        = std::min<T>(local_minDiffDist, diffDistance);
                    local_maxDiffDist
                        = std::max<T>(local_maxDiffDist, diffDistance);
                    local_allDifsDist.push_back(
                        static_cast<long double>(diffDistance));
                    classifyResult(local_stats, dist, toi, hasIntersection);
                }

                // Merge local results into shared variables
                std::lock_guard<std::mutex> lock(mutex);
                thread_allComputedTimers.insert(thread_allComputedTimers.end(),
                                                local_allComputedTimers.begin(),
                                                local_allComputedTimers.end());
                thread_minDiff = std::min(thread_minDiff, local_minDiff);
                thread_maxDiff = std::max(thread_maxDiff, local_maxDiff);
                thread_allDifs.insert(thread_allDifs.end(),
                                      local_allDifs.begin(),
                                      local_allDifs.end());
                thread_minDiffDist
                    = std::min(thread_minDiffDist, local_minDiffDist);
                thread_maxDiffDist
                    = std::max(thread_maxDiffDist, local_maxDiffDist);
                thread_allDifsDist.insert(thread_allDifsDist.end(),
                                          local_allDifsDist.begin(),
                                          local_allDifsDist.end());

                // Merge statistics
                thread_stats.almostAtSolution += local_stats.almostAtSolution;
                thread_stats.outsideButTOINotEnd
                    += local_stats.outsideButTOINotEnd;
                thread_stats.outsideAndTOIEqualsEnd
                    += local_stats.outsideAndTOIEqualsEnd;
                thread_stats.insideButTOINotStart
                    += local_stats.insideButTOINotStart;
                thread_stats.insideAndTOIEqualsStart
                    += local_stats.insideAndTOIEqualsStart;
                thread_stats.otherCases += local_stats.otherCases;
            };

            // Calculate work distribution
            int chunkSize = numTriangles / numThreads;
            int remainder = numTriangles % numThreads;
            int start = 0;

            for (int t = 0; t < numThreads; t++)
            {
                int end = start + chunkSize + (t < remainder ? 1 : 0);
                threads.emplace_back(worker, start, end);
                start = end;
            }

            // Wait for all threads to complete
            for (auto& thread : threads) { thread.join(); }

            // Copy thread results back to original variables
            allComputedTimers = std::move(thread_allComputedTimers);
            minDiff = thread_minDiff;
            maxDiff = thread_maxDiff;
            allDifs = std::move(thread_allDifs);
            minDiffDist = thread_minDiffDist;
            maxDiffDist = thread_maxDiffDist;
            allDifsDist = std::move(thread_allDifsDist);
            stats = thread_stats;
        }

        stats.printStatistics();
        stats.printAverages(allComputedTimers);

        {
            // --- Mean ---
            long double mean
                = std::accumulate(allDifs.begin(), allDifs.end(), 0.0L)
                / allDifs.size();

            // --- Median ---
            std::vector<long double> sorted = allDifs;
            std::sort(sorted.begin(), sorted.end());
            long double median;
            long double two = 2.0;
            long double one = 1.0;
            if (sorted.size() % 2 == 0)
                median = (sorted[sorted.size() / two - one]
                          + sorted[sorted.size() / two])
                       / 2.0L;
            else
                median = sorted[sorted.size() / two];

            // --- Standard Deviation ---
            long double sumSqDiff = 0.0L;
            for (auto x : allDifs) sumSqDiff += (x - mean) * (x - mean);
            long double stddev = std::sqrtl(sumSqDiff / (allDifs.size() - 1));

            // --- Output ---
            std::cout << "Mean diff of toi:   " << mean << "\n";
            std::cout << "Median diff of toi: " << std::setprecision(20)
                      << median << "\n";
            std::cout << "StdDev diff of toi: " << stddev << "\n";
            std::cout << "min diff of toi: " << minDiff << "\n";
            std::cout << "max diff of toi: " << maxDiff << "\n";
        }
        {
            // --- Mean ---
            long double mean
                = std::accumulate(allDifsDist.begin(), allDifsDist.end(), 0.0L)
                / allDifsDist.size();

            // --- Median ---
            std::vector<long double> sorted = allDifsDist;
            std::sort(sorted.begin(), sorted.end());
            long double median;
            long double two = 2.0;
            long double one = 1.0;
            if (sorted.size() % 2 == 0)
                median = (sorted[sorted.size() / two - one]
                          + sorted[sorted.size() / two])
                       / two;
            else
                median = sorted[sorted.size() / two];

            // --- Standard Deviation ---
            long double sumSqDiff = 0.0L;
            for (auto x : allDifsDist) sumSqDiff += (x - mean) * (x - mean);
            long double stddev
                = std::sqrtl(sumSqDiff / (allDifsDist.size() - 1));

            // --- Output ---
            std::cout << "Mean diff (phi (x)) of toi:   " << mean << "\n";
            std::cout << "Median diff (phi (x)) of toi: "
                      << std::setprecision(20) << median << "\n";
            std::cout << "StdDev diff (phi (x)) of toi: " << stddev << "\n";
            std::cout << "min diff (phi (x)) of toi: " << minDiffDist << "\n";
            std::cout << "max diff (phi (x)) of toi: " << maxDiffDist << "\n";
        }
    }

private:
    EigenVector3<T> randomPointOnSphere(std::mt19937& gen, T radius)
    {
        std::uniform_real_distribution<T> dist(-1.0, 1.0);
        EigenVector3<T> point;
        do {
            point = EigenVector3<T>(dist(gen), dist(gen), dist(gen));
        } while (point.norm() < 0.001); // Avoid division by zero
        return point.normalized() * radius;
    }

    std::tuple<EigenVector3<T>, EigenVector3<T>, EigenVector3<T>>
    generateTriangleFacingOrigin(const EigenVector3<T>& center,
                                 std::mt19937& gen)
    {
        std::uniform_real_distribution<T> sizeDist(0.1, 1.0);
        std::uniform_real_distribution<T> angleDist(-0.3, 0.3);

        T triangleSize = sizeDist(gen);

        // Create basis vectors for triangle plane
        EigenVector3<T> toOrigin = -center.normalized();
        EigenVector3<T> tangent1 = toOrigin.cross(EigenVector3<T>(1, 0, 0));
        if (tangent1.norm() < 0.001)
        {
            tangent1 = toOrigin.cross(EigenVector3<T>(0, 1, 0));
        }
        tangent1.normalize();
        EigenVector3<T> tangent2 = toOrigin.cross(tangent1).normalized();

        // Generate triangle vertices
        EigenVector3<T> v0
            = center + tangent1 * triangleSize + tangent2 * triangleSize;
        EigenVector3<T> v1
            = center - tangent1 * triangleSize + tangent2 * triangleSize;
        EigenVector3<T> v2
            = center - tangent1 * triangleSize - tangent2 * triangleSize;

        // Add some random perturbation
        std::uniform_real_distribution<T> perturbDist(-0.1, 0.1);
        v0 += EigenVector3<T>(perturbDist(gen), perturbDist(gen),
                              perturbDist(gen));
        v1 += EigenVector3<T>(perturbDist(gen), perturbDist(gen),
                              perturbDist(gen));
        v2 += EigenVector3<T>(perturbDist(gen), perturbDist(gen),
                              perturbDist(gen));

        return std::make_tuple(v0, v1, v2);
    }

    void classifyResult(CCDStatistics<T>& stats, const T distance, T toi,
                        bool hasIntersection)
    {
        //        T distance = sdf(intersectPoint);
        const T epsilon = T(1e-5);
        //std::cerr << "DIST: " << distance << "\n";

        //Almost at solution (small interpenetration)
        if (std::abs<T>(distance) < epsilon) { stats.almostAtSolution++; }
        //Outside solution, but TOI != tEnd
        else if (distance > epsilon && toi != tEnd)
        {
            stats.outsideButTOINotEnd++;
        }
        //TOI == tEnd and outside solution
        else if (distance > epsilon && toi == tEnd)
        {
            stats.outsideAndTOIEqualsEnd++;
        }
        //Inside solution, but TOI != tStart
        else if (distance < -epsilon && toi > tStart + 0.000001)
        {
            stats.insideButTOINotStart++;
        }
        //TOI == tStart and inside solution
        else if (distance < -epsilon && toi == tStart)
        {
            stats.insideAndTOIEqualsStart++;
        }
        //Other cases
        else { stats.otherCases++; }
    }
};
BOOST_AUTO_TEST_CASE(grid_local_strategy)
{
    {
        using T = double;
        TriangleCCDTester<T> triangleTester;
        triangleTester.runTestsParallel(2000, 5);
        return;
        //        using T = double;
        using D = T;

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
        Eigen::MatrixXd P((Eigen::Index)total, 0);
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

        T dt = 0.00916198;
        EigenVector3<T> angularVelA(0, 0, 0);
        EigenVector3<T> angularVelB(0, 0, 0);
        EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
        EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
        EigenVector3<T> transA(2.0, 0.0667, 0.0667);
        EigenVector3<T> transB(0, 0, 0);
        EigenVector3<T> linearA(-705.707, -214.857, -675.14);
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

        rInfo.A_p0 = EigenVector3<T>(6.41484, 1.53599, 7.5606);
        rInfo.A_p1 = EigenVector3<T>(6.41894, 3.14409, 7.0504);
        rInfo.A_p2 = EigenVector3<T>(7.61829, 2.69284, 5.86892);

        TriangleAtTimeInfo<T> tri = getTriangleAtTime(dt, rInfo);
        EigenVector3<T> p0AtTimeTi = tri.A_p0;
        EigenVector3<T> p1AtTimeTi = tri.A_p1;
        EigenVector3<T> p2AtTimeTi = tri.A_p2;

        EigenVector3<T> contactPoint;
        EigenVector3<T> lastContactPoint;
        EigenVector3<T> normalDummy;
        T penetration;
        optimizeTriangleFW_Working(p0AtTimeTi, p1AtTimeTi, p2AtTimeTi,
                                   *(rInfo.B_sdf), contactPoint, normalDummy,
                                   penetration, 1000);
        T penetrations = valueAtProjection(*(rInfo.B_sdf), contactPoint,
                                           *(rInfo.B_centerTranslation),
                                           *(rInfo.B_centerRotation));
        std::vector<PointPenetration<T>> pens = getPenetrationForTris(
            p0AtTimeTi, p1AtTimeTi, p2AtTimeTi, rInfo, 846);
        T minPenetration = std::numeric_limits<T>::max();
        EigenVector3<T> minPoint;
        for (size_t i = 0; i < pens.size(); ++i)
        {
            PointPenetration<T> currPen = pens[i];
            if (currPen.distToSDF < minPenetration)
            {
                minPoint = currPen.point;
                minPenetration = currPen.distToSDF;
            }
        }
        lastContactPoint = minPoint;
        std::cerr << "Lengt of peentrations: " << pens.size() << "\n";
        std::cerr << "Penetration for FWGD vs 100k point penetration: "
                  << penetrations << " vs " << minPenetration << "\n";
        lastContactPoint = contactPoint;

        EigenVector3<T> firstIntersectPoint;
        std::vector<T> minimizerTimes;
        T toi = FrankWolfeBRENT_BENCHMARK_TIME<T>(
            0.0, 0.01, rInfo, firstIntersectPoint, minimizerTimes);
        T penetrationsBRENT = valueAtProjection(
            *(rInfo.B_sdf), firstIntersectPoint, *(rInfo.B_centerTranslation),
            *(rInfo.B_centerRotation));
        std::cerr << "Ended up with Brent dist vs GT dist: "
                  << penetrationsBRENT << " vs " << minPenetration << "\n";
        std::cerr << "NOTE TOI=" << toi << "\n";

        //if (penetrations <= 0.0) { break; }
    }
}

BOOST_AUTO_TEST_SUITE_END();
