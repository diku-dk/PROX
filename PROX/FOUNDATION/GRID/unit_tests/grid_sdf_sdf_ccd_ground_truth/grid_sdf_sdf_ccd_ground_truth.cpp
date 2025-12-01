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
#include <grid_sdf_sdf_ccd_gradient_descent.hpp>

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

//Note for this function we assume one satatic SDF and one moving towards the static one.
//As soon as it returns true we can stop
template <typename T>
bool getContactInTimeInstance(
    const std::vector<EigenVector3<T>>& points,
    const grid::Grid<T, T>& staticSDF,
    const SDFSDFContact::SingleRigidBodyInfo<T>& movingSDF,
    T currentTimeInstance, std::vector<EigenVector3<T>>& outContactPoints)
{
    for (size_t i = 0; i < points.size(); ++i)
    {
        EigenVector3<T> pointWorldSpace
            = *(movingSDF.A_centerRotation) * points[i]
            + *(movingSDF.A_centerTranslation);
        EigenVector3<T> pointAtNewTime = SDFSDFContact::getVertexPosAtMat(
            *(movingSDF.A_centerTranslation), *(movingSDF.A_linearVel),
            *(movingSDF.A_angularVel), pointWorldSpace, currentTimeInstance);
        /*        T distToSolution = grid::valueAtProjection(
            grid, pointAtNewTime, *(movingSDF.A_centerTranslation),
            *(movingSDF.A_centerRotation));*/
        T distToSolution = grid::value_at_2(staticSDF, pointAtNewTime);
        if (distToSolution <= 0) { outContactPoints.push_back(pointAtNewTime); }
    }
    return outContactPoints.size() > 0;
}

template <typename T> struct StartConfigurations
{
    EigenVector3<T> translationFromZero;
    EigenVector3<T> linearVelocity;
};

template <typename T>
EigenVector3<T> rotateAroundAxis(const EigenVector3<T>& v,
                                 const EigenVector3<T>& axis, T angle)
{
    // axis assumed normalized
    const T ca = std::cos(angle);
    const T sa = std::sin(angle);
    return v * ca + axis.cross(v) * sa + axis * (axis.dot(v)) * (T(1) - ca);
}

template <typename T = double>
std::vector<StartConfigurations<T>> generateStartConfigurations(
    size_t N, T static_radius, T moving_radius,
    T max_start_distance = T(4.5), // default keep <5 so collisions possible
    T min_start_distance
    = T((T)0.0), // if 0, we'll bump it to safe minimum inside the function
    T min_speed = T(0.1), T max_speed = T(50.0),
    T max_perturb_angle_deg = T(10.0), unsigned int rng_seed = 0u)
{
    // Safety margin to avoid exact contact / penetration at start
    const T start_margin = T(1e-3);

    // If user didn't set min_start_distance, enforce non-penetration lower bound
    T d_min_safe = static_radius + moving_radius + start_margin;
    if (min_start_distance < d_min_safe) min_start_distance = d_min_safe;

    const T hard_collision_limit = T(5.0);
    if (max_start_distance > hard_collision_limit)
        max_start_distance = hard_collision_limit;

    // Also ensure that configurations are reachable within max_speed:
    // required_distance_to_contact = d - (static_radius + moving_radius)
    // we need required_distance_to_contact <= max_speed (since t in [0,1])
    // => d <= max_speed + static_radius + moving_radius
    T reachability_cap = max_speed + static_radius + moving_radius;
    if (max_start_distance > reachability_cap)
        max_start_distance = reachability_cap;

    // If after constraints the interval is invalid, try to fix or fail gracefully:
    if (max_start_distance <= min_start_distance)
    {
        // Expand max_start_distance slightly so the loop can produce something
        max_start_distance = min_start_distance + T(1e-2);
    }

    // RNG setup
    std::random_device rd;
    std::mt19937_64 rng(rng_seed == 0u ? rd() : rng_seed);
    std::uniform_real_distribution<T> uni01((T)0.0, (T)1.0);
    std::uniform_real_distribution<T> dist_radius(min_start_distance,
                                                  max_start_distance);
    std::uniform_real_distribution<T> dist_azimuth((T)0.0, (T)2.0 * M_PI);
    std::uniform_real_distribution<T> dist_cos_theta((T)-1.0, (T)1.0);
    std::uniform_real_distribution<T> dist_angle(
        (T)0.0, (T)(max_perturb_angle_deg * M_PI / 180.0));

    std::vector<StartConfigurations<T>> out;
    out.reserve(N);

    for (size_t i = 0; i < N; ++i)
    {
        // sample a point uniformly on spherical shell [rmin, rmax]
        T r = dist_radius(rng);
        // uniform direction on sphere
        T phi = dist_azimuth(rng);
        T cos_theta = dist_cos_theta(rng);
        T sin_theta
            = std::sqrt(std::max((T)0.0, (T)1.0 - cos_theta * cos_theta));

        EigenVector3<T> dir;
        dir.x() = sin_theta * std::cos(phi);
        dir.y() = sin_theta * std::sin(phi);
        dir.z() = cos_theta;

        EigenVector3<T> position = dir * r; // position measured from origin

        // Compute required distance to contact along direct path toward origin:
        // distance from moving center to static surface along line = r - (static_radius + moving_radius)
        T required_distance_to_contact = r - (static_radius + moving_radius);
        if (required_distance_to_contact < (T)0.0)
        {
            // defensive: if by floating rounding the sampled r was too small, push it out
            required_distance_to_contact = (T)0.0;
            r = static_radius + moving_radius + start_margin;
            position = dir * r;
        }

        // Ensure speed is large enough to cover required_distance_to_contact within 1s.
        // Choose speed = max(required_distance_to_contact + small_margin, min_speed), but <= max_speed
        const T speed_safety_margin = (T)0.01;
        T chosen_speed = required_distance_to_contact + speed_safety_margin;
        if (chosen_speed < min_speed) chosen_speed = min_speed;
        if (chosen_speed > max_speed)
        {
            // This should not happen because we clamped max_start_distance earlier,
            // but be defensive: set chosen_speed to max_speed and optionally reduce r so it's reachable.
            chosen_speed = max_speed;
            // Optionally reduce r to be reachable within max_speed (keeps distribution reasonable)
            T reachable_r = max_speed + static_radius + moving_radius;
            if (r > reachable_r)
            {
                r = reachable_r;
                position = dir * r;
                required_distance_to_contact
                    = r - (static_radius + moving_radius);
            }
        }

        // Compute nominal direction toward the static body (origin)
        // vector from moving center to origin is (-position)
        EigenVector3<T> toward_origin = -position;
        T toward_norm = toward_origin.norm();
        if (toward_norm <= std::numeric_limits<T>::epsilon())
        {
            // extremely unlikely: sampled position exactly at origin; nudge slightly outward
            toward_origin = EigenVector3<T>(T(1.0), T(0.0), T(0.0));
            toward_norm = (T)1.0;
            position = EigenVector3<T>(
                (static_radius + moving_radius + start_margin), (T)0.0, (T)0.0);
        }
        toward_origin /= toward_norm; // unit

        // apply angular perturbation up to max_perturb_angle_deg to avoid perfectly centered hits
        T perturb_angle = dist_angle(rng);
        EigenVector3<T> perturbed_direction = toward_origin;
        if (perturb_angle > (T)1e-8)
        {
            // pick a random axis perpendicular to toward_origin: generate random unit vector and orthonormalize
            EigenVector3<T> random_vec;
            // sample random vector uniformly in cube, then orthonormalize
            random_vec.x() = uni01(rng) * 2 - 1;
            random_vec.y() = uni01(rng) * 2 - 1;
            random_vec.z() = uni01(rng) * 2 - 1;
            // make sure not colinear
            if (random_vec.norm() < (T)1e-6)
                random_vec = EigenVector3<T>((T)1.0, (T)0.0, (T)0.0);

            // axis = normalized( random_vec - proj_random_on_toward )
            EigenVector3<T> axis
                = random_vec - toward_origin * (toward_origin.dot(random_vec));
            T axis_norm = axis.norm();
            if (axis_norm < (T)1e-6)
            {
                // fallback axis
                axis = EigenVector3<T>(toward_origin.y(), -toward_origin.x(),
                                       (T)0.0);
                axis_norm = axis.norm();
                if (axis_norm < (T)1e-6)
                    axis = EigenVector3<T>((T)0.0, (T)1.0, (T)0.0),
                    axis_norm = (T)1.0;
            }
            axis /= axis_norm;

            // rotate toward_origin around axis by +perturb_angle (random sign)
            std::uniform_int_distribution<int> sign01(0, 1);
            T signed_angle = perturb_angle * (sign01(rng) ? (T)1.0 : (T)-1.0);
            perturbed_direction
                = rotateAroundAxis<T>(toward_origin, axis, signed_angle);
            // normalize safe
            perturbed_direction.normalize();
        }

        // final linear velocity is in direction 'perturbed_direction' scaled by chosen_speed
        EigenVector3<T> linear_velocity = perturbed_direction * chosen_speed;

        StartConfigurations<T> s;
        s.translationFromZero = position;
        s.linearVelocity = linear_velocity;
        out.push_back(s);
    }

    return out;
}

template <typename T> struct ComputedTimers
{
    std::vector<T> minimizevals;
    T wholeFuncVal;
};

template <typename T> class CCDStatistics2
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

template <typename T>
void classifyResult(CCDStatistics2<T>& stats, const T distance, T toi,
                    bool hasIntersection, T tStart, T tEnd)
{
    //        T distance = sdf(intersectPoint);
    const T epsilon = T(1e-5);
    //std::cerr << "DIST: " << distance << "\n";

    //Almost at solution (small interpenetration)
    if (std::abs<T>(distance) < epsilon) { stats.almostAtSolution++; }
    //Outside solution, but TOI != tEnd
    else if (distance > epsilon && toi != tEnd) { stats.outsideButTOINotEnd++; }
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

template <typename T>
void investigateStartConfigsParallel(
    std::vector<StartConfigurations<T>>& startConfigs,
    const grid::Grid<T, T>& SDFA, const grid::Grid<T, T>& SDFB,
    const std::vector<EigenVector3<T>>& final_points,
    std::vector<SDFSDFContact::SDFVoxel<T>>& finishedVoxelsA,
    std::vector<SDFSDFContact::SDFVoxel<T>>& finishedVoxelsB)
{
    // Global summary containers (will be filled after threads finish)
    T minDiff = std::numeric_limits<T>::max();
    T maxDiff = std::numeric_limits<T>::min();
    std::vector<long double> allDifs;
    allDifs.reserve(startConfigs.size());

    T minDiffDist = std::numeric_limits<T>::max();
    T maxDiffDist = std::numeric_limits<T>::min();
    std::vector<long double> allDifsDist;
    allDifsDist.reserve(startConfigs.size());

    CCDStatistics2<T>
        stats; // will be filled by calling classifyResult serially after threads finish

    const size_t n = startConfigs.size();
    if (n == 0)
    {
        stats.printStatistics();
        return;
    }

    // Decide number of threads
    unsigned int hw = std::thread::hardware_concurrency();
    //size_t num_threads = (hw == 0) ? 4u : std::min<unsigned int>(std::max<unsigned int>(1u, hw), static_cast<unsigned int>(n));
    size_t num_threads = 14;
    // Partition
    size_t chunk = (n + num_threads - 1) / num_threads;

    // Per-thread storage
    struct ClassificationRecord
    {
        T dist;
        T toi;
        bool hasIntersection;
    };

    std::vector<std::thread> workers;
    workers.reserve(num_threads);

    // Each thread writes into its own locals to avoid synchronization during heavy compute
    struct ThreadLocal
    {
        std::vector<long double> diffs;
        std::vector<long double> diffsDist;
        std::vector<ClassificationRecord> classifications;
        T localMinDiff = std::numeric_limits<T>::max();
        T localMaxDiff = std::numeric_limits<T>::min();
        T localMinDiffDist = std::numeric_limits<T>::max();
        T localMaxDiffDist = std::numeric_limits<T>::min();
        T localTimeTaken = T(0);
        std::string stderrBuffer;
    };

    std::vector<ThreadLocal> locals(num_threads);

    // Launch workers
    for (size_t t = 0; t < num_threads; ++t)
    {
        size_t begin = t * chunk;
        size_t end = std::min(begin + chunk, n);

        workers.emplace_back(
            [t, begin, end, &startConfigs, &SDFA, &SDFB, &final_points,
             &finishedVoxelsA, &finishedVoxelsB, &locals]()
            {
                ThreadLocal& local = locals[t];
                local.diffs.reserve((end > begin) ? (end - begin) : 0);
                local.diffsDist.reserve((end > begin) ? (end - begin) : 0);
                local.classifications.reserve((end > begin) ? (end - begin)
                                                            : 0);

                std::ostringstream erross;

                for (size_t i = begin; i < end; ++i)
                {
                    EigenVector3<T> angularVelocity = EigenVector3<T>(0, 0, 0);
                    EigenVector3<T> linearVelocity
                        = startConfigs[i].linearVelocity;
                    EigenVector3<T> centerTranslation
                        = startConfigs[i].translationFromZero;
                    EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
                    SDFSDFContact::SingleRigidBodyInfo<T> bodyInfo;
                    bodyInfo.A_angularVel = &angularVelocity;
                    bodyInfo.A_linearVel = &linearVelocity;
                    bodyInfo.A_centerTranslation = &centerTranslation;
                    bodyInfo.A_centerRotation = &rotA;
                    bodyInfo.sdf = &SDFA;

                    //NOW TODO, MAKE ACTUAL GROUND TRUTH!
                    EigenVector3<T> lastContactPoint;
                    T dt = 0.0;
                    T maxDT = 1.0;
                    T smallStep = 1e-5;
                    T evenSmallerStep = 1e-8;
                    std::vector<EigenVector3<T>> outContactPoints;
                    while (dt <= maxDT)
                    {
                        bool hasPenetrated = getContactInTimeInstance<T>(
                            final_points, SDFB, bodyInfo, dt, outContactPoints);
                        if (hasPenetrated)
                        {
                            T stepBackDT = dt - smallStep;
                            while (stepBackDT <= dt)
                            {
                                bool hasPenetrated2
                                    = getContactInTimeInstance<T>(
                                        final_points, SDFB, bodyInfo, dt,
                                        outContactPoints);
                                if (hasPenetrated2)
                                {
                                    dt = stepBackDT;
                                    break;
                                }
                                stepBackDT += evenSmallerStep;
                            }

                            break;
                        }

                        dt += smallStep;
                    }
                    dt = std::min<T>(maxDT, dt);
                    //For now select one point only.
                    if (dt != maxDT && false)
                    {
                        lastContactPoint = outContactPoints[0];
                    }
                    else
                    {
                        lastContactPoint
                            = EigenVector3<T>(-10000.0, -10000.0, -10000.0);
                    }

                    //Now run our algo's ground truth!
                    EigenVector3<T> angularVelB(0, 0, 0);
                    EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
                    EigenVector3<T> transB(0.0, 0.0, 0.0);
                    EigenVector3<T> linearB(0, 0, 0);
                    SDFSDFContact::SingleRigidBodyInfo<T> rInfoA = bodyInfo;

                    SDFSDFContact::SingleRigidBodyInfo<T> rInfoB;
                    rInfoB.A_angularVel = &angularVelB;
                    rInfoB.A_centerRotation = &rotB;
                    rInfoB.A_centerTranslation = &transB;
                    rInfoB.A_linearVel = &linearB;
                    rInfoB.sdf = &SDFB;

                    std::vector<EigenVector3<T>> outContacts;
                    auto gss_start = std::chrono::high_resolution_clock::now();

                    T toi = SDFSDFContact::getSDFSDFTOI(
                        finishedVoxelsA, finishedVoxelsB, rInfoA, rInfoB,
                        T(0.0), T(1.0), outContacts);
                    auto gss_end = std::chrono::high_resolution_clock::now();
                    auto gss_us
                        = std::chrono::duration_cast<std::chrono::nanoseconds>(
                              gss_end - gss_start)
                              .count();
                    T totalAlgorithmTime = T(gss_us);

                    std::cerr << "outContacts: " << outContacts[0];
                    toi = std::min<T>(toi, maxDT);
                    EigenVector3<T> firstIntersectPoint;
                    if (toi != maxDT && false)
                    {
                        firstIntersectPoint = outContacts[0];
                    }
                    else
                    {
                        firstIntersectPoint
                            = EigenVector3<T>(-10000.0, -10000.0, -10000.0);
                    }

                    // Buffer some debug output — we'll flush it thread-safely later
                    erross << "GOT DT VS TOI: " << dt << " vs " << toi << "\n";
                    erross << "Center translation for SDF A: "
                           << *(rInfoA.A_centerTranslation) << "\n";
                    erross << "Linear velocity for SDF A: "
                           << *(rInfoA.A_linearVel) << "\n";
                    erross << "END translation: "
                           << (*(rInfoA.A_centerTranslation)
                               + *(rInfoA.A_linearVel) * toi)
                           << "\n";
                    erross << "END translation GT: "
                           << (*(rInfoA.A_centerTranslation)
                               + *(rInfoA.A_linearVel) * dt)
                           << "\n\n";

                    T dist
                        = grid::valueAtProjection(SDFB, firstIntersectPoint,
                                                  *(rInfoB.A_centerTranslation),
                                                  *(rInfoB.A_centerRotation));
                    T distGT = grid::valueAtProjection(
                        SDFB, lastContactPoint, *(rInfoB.A_centerTranslation),
                        *(rInfoB.A_centerRotation));

                    T diffDistance = dist - distGT;

                    if (dist < 0.0 && dist < -T(1e-3))
                    {
                        erross << "XTI: " << firstIntersectPoint << "\n";
                        erross << "TOI: " << toi << "\n";
                        erross << "DIST: " << dist << "\n";
                    }

                    bool hasIntersection = (toi >= 0.01);
                    T diff = toi - dt;

                    if (diffDistance <= -1e-5)
                    {
                        erross << "=========================================\n";
                        erross << "Note diffDistance is" << diffDistance
                               << ", diffGT is " << distGT << " and diff is "
                               << dist << "\n";
                        erross << "=========================================\n";
                    }
                    if (diffDistance >= 1e-2)
                    {
                        erross << "=========================================\n";
                        erross << "Note diffDistance is" << diffDistance
                               << ", diffGT is " << distGT << " and diff is "
                               << dist << "\n";
                        erross << "=========================================\n";
                    }

                    if (toi > dt + 1e-5)
                    {
                        erross << "=========================================\n";
                        erross << "WE GOT A TOI that is larger thann DT! (1e5 "
                                  "prec)\n";
                        erross << "=========================================\n";
                    }
                    else if (toi > dt + 1e-8)
                    {
                        erross << "WE GOT A TOI that is larger thann DT! (1e8 "
                                  "prec)\n";
                    }
                    else if (toi > dt + 1e-11)
                    {
                        erross << "WE GOT A TOI that is larger thann DT! (no "
                                  "prec)\n";
                    }

                    // Save results locally
                    local.diffs.push_back(static_cast<long double>(diff));
                    local.diffsDist.push_back(
                        static_cast<long double>(diffDistance));
                    local.classifications.push_back(
                        ClassificationRecord{dist, toi, hasIntersection});

                    // Update local mins/maxs
                    if (diff < local.localMinDiff) local.localMinDiff = diff;
                    if (diff > local.localMaxDiff) local.localMaxDiff = diff;
                    if (diffDistance < local.localMinDiffDist)
                        local.localMinDiffDist = diffDistance;
                    if (diffDistance > local.localMaxDiffDist)
                        local.localMaxDiffDist = diffDistance;
                    local.localTimeTaken += totalAlgorithmTime;

                    // --- end per-config code ---
                } // end for each i

                // store buffered stderr
                local.stderrBuffer = erross.str();
            });
    } // end launch threads

    // join
    for (auto& th : workers) th.join();

    // Merge per-thread results into globals
    T allTime = 0.0;
    for (const auto& local : locals)
    {
        // append diffs
        allDifs.insert(allDifs.end(), local.diffs.begin(), local.diffs.end());
        allDifsDist.insert(allDifsDist.end(), local.diffsDist.begin(),
                           local.diffsDist.end());

        // combine min/max
        if (local.localMinDiff < minDiff) minDiff = local.localMinDiff;
        if (local.localMaxDiff > maxDiff) maxDiff = local.localMaxDiff;
        if (local.localMinDiffDist < minDiffDist)
            minDiffDist = local.localMinDiffDist;
        if (local.localMaxDiffDist > maxDiffDist)
            maxDiffDist = local.localMaxDiffDist;
        allTime += local.localTimeTaken;

        // print buffered stderr in sequence to avoid interleaving
        if (!local.stderrBuffer.empty()) { std::cerr << local.stderrBuffer; }
    }

    // Now call classifyResult serially for every collected classification record (this updates stats)
    for (const auto& local : locals)
    {
        for (const auto& rec : local.classifications)
        {
            classifyResult<T>(stats, rec.dist, rec.toi, rec.hasIntersection,
                              T(0.0), T(1.0));
        }
    }

    // Print stats gathered by classifyResult
    stats.printStatistics();
    std::cerr << "The gradient descent function took "
              << allTime / T(startConfigs.size()) << " ns.\n";

    // Compute aggregated summary stats for allDifs (same as original)
    {
        // --- Mean ---
        long double mean = std::accumulate(allDifs.begin(), allDifs.end(), 0.0L)
                         / allDifs.size();

        // --- Median ---
        std::vector<long double> sorted = allDifs;
        std::sort(sorted.begin(), sorted.end());
        long double median;
        if (sorted.size() % 2 == 0)
            median = (sorted[sorted.size() / 2 - 1] + sorted[sorted.size() / 2])
                   / 2.0L;
        else
            median = sorted[sorted.size() / 2];

        // --- Standard Deviation ---
        long double sumSqDiff = 0.0L;
        for (auto x : allDifs) sumSqDiff += (x - mean) * (x - mean);
        long double stddev = std::sqrtl(sumSqDiff / (allDifs.size() - 1));

        // --- Output ---
        std::cout << "Mean diff of toi:   " << mean << "\n";
        std::cout << "Median diff of toi: " << std::setprecision(20) << median
                  << "\n";
        std::cout << "StdDev diff of toi: " << stddev << "\n";
        std::cout << "min diff of toi: " << minDiff << "\n";
        std::cout << "max diff of toi: " << maxDiff << "\n";
    }

    // Compute aggregated summary stats for allDifsDist (same as original)
    {
        long double mean
            = std::accumulate(allDifsDist.begin(), allDifsDist.end(), 0.0L)
            / allDifsDist.size();

        std::vector<long double> sorted = allDifsDist;
        std::sort(sorted.begin(), sorted.end());
        long double median;
        if (sorted.size() % 2 == 0)
            median = (sorted[sorted.size() / 2 - 1] + sorted[sorted.size() / 2])
                   / 2.0L;
        else
            median = sorted[sorted.size() / 2];

        long double sumSqDiff = 0.0L;
        for (auto x : allDifsDist) sumSqDiff += (x - mean) * (x - mean);
        long double stddev = std::sqrtl(sumSqDiff / (allDifsDist.size() - 1));

        std::cout << "Mean diff (phi (x)) of toi:   " << mean << "\n";
        std::cout << "Median diff (phi (x)) of toi: " << std::setprecision(20)
                  << median << "\n";
        std::cout << "StdDev diff (phi (x)) of toi: " << stddev << "\n";
        std::cout << "min diff (phi (x)) of toi: " << minDiffDist << "\n";
        std::cout << "max diff (phi (x)) of toi: " << maxDiffDist << "\n";
    }
}

template <typename T>
void investigateStartConfigs(
    std::vector<StartConfigurations<T>>& startConfigs,
    const grid::Grid<T, T>& SDFA, const grid::Grid<T, T>& SDFB,
    const std::vector<EigenVector3<T>>& final_points,
    std::vector<SDFSDFContact::SDFVoxel<T>>& finishedVoxelsA,
    std::vector<SDFSDFContact::SDFVoxel<T>>& finishedVoxelsB)
{

    T minDiff = std::numeric_limits<T>::max();
    T maxDiff = std::numeric_limits<T>::min();
    std::vector<long double> allDifs(startConfigs.size());

    T minDiffDist = std::numeric_limits<T>::max();
    T maxDiffDist = std::numeric_limits<T>::min();
    std::vector<long double> allDifsDist(startConfigs.size());
    CCDStatistics2<T> stats;

    for (size_t i = 0; i < startConfigs.size(); ++i)
    {
        /*EigenVector3<T> angularVelocity(0.0, 0.0, 0.0);
                EigenVector3<T> linearVelocity(0.0, 0.0, 0.0);
                EigenVector3<T> centerTranslation(0.0, 0.0, 0.0);*/
        EigenVector3<T> angularVelocity = EigenVector3<T>(0, 0, 0);
        EigenVector3<T> linearVelocity = startConfigs[i].linearVelocity;
        EigenVector3<T> centerTranslation = startConfigs[i].translationFromZero;
        EigenQuaternion<T> rotA = EigenQuaternion<T>::Identity();
        SDFSDFContact::SingleRigidBodyInfo<T> bodyInfo;
        bodyInfo.A_angularVel = &angularVelocity;
        bodyInfo.A_linearVel = &linearVelocity;
        bodyInfo.A_centerTranslation = &centerTranslation;
        bodyInfo.A_centerRotation = &rotA;
        bodyInfo.sdf = &SDFA;

        //NOW TODO, MAKE ACTUAL GROUND TRUTH!
        EigenVector3<T> lastContactPoint;
        T dt = 0.0;
        T maxDT = 1.0;
        T smallStep = 1e-5;
        T evenSmallerStep = 1e-8;
        std::vector<EigenVector3<T>> outContactPoints;
        while (dt <= maxDT)
        {
            bool hasPenetrated = getContactInTimeInstance<T>(
                final_points, SDFB, bodyInfo, dt, outContactPoints);
            if (hasPenetrated)
            {
                T stepBackDT = dt - smallStep;
                while (stepBackDT <= dt)
                {
                    bool hasPenetrated2 = getContactInTimeInstance<T>(
                        final_points, SDFB, bodyInfo, dt, outContactPoints);
                    if (hasPenetrated2)
                    {
                        dt = stepBackDT;
                        break;
                    }
                    stepBackDT += evenSmallerStep;
                }

                break;
            }

            dt += smallStep;
        }

        dt = std::min<T>(maxDT, dt);
        //For now select one point only.
        if (dt != maxDT && false) { lastContactPoint = outContactPoints[0]; }
        else
        {
            lastContactPoint = EigenVector3<T>(-10000.0, -10000.0, -10000.0);
        }

        //Now run our algo's ground truth!
        EigenVector3<T> angularVelB(0, 0, 0);
        EigenQuaternion<T> rotB = EigenQuaternion<T>::Identity();
        EigenVector3<T> transB(0.0, 0.0, 0.0);
        EigenVector3<T> linearB(0, 0, 0);
        SDFSDFContact::SingleRigidBodyInfo<T> rInfoA = bodyInfo;

        SDFSDFContact::SingleRigidBodyInfo<T> rInfoB;
        rInfoB.A_angularVel = &angularVelB;
        rInfoB.A_centerRotation = &rotB;
        rInfoB.A_centerTranslation = &transB;
        rInfoB.A_linearVel = &linearB;
        rInfoB.sdf = &SDFB;

        std::vector<EigenVector3<T>> outContacts;
        T toi = SDFSDFContact::getSDFSDFTOI(finishedVoxelsA, finishedVoxelsB,
                                            rInfoA, rInfoB, T(0.0), T(1.0),
                                            outContacts);

        toi = std::min<T>(toi, maxDT);
        EigenVector3<T> firstIntersectPoint;
        if (toi != maxDT && false) { firstIntersectPoint = outContacts[0]; }
        else
        {
            firstIntersectPoint = EigenVector3<T>(-10000.0, -10000.0, -10000.0);
        }

        std::cerr << "OUTCONTACTS SIZE" << outContacts.size();
        std::cerr << "\nFIRST OUT CIONTACTS:" << outContacts[0] << "\n";
        std::cerr << "GOT DT VS TOI: " << dt << " vs " << toi << "\n";
        std::cerr << "Center translation for SDF A: "
                  << *(rInfoA.A_centerTranslation) << "\n";
        std::cerr << "Linear velocity for SDF A: " << *(rInfoA.A_linearVel)
                  << "\n";
        std::cerr << "END translation: "
                  << (*(rInfoA.A_centerTranslation)
                      + *(rInfoA.A_linearVel) * toi)
                  << "\n";

        std::cerr << "END translation GT: "
                  << (*(rInfoA.A_centerTranslation)
                      + *(rInfoA.A_linearVel) * dt)
                  << "\n";
        std::cerr << "\n\n";

        T dist = grid::valueAtProjection(SDFB, firstIntersectPoint,
                                         *(rInfoB.A_centerTranslation),
                                         *(rInfoB.A_centerRotation));
        T distGT = grid::valueAtProjection(SDFB, lastContactPoint,
                                           *(rInfoB.A_centerTranslation),
                                           *(rInfoB.A_centerRotation));

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
            std::cerr << "========================================="
                         "====\n";
            std::cerr << "Note diffDistance is" << diffDistance
                      << ", diffGT is " << distGT << " and diff is " << dist
                      << "\n";
            std::cerr << "========================================="
                         "====\n";
        }
        if (diffDistance >= 1e-2)
        {
            std::cerr << "========================================="
                         "====\n";

            std::cerr << "Note diffDistance is" << diffDistance
                      << ", diffGT is " << distGT << " and diff is " << dist
                      << "\n";
            std::cerr << "========================================="
                         "====\n";
        }

        if (toi > dt + 1e-5)
        {
            std::cerr << "========================================="
                         "====\n";
            std::cerr << "WE GOT A TOI that is larger thann DT! "
                         "(1e5 prec)\n";
            std::cerr << "========================================="
                         "====\n";
        }
        else if (toi > dt + 1e-8)
        {
            std::cerr << "WE GOT A TOI that is larger thann DT! "
                         "(1e8 prec)\n";
        }
        else if (toi > dt + 1e-11)
        {
            std::cerr << "WE GOT A TOI that is larger thann DT! "
                         "(no prec)\n";
        }

        minDiff = std::min<T>(minDiff, diff);
        maxDiff = std::max<T>(maxDiff, diff);
        allDifs.push_back(static_cast<long double>(diff));

        minDiffDist = std::min<T>(minDiffDist, diffDistance);
        maxDiffDist = std::max<T>(maxDiffDist, diffDistance);
        allDifsDist.push_back(static_cast<long double>(diffDistance));
        classifyResult<T>(stats, dist, toi, hasIntersection, T(0.0), T(1.0));
    }
    stats.printStatistics();

    {
        // --- Mean ---
        long double mean = std::accumulate(allDifs.begin(), allDifs.end(), 0.0L)
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
        std::cout << "Median diff of toi: " << std::setprecision(20) << median
                  << "\n";
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
        long double stddev = std::sqrtl(sumSqDiff / (allDifsDist.size() - 1));

        // --- Output ---
        std::cout << "Mean diff (phi (x)) of toi:   " << mean << "\n";
        std::cout << "Median diff (phi (x)) of toi: " << std::setprecision(20)
                  << median << "\n";
        std::cout << "StdDev diff (phi (x)) of toi: " << stddev << "\n";
        std::cout << "min diff (phi (x)) of toi: " << minDiffDist << "\n";
        std::cout << "max diff (phi (x)) of toi: " << maxDiffDist << "\n";
    }
}

BOOST_AUTO_TEST_CASE(grid_local_strategy)
{
    {
        using D = float;
        using T = float;
        {
            grid::Grid<T, T> SDFA;
            std::vector<SDFSDFContact::SDFVoxel<T>> finishedVoxelsA;
            std::vector<SDFSDFContact::SDFVoxel<T>> voxelsA;
            makeSDF(SDFA, "blender_star.obj", finishedVoxelsA, voxelsA, 64, 8);
            std::cerr << "Finished filtering SDF 1/2!" << "\n";
            //This is static for tihs case!
            grid::Grid<T, T> SDFB;
            std::vector<SDFSDFContact::SDFVoxel<T>> voxelsB;
            std::vector<SDFSDFContact::SDFVoxel<T>> finishedVoxelsB;
            makeSDF(SDFB, "blender_star.obj", finishedVoxelsB, voxelsB, 64, 8);
            std::cerr << "Finished filtering SDF 2/2!" << "\n";
            std::vector<Eigen::Matrix<T, 3, 1>>
                final_points; // collects points to push

            const T eps = static_cast<T>(1e-12); // "almost zero" threshold
            const std::array<std::array<int, 3>, 6> offs = {
                {{{1, 0, 0}},
                 {{-1, 0, 0}},
                 {{0, 1, 0}},
                 {{0, -1, 0}},
                 {{0, 0, 1}},
                 {{0, 0, -1}}}
            };

            for (size_t k = 0; k < SDFA.K(); ++k)
            {
                for (size_t j = 0; j < SDFA.J(); ++j)
                {
                    for (size_t i = 0; i < SDFA.I(); ++i)
                    {
                        Eigen::Matrix<size_t, 3, 1> idx(i, j, k);
                        Eigen::Matrix<T, 3, 1> p;
                        grid::node_position(SDFA, idx, p);
                        const T val = grid::value_at_2(SDFA, p);
                        //const T val = p(0);

                        // 1) nearly zero -> push and skip neighbour checks
                        if (std::abs(val) <= eps)
                        {
                            final_points.push_back(p);
                            continue;
                        }

                        // use signed arithmetic for neighbour indices
                        const int ii = static_cast<int>(i);
                        const int jj = static_cast<int>(j);
                        const int kk = static_cast<int>(k);
                        const int Ii = static_cast<int>(SDFA.I());
                        const int Jj = static_cast<int>(SDFA.J());
                        const int Kk = static_cast<int>(SDFA.K());

                        // 2 & 3) check neighbours for sign changes
                        for (const auto& d : offs)
                        {
                            const int ni = ii + d[0];
                            const int nj = jj + d[1];
                            const int nk = kk + d[2];

                            if (ni < 0 || ni >= Ii || nj < 0 || nj >= Jj
                                || nk < 0 || nk >= Kk)
                            {
                                continue; // neighbour outside grid
                            }

                            Eigen::Matrix<size_t, 3, 1> nidx(
                                static_cast<size_t>(ni),
                                static_cast<size_t>(nj),
                                static_cast<size_t>(nk));
                            Eigen::Matrix<T, 3, 1> np;
                            grid::node_position(SDFA, nidx, np);
                            const T nval = grid::value_at_2(SDFA, np);
                            //const T nval = np(0);

                            // if current > 0 and neighbor < 0, push
                            if (val > static_cast<T>(0)
                                && nval < static_cast<T>(0))
                            {
                                final_points.push_back(p);
                                break;
                            }

                            // if current < 0 and neighbor > 0, push
                            if (val < static_cast<T>(0)
                                && nval > static_cast<T>(0))
                            {
                                final_points.push_back(p);
                                break;
                            }
                        } // neighbours
                    }
                }
            }
            //Now project all points to isosurface such that we are sure they cover the isosurface
            for (size_t i = 0; i < final_points.size(); ++i)
            {
                final_points[i] = final_points[i]
                                - grid::value_at_2(SDFA, final_points[i])
                                      * grid::computeGradient_Working(
                                          final_points[i], SDFA);
            }

            std::cerr << "SIZE: " << final_points.size() << "\n";

            std::vector<StartConfigurations<T>> startConfigs
                = generateStartConfigurations<T>(50, 1.0, 1.0, 6.0, 3.0, 15.0,
                                                 50);
            investigateStartConfigsParallel<T>(startConfigs, SDFA, SDFB,
                                               final_points, finishedVoxelsA,
                                               finishedVoxelsB);
        }
    }
}
BOOST_AUTO_TEST_SUITE_END();
