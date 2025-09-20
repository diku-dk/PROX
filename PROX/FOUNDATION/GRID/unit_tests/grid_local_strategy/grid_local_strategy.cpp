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

#include "igl/signed_distance.h"
#include "igl/marching_cubes.h"
#include "igl/writeOBJ.h"
//#include <libigl/signed_distance.h>
//#include <libigl/sphere.h>
//#include <libigl/marching_cubes.h>
//#include <libigl/writeOBJ.h>

#include <igl/read_triangle_mesh.h>

BOOST_AUTO_TEST_SUITE(grid);

void write_obj(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
               const std::string& filename)
{
    std::ofstream out(filename);
    if (!out)
    {
        std::cerr << "Failed to open " << filename << "\n";
        return;
    }

    // Write vertices (OBJ is 1-based indexing)
    for (int i = 0; i < V.rows(); ++i)
    {
        out << "v " << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << "\n";
    }

    // Optionally write faces. Faces expect 1-based indices.
    for (int i = 0; i < F.rows(); ++i)
    {
        out << "f " << (F(i, 0) + 1) << " " << (F(i, 1) + 1) << " "
            << (F(i, 2) + 1) << "\n";
    }

    out.close();
}

// Write a single triangle (3 vertices) as an OBJ (3 vertices + 1 face)
template <typename T>
void write_triangle_obj(const Eigen::Matrix<T, 3, 1>& p,
                        const Eigen::Matrix<T, 3, 1>& q,
                        const Eigen::Matrix<T, 3, 1>& r,
                        const std::string& filename)
{
    std::ofstream out(filename);
    if (!out)
    {
        std::cerr << "Failed to open " << filename << "\n";
        return;
    }

    out << "v " << p.x() << " " << p.y() << " " << p.z() << "\n";
    out << "v " << q.x() << " " << q.y() << " " << q.z() << "\n";
    out << "v " << r.x() << " " << r.y() << " " << r.z() << "\n";
    out << "f 1 2 3\n";
    out.close();
}

// Write a single contact point to a simple .xyz file (one line: x y z)
template <typename T>
void write_contact_xyz(const Eigen::Matrix<T, 3, 1>& contact,
                       const std::string& filename)
{
    std::ofstream out(filename);
    if (!out)
    {
        std::cerr << "Failed to open " << filename << "\n";
        return;
    }
    out << contact.x() << " " << contact.y() << " " << contact.z() << "\n";
    out.close();
}

template <typename T> T sphereSdf(const Eigen::Matrix<T, 3, 1>& p, const Eigen::Matrix<T, 3, 1>& c, T r)
{
    return (p - c).norm() - r;
}

template <typename T>
T coneSdf(const Eigen::Matrix<T, 3, 1>& p, const Eigen::Matrix<T, 3, 1>& c,
          const Eigen::Matrix<T, 3, 1>& v, T h, T r)
{
    // Normalize the axis vector v
    Eigen::Matrix<T, 3, 1> axis = v.normalized();

    // Vector from the apex to the point
    Eigen::Matrix<T, 3, 1> dp = p - c;

    // Project the vector dp onto the cone axis
    T d = dp.dot(axis);

    // Compute the perpendicular distance to the cone axis
    Eigen::Matrix<T, 3, 1> radial = dp - d * axis;
    T radialDistance = radial.norm();

    // Compute the radius at the height d
    T coneRadiusAtHeight = (d / h) * r;

    // If the point is below the apex or above the cone's top, it's outside the cone
    if (d < 0) return dp.norm(); // Distance to the apex (cone vertex)
    if (d > h)
        return (p - (c + h * axis)).norm(); // Distance to the cone's base

    // Otherwise, return the signed distance to the cone's side
    return radialDistance - coneRadiusAtHeight;
}

template <typename D, typename T>
Eigen::Matrix<T, 3, 1> computeGradient(const Eigen::Matrix<T, 3, 1>& p,
                                       const grid::Grid<D, T>& grid)
{
    float h = 0.0075;
    float dx = (grid::value_at_2(
                    grid, Eigen::Matrix<T, 3, 1>(p.x() + h, p.y(), p.z()))
                - grid::value_at_2(
                    grid, Eigen::Matrix<T, 3, 1>(p.x() - h, p.y(), p.z())))
             / (2 * h);
    float dy = (grid::value_at_2(
                    grid, Eigen::Matrix<T, 3, 1>(p.x(), p.y() + h, p.z()))
                - grid::value_at_2(
                    grid, Eigen::Matrix<T, 3, 1>(p.x(), p.y() - h, p.z())))
             / (2 * h);
    float dz = (grid::value_at_2(
                    grid, Eigen::Matrix<T, 3, 1>(p.x(), p.y(), p.z() + h))
                - grid::value_at_2(
                    grid, Eigen::Matrix<T, 3, 1>(p.x(), p.y(), p.z() - h)))
             / (2 * h);
    return Eigen::Matrix<T, 3, 1>(dx, dy, dz).normalized();
}

// Frank-Wolfe optimization for triangle-SDF collision detection
template <typename D, typename T>
bool optimizeTriangleFWBad(const Eigen::Matrix<T, 3, 1>& p,
                           const Eigen::Matrix<T, 3, 1>& q,
                           const Eigen::Matrix<T, 3, 1>& r,
                           const grid::Grid<D, T>& cone,
                           Eigen::Matrix<T, 3, 1>& contactPoint,
                           Eigen::Matrix<T, 3, 1>& normal, T& penetration,
                           int maxIterations = 20)
{
    // Initialize to triangle centroid
    Eigen::Matrix<T, 3, 1> x = (p + q + r) / 3.0f;

    for (int i = 0; i < maxIterations; i++)
    {
        // Compute gradient at current position
        Eigen::Matrix<T, 3, 1> gradPhi = computeGradient<D, T>(x, cone);

        // Find the vertex that minimizes s^T * ∇φ(x_i)
        T pDot = p.dot(gradPhi);
        T qDot = q.dot(gradPhi);
        T rDot = r.dot(gradPhi);

        Eigen::Matrix<T, 3, 1> s;
        if (pDot <= qDot && pDot <= rDot) { s = p; }
        else if (qDot <= pDot && qDot <= rDot) { s = q; }
        else { s = r; }

        // Frank-Wolfe step size
        T alpha = 2.0f / (i + 2.0f);

        // Update position
        x = x + alpha * (s - x);
    }

    // Calculate final results
    contactPoint = x;
    std::cerr << "CONTACT POINT" << contactPoint << "\n";

    penetration = grid::value_at_2(cone, contactPoint);
    normal = computeGradient<D, T>(contactPoint, cone);

    // Return true if penetration is negative (inside the object)
    return penetration < 0.0f;
}

template <typename T> void project_onto_simplex(Eigen::Matrix<T, 3, 1>& v)
{
    // Sort the vector in descending order
    Eigen::Matrix<T, 3, 1> u = v;
    std::sort(u.data(), u.data() + 3, std::greater<T>());

    T sum = 0.0;
    int rho = -1;
    for (int j = 0; j < 3; j++)
    {
        sum += u(j);
        if (u(j) + (1.0f - sum) / (j + 1) > 0) { rho = j; }
    }

    T lambda = (1.0f - u.head(rho + 1).sum()) / (rho + 1);
    v = (v.array() + lambda).cwiseMax(0.0f);
}

template <typename D, typename T>
Eigen::Matrix<T, 3, 1> computeGradient2(const Eigen::Matrix<T, 3, 1>& p,
                                        const grid::Grid<D, T>& grid)
{
    // Calculate grid cell size based on resolution and bounds
    Eigen::Matrix<T, 3, 1> diff = (grid.m_max - grid.m_min);
    Eigen::Matrix<T, 3, 1> cell_size = Eigen::Matrix<T, 3, 1>(
        diff.x() / (grid.m_nodes.x() - 1), diff.y() / (grid.m_nodes.y() - 1),
        diff.z() / (grid.m_nodes.z() - 1));

    // Use cell size for finite differences
    T hx = cell_size.x();
    T hy = cell_size.y();
    T hz = cell_size.z();

    // Central difference for gradient approximation
    T dx = (grid::value_at_2(grid,
                             Eigen::Matrix<T, 3, 1>(p.x() + hx, p.y(), p.z()))
            - grid::value_at_2(
                grid, Eigen::Matrix<T, 3, 1>(p.x() - hx, p.y(), p.z())))
         / (2 * hx);

    T dy = (grid::value_at_2(grid,
                             Eigen::Matrix<T, 3, 1>(p.x(), p.y() + hy, p.z()))
            - grid::value_at_2(
                grid, Eigen::Matrix<T, 3, 1>(p.x(), p.y() - hy, p.z())))
         / (2 * hy);

    T dz = (grid::value_at_2(grid,
                             Eigen::Matrix<T, 3, 1>(p.x(), p.y(), p.z() + hz))
            - grid::value_at_2(
                grid, Eigen::Matrix<T, 3, 1>(p.x(), p.y(), p.z() - hz)))
         / (2 * hz);

    return Eigen::Matrix<T, 3, 1>(dx, dy, dz);
}

// Frank-Wolfe optimization for triangle-SDF collision detection
template <typename D, typename T>
bool optimizeTriangleFW2(const Eigen::Matrix<T, 3, 1>& p,
                         const Eigen::Matrix<T, 3, 1>& q,
                         const Eigen::Matrix<T, 3, 1>& r,
                         const grid::Grid<D, T>& cone,
                         Eigen::Matrix<T, 3, 1>& contactPoint,
                         Eigen::Matrix<T, 3, 1>& normal, T& penetration,
                         int maxIterations = 32, T tolerance = 1e-4f)
{
    // Better initialization: evaluate at vertices and choose the one with smallest SDF value
    T phi_p = grid::value_at_2<D, T>(cone, p);
    T phi_q = grid::value_at_2<D, T>(cone, q);
    T phi_r = grid::value_at_2<D, T>(cone, r);

    Eigen::Matrix<T, 3, 1> x;
    if (phi_p <= phi_q && phi_p <= phi_r) { x = p; }
    else if (phi_q <= phi_p && phi_q <= phi_r) { x = q; }
    else { x = r; }

    T prev_obj_value = std::numeric_limits<T>::max();

    for (int i = 0; i < maxIterations; i++)
    {
        // Compute gradient at current position
        Eigen::Matrix<T, 3, 1> gradPhi = computeGradient2<D, T>(x, cone);

        // Find the vertex that minimizes s^T * ∇φ(x_i)
        T pDot = p.dot(gradPhi);
        T qDot = q.dot(gradPhi);
        T rDot = r.dot(gradPhi);

        Eigen::Matrix<T, 3, 1> s;
        if (pDot <= qDot && pDot <= rDot) { s = p; }
        else if (qDot <= pDot && qDot <= rDot) { s = q; }
        else { s = r; }

        // Frank-Wolfe step size (decreasing as suggested in the paper)
        T alpha = static_cast<T>(2.0 / (i + 2.0));

        // Update position
        Eigen::Matrix<T, 3, 1> new_x = x + alpha * (s - x);

        // Check for convergence
        T obj_value = grid::value_at_2<D, T>(cone, new_x);
        T improvement = std::abs(prev_obj_value - obj_value);

        if (improvement < tolerance)
        {
            x = new_x;
            break;
        }

        x = new_x;
        prev_obj_value = obj_value;
    }

    // Calculate final results
    contactPoint = x;
    penetration = grid::value_at_2<D, T>(cone, contactPoint);
    normal = computeGradient2<D, T>(contactPoint, cone);

    // Normalize the normal (but only if it's not zero)
    T norm = normal.norm();
    if (norm > std::numeric_limits<T>::epsilon()) { normal /= norm; }

    // Return true if penetration is negative (inside the object)
    return penetration < static_cast<T>(0);
}

template <typename D, typename T>
bool optimizeTriangleFW3(const Eigen::Matrix<T, 3, 1>& p,
                         const Eigen::Matrix<T, 3, 1>& q,
                         const Eigen::Matrix<T, 3, 1>& r,
                         const grid::Grid<D, T>& sdf,
                         Eigen::Matrix<T, 3, 1>& contactPoint,
                         Eigen::Matrix<T, 3, 1>& normal, T& penetration,
                         size_t maxIterations = 320)
{
    // Better initialization: evaluate at vertices and choose the one with smallest SDF value
    T phi_p = grid::value_at_2<D, T>(sdf, p);
    T phi_q = grid::value_at_2<D, T>(sdf, q);
    T phi_r = grid::value_at_2<D, T>(sdf, r);

    Eigen::Matrix<T, 3, 1> x;
    if (phi_p <= phi_q && phi_p <= phi_r) { x = p; }
    else if (phi_q <= phi_p && phi_q <= phi_r) { x = q; }
    else { x = r; }

    for (size_t i = 0; i < maxIterations; ++i)
    {
        Eigen::Matrix<T, 3, 1> gradient = computeGradient2(x, sdf);
        Eigen::Matrix<T, 3, 1> gradientTransposed = gradient.transpose();
        T Lp = gradientTransposed.dot(p);
        T Lq = gradientTransposed.dot(q);
        T Lr = gradientTransposed.dot(r);
        Eigen::Matrix<T, 3, 1> si;
        if (Lp <= Lq && Lp <= Lr) { si = p; }
        else if (Lq <= Lp && Lq <= Lr) { si = q; }
        else { si = r; }

        //Eigen::Matrix<T, 3, 1> sitmp = (si.transpose().eval()).dot(gradient);
        T alpha = T(2) / (T(i) + T(2));
        //xi+1 = xi...
        x = x + alpha * (si - x);
    }

    contactPoint = x;
    normal = computeGradient2(x, sdf);
    penetration = grid::value_at_2<D, T>(sdf, x);
    return penetration <= T(0);
}

// Projected Gradient Descent optimization for triangle-SDF collision detection
template <typename D, typename T>
bool optimizeTrianglePGD(const Eigen::Matrix<T, 3, 1>& p,
                         const Eigen::Matrix<T, 3, 1>& q,
                         const Eigen::Matrix<T, 3, 1>& r,
                         const grid::Grid<D, T>& sdfGrid,
                         Eigen::Matrix<T, 3, 1>& contactPoint,
                         Eigen::Matrix<T, 3, 1>& normal, T& penetration,
                         int maxIterations = 1000, T initialAlpha = 0.001f)
{
    // Initialize to the vertex with minimum SDF value
    T phi_p = grid::value_at_2<D, T>(sdfGrid, p);
    T phi_q = grid::value_at_2<D, T>(sdfGrid, q);
    T phi_r = grid::value_at_2<D, T>(sdfGrid, r);
    T phi_centroid = grid::value_at_2<D, T>(sdfGrid, (p + q + r) / T(3));

    Eigen::Matrix<T, 3, 1> c;
    if (phi_p <= phi_q && phi_p <= phi_r)
    {
        c = Eigen::Matrix<T, 3, 1>(1.0f, 0.0f, 0.0f);
    }
    else if (phi_q <= phi_p && phi_q <= phi_r)
    {
        c = Eigen::Matrix<T, 3, 1>(0.0f, 1.0f, 0.0f);
    }
    else { c = Eigen::Matrix<T, 3, 1>(0.0f, 0.0f, 1.0f); }
    if (phi_centroid <= phi_p && phi_centroid <= phi_q && phi_centroid <= phi_r)
    {
        c = (p + q + r) / T(3);
    }

    T alpha = initialAlpha;
    T prev_phi = std::numeric_limits<T>::max();

    for (int i = 0; i < maxIterations; i++)
    {
        // Compute spatial point from barycent ric coordinates
        Eigen::Matrix<T, 3, 1> x = c[0] * p + c[1] * q + c[2] * r;

        // Compute SDF value and gradient at current position
        T current_phi = grid::value_at_2(sdfGrid, x);
        Eigen::Matrix<T, 3, 1> gradPhi = computeGradient<D, T>(x, sdfGrid);

        // Check for convergence
        if (std::abs(current_phi - prev_phi) < 1e-6f) { break; }
        prev_phi = current_phi;

        // Compute gradient with respect to barycentric coordinates
        Eigen::Matrix<T, 3, 1> d;
        d[0] = gradPhi.dot(p);
        d[1] = gradPhi.dot(q);
        d[2] = gradPhi.dot(r);

        // Normalize the gradient to handle large triangles
        T grad_magnitude = d.norm();
        if (grad_magnitude > 1e-6f) { d = d / grad_magnitude; }

        // Adaptive step size based on triangle size
        T triangle_size = (p - q).norm() + (q - r).norm() + (r - p).norm();
        T adaptive_alpha = alpha / (1.0f + triangle_size * 0.1f);

        // Update barycentric coordinates
        Eigen::Matrix<T, 3, 1> c_new = c - adaptive_alpha * d;

        // Project onto simplex constraint
        project_onto_simplex(c_new);

        // Check if we're making progress
        Eigen::Matrix<T, 3, 1> x_new
            = c_new[0] * p + c_new[1] * q + c_new[2] * r;
        T new_phi = grid::value_at_2(sdfGrid, x_new);

        if (new_phi < current_phi)
        {
            c = c_new;
            // Increase step size slightly if we're making good progress
            alpha *= 1.1f;
        }
        else
        {
            // Reduce step size if not making progress
            alpha *= 0.5f;
        }

        // Clamp step size to reasonable bounds
        alpha = std::max<T>(1e-5f, std::min<T>(alpha, 1.0f));
    }

    // Calculate final results
    contactPoint = c[0] * p + c[1] * q + c[2] * r;
    penetration = grid::value_at_2(sdfGrid, contactPoint);
    normal = computeGradient<D, T>(contactPoint, sdfGrid);

    // Return true if penetration is negative (inside the object)
    return penetration < 0.05f;
}

void create_cone(double r, double h, int s, Eigen::MatrixXd& V,
                 Eigen::MatrixXi& F)
{
    V.resize(s + 2, 3);
    F.resize(2 * s, 3);

    // Tip vertex (index 0)
    V.row(0) << 0, 0, h;

    // Base vertices (indices 1 to s)
    for (int i = 0; i < s; i++)
    {
        double angle = 2 * M_PI * i / s;
        V.row(i + 1) << r * cos(angle), r * sin(angle), 0;
    }

    // Base center vertex (index s+1)
    V.row(s + 1) << 0, 0, 0;

    // Side faces (triangles connecting tip to base vertices)
    for (int i = 0; i < s; i++)
    {
        int j = (i + 1) % s;
        F.row(i) << 0, i + 1, j + 1;
    }

    // Base faces (triangles connecting base center to base vertices)
    for (int i = 0; i < s; i++)
    {
        int j = (i + 1) % s;
        F.row(s + i) << s + 1, j + 1, i + 1;
    }
}

BOOST_AUTO_TEST_CASE(grid_local_strategy)
{
    {
        using D = double;
        using T = double;

        int res = 32;

        //Read mesh
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        create_cone(1.0, 2.0, 64, V, F); // Example parameters
        std::cout << "Read mesh: #V = " << V.rows() << "  #F = " << F.rows()
                  << "\n";

        namespace fs = std::filesystem;
        fs::path source_dir = fs::path(__FILE__).parent_path();
        fs::path bunnyRelative = "../../../../../bin/resources/objs/bunny.obj";
        fs::path bunnyFull = source_dir / bunnyRelative;
        fs::path bunnyNormalized = bunnyFull.lexically_normal();

        fs::path bunnySDFRelative = "../../../../../bin/output/sdftest/sdf.csv";
        fs::path bunnySDFFull = source_dir / bunnySDFRelative;
        fs::path bunnySDFNormalized = bunnySDFFull.lexically_normal();
        fs::path conePathFull
            = source_dir / "../../../../../bin/output/sdftest/cone.obj";

        std::string conePathNormal = conePathFull.lexically_normal().string();
        // Suppose V (Eigen::MatrixXd), F (Eigen::MatrixXi) are filled by create_cone(...)
        write_obj(V, F, conePathNormal);

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

        //Copy into grid storage
        for (size_t s = 0; s < total; ++s)
        {
            G.data()[s] = static_cast<D>(S((Eigen::Index)s));
        }

        const std::string outCsv = bunnySDFNormalized.string();

        std::ofstream csv(outCsv);
        if (!csv)
        {
            std::cerr << "ERROR: cannot open " << outCsv << " for writing\n";
        }
        csv << "x,y,z,sdf\n";

        D maxZ = 0.0;
        for (size_t k = 0; k < G.K(); ++k)
        {
            for (size_t j = 0; j < G.J(); ++j)
            {
                for (size_t i = 0; i < G.I(); ++i)
                {
                    Eigen::Matrix<size_t, 3, 1> idx(i, j, k);
                    Eigen::Matrix<T, 3, 1> p;
                    grid::node_position(G, idx, p);
                    csv << std::setprecision(9) << p.x() << "," << p.y() << ","
                        << p.z() << "," << G(idx) << "\n";
                    maxZ = std::max<D>(p.z(), maxZ);
                }
            }
        }
        std::cerr << "MAXZ: " << maxZ << "\n";
        csv.close();
        std::cout << "Wrote: " << outCsv << " (" << total << " rows)\n";

        {
            BOOST_CHECK_CLOSE(-1.0, -1.0, 0.01);
        }
        //=========================================
        //NEW COLLISION LOOP
        float triangleSize = 10.0f;
        Eigen::Matrix<T, 3, 1> p(-triangleSize, -triangleSize, 2.0f);
        Eigen::Matrix<T, 3, 1> q(triangleSize, -triangleSize, 3.25f);
        Eigen::Matrix<T, 3, 1> r(0.0f, triangleSize, 4.5f);
        Eigen::Matrix<T, 3, 1> contactPoint, normal;
        T penetration;

        bool colliding = false;

        float startHeight = 2.5f;
        float endHeight = -1.0f;
        float step = -0.005f;
        for (float y = startHeight; y > endHeight; y += step)
        {
            // Update triangle position
            Eigen::Matrix<T, 3, 1> offset(0, 0, y - p.z());
            Eigen::Matrix<T, 3, 1> p_new = p + offset;
            Eigen::Matrix<T, 3, 1> q_new = q + offset;
            Eigen::Matrix<T, 3, 1> r_new = r + offset;

            // Check for collision
            colliding = optimizeTriangleFW3<D, T>(
                p_new, q_new, r_new, G, contactPoint, normal, penetration);
            std::cerr << "PENETRATION" << penetration << "\n";

            if (colliding)
            {
                std::cout << "Collision detected at z = " << y << std::endl;
                std::cout << "Contact point: " << contactPoint.transpose()
                          << std::endl;
                std::cout << "Penetration: " << penetration << std::endl;
                std::cout << "SDF value at contact: "
                          << grid::value_at(G, contactPoint) << std::endl;
                p = p_new;
                q = q_new;
                r = r_new;
                break;
            }
            else { std::cerr << "No collision for z = " << y << "\n"; }
        }

        conePathFull
            = source_dir / "../../../../../bin/output/sdftest/triangle.obj";
        conePathNormal = conePathFull.lexically_normal().string();
        // triangle: p, q, r as Eigen::Matrix<T,3,1>
        write_triangle_obj(p, q, r, conePathNormal);

        conePathFull
            = source_dir / "../../../../../bin/output/sdftest/contact.xyz";
        conePathNormal = conePathFull.lexically_normal().string();
        // contactPoint: Eigen::Matrix<T,3,1>
        write_contact_xyz(contactPoint, conePathNormal);
    }
}

BOOST_AUTO_TEST_SUITE_END();
