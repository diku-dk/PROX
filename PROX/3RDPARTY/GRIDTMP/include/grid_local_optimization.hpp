#ifndef GRID_LOCAL_OPTIMIZATION_HPP
#define GRID_LOCAL_OPTIMIZATION_HPP

#include <grid_grid.h>
#include <grid_value_at.h>
#include <eigenhelperall.h>

namespace grid
{
template <typename D, typename T>
Eigen::Matrix<T, 3, 1> computeGradient(const Eigen::Matrix<T, 3, 1>& p,
                                       const grid::Grid<D, T>& grid)
{
    float h = 0.075;
    float dx
        = (grid::value_at(grid, Eigen::Matrix<T, 3, 1>(p.x() + h, p.y(), p.z()))
           - grid::value_at(grid,
                            Eigen::Matrix<T, 3, 1>(p.x() - h, p.y(), p.z())))
        / (2 * h);
    float dy
        = (grid::value_at(grid, Eigen::Matrix<T, 3, 1>(p.x(), p.y() + h, p.z()))
           - grid::value_at(grid,
                            Eigen::Matrix<T, 3, 1>(p.x(), p.y() - h, p.z())))
        / (2 * h);
    float dz
        = (grid::value_at(grid, Eigen::Matrix<T, 3, 1>(p.x(), p.y(), p.z() + h))
           - grid::value_at(grid,
                            Eigen::Matrix<T, 3, 1>(p.x(), p.y(), p.z() - h)))
        / (2 * h);
    return Eigen::Matrix<T, 3, 1>(dx, dy, dz).normalized();
}

template <typename D, typename T>
bool optimizeTriangleFW(const Eigen::Matrix<T, 3, 1>& p,
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
        float alpha = 2.0f / (i + 2.0f);

        // Update position
        x = x + alpha * (s - x);
    }

    // Calculate final results
    contactPoint = x;
    penetration = grid::value_at(cone, contactPoint);
    std::cerr << "CONTACT POINT" << contactPoint << "\n";
    std::cerr << "PENETRATION: " << penetration << "\n";
    normal = computeGradient<D, T>(contactPoint, cone);

    // Return true if penetration is negative (inside the object)
    return penetration <= 0.0f;
}

template <typename D, typename T>
Eigen::Matrix<T, 3, 1> computeGradient_Working(const Eigen::Matrix<T, 3, 1>& p,
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

template <typename D, typename T>
bool optimizeTriangleFW_Working(const Eigen::Matrix<T, 3, 1>& p,
                                const Eigen::Matrix<T, 3, 1>& q,
                                const Eigen::Matrix<T, 3, 1>& r,
                                const grid::Grid<D, T>& sdf,
                                Eigen::Matrix<T, 3, 1>& contactPoint,
                                Eigen::Matrix<T, 3, 1>& normal, T& penetration,
                                size_t maxIterations = 32)
{
    //Early out: If sphere-radius r of tri is less than the sampled centroid x_c in our SDF, i.e. \phi(x_c) < r, ignore!
    Eigen::Matrix<T, 3, 1> centroid = (p + q + r) * (T(1) / T(3));

    T d0 = (p - centroid).squaredNorm();
    T d1 = (q - centroid).squaredNorm();
    T d2 = (r - centroid).squaredNorm();

    //pairwise max without initializer-list temporaries
    T max_sq = d0 > d1 ? d0 : d1;
    max_sq = d2 > max_sq ? d2 : max_sq;

    // Avoid sqrt: radius = sqrt(max_sq).
    // Condition phi_centroid >= radius  <=>  phi_centroid >= 0 && phi_centroid*phi_centroid >= max_sq
    T phi_centroid = grid::value_at_2<D, T>(sdf, centroid);
    if (phi_centroid >= T(0) && (phi_centroid * phi_centroid) >= max_sq)
    {
        //return false;
    }

    //Better initialization: evaluate at vertices and choose the one with smallest SDF value
    T phi_p = grid::value_at_2<D, T>(sdf, p);
    T phi_q = grid::value_at_2<D, T>(sdf, q);
    T phi_r = grid::value_at_2<D, T>(sdf, r);

    Eigen::Matrix<T, 3, 1> x;
    if (phi_p <= phi_q && phi_p <= phi_r) { x = p; }
    else if (phi_q <= phi_p && phi_q <= phi_r) { x = q; }
    else { x = r; }

    T threshold = 1e-12;
    for (size_t i = 0; i < maxIterations; ++i)
    {
        Eigen::Matrix<T, 3, 1> gradient = computeGradient_Working(x, sdf);
        Eigen::Matrix<T, 3, 1> gradientTransposed = gradient.transpose();
        T Lp = gradientTransposed.dot(p);
        T Lq = gradientTransposed.dot(q);
        T Lr = gradientTransposed.dot(r);
        Eigen::Matrix<T, 3, 1> si;
        if (Lp <= Lq && Lp <= Lr) { si = p; }
        else if (Lq <= Lp && Lq <= Lr) { si = q; }
        else { si = r; }

        /*        T subproblem = (si.transpose()).dot(- gradient);
        if (subproblem < threshold)
        {
            //Converged
            break;
        }*/

        //Eigen::Matrix<T, 3, 1> sitmp = (si.transpose().eval()).dot(gradient);
        T alpha = T(2) / (T(i) + T(2));
        //xi+1 = xi...
        x = x + alpha * (si - x);
    }

    contactPoint = x;
    normal = computeGradient_Working(x, sdf).normalized();
    penetration = grid::value_at_2<D, T>(sdf, x);
    return penetration <= T(0);
}

template <typename D, typename T>
bool optimizeTriangleFWTransform(
    const Eigen::Matrix<T, 3, 1>& p, const Eigen::Matrix<T, 3, 1>& q,
    const Eigen::Matrix<T, 3, 1>& r, const Eigen::Matrix<T, 3, 1>& translation,
    const Eigen::Quaternion<T> rotation, const grid::Grid<D, T>& cone,
    Eigen::Matrix<T, 3, 1>& contactPoint, Eigen::Matrix<T, 3, 1>& normal,
    T& penetration, int maxIterations = 20)
{
    /*Eigen::Quaternion<T> conj = rotation.conjugate();
    Eigen::Matrix<T, 3, 1> newP = (conj) * (p - translation);
    Eigen::Matrix<T, 3, 1> newQ = (conj) * (q - translation);
    Eigen::Matrix<T, 3, 1> newR = (conj) * (r - translation);*/

    //CURR
    Eigen::Matrix3<T> R = rotation.toRotationMatrix();
    Eigen::Matrix3<T> RTrans = R.transpose().eval();
    Eigen::Matrix<T, 3, 1> newP = (RTrans) * (p - translation);
    Eigen::Matrix<T, 3, 1> newQ = (RTrans) * (q - translation);
    Eigen::Matrix<T, 3, 1> newR = (RTrans) * (r - translation);

    /*    Eigen::Matrix<T, 3, 1> newP = p;
    Eigen::Matrix<T, 3, 1> newQ = q;
    Eigen::Matrix<T, 3, 1> newR = r;*/
    T a = value_at(cone, newP);
    T b = value_at(cone, newQ);
    T c = value_at(cone, newR);

    std::cerr << "MIN BOUNDS " << cone.min() << "\n";
    std::cerr << "MAX BOUNDS " << cone.max() << "\n";

    if (a <= 0.0)
    {
        normal = computeGradient<D, T>(newP, cone);
        contactPoint = R * newP + translation;
        penetration = a;
        return true;
    }

    if (b <= 0.0)
    {
        std::cerr << "DATA" << "\n";
        std::cerr << cone.data();
        normal = computeGradient<D, T>(newQ, cone);
        contactPoint = R * newQ + translation;
        penetration = b;
        return true;
    }

    if (c <= 0.0)
    {
        normal = computeGradient<D, T>(newR, cone);
        contactPoint = R * newR + translation;
        penetration = c;
        return true;
    }
    return false;

    /*    return optimizeTriangleFW<D, T>(newP, newQ, newR, cone, contactPoint,
                                    normal, penetration, maxIterations);*/
}
} // namespace grid

#endif // GRID_LOCAL_OPTIMIZATION_HPP
