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
    std::cerr << "CONTACT POINT" << contactPoint << "\n";
    std::cerr << "PENETRATION: " << penetration << "\n";
    penetration = grid::value_at(cone, contactPoint);
    normal = computeGradient<D, T>(contactPoint, cone);

    // Return true if penetration is negative (inside the object)
    return penetration <= 0.0f;
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
    R = R.inverse().eval();
    Eigen::Matrix<T, 3, 1> newP = (R) * (p - translation);
    Eigen::Matrix<T, 3, 1> newQ = (R) * (q - translation);
    Eigen::Matrix<T, 3, 1> newR = (R) * (r - translation);

/*    Eigen::Matrix<T, 3, 1> newP = p;
    Eigen::Matrix<T, 3, 1> newQ = q;
    Eigen::Matrix<T, 3, 1> newR = r;*/
    T a = value_at(cone, newP);
    T b = value_at(cone, newQ);
    T c = value_at(cone, newR);

    if (a <= 0.0)
    {
        normal = computeGradient<D, T>(newP, cone);
        contactPoint = newP;
        penetration = a;
        return true;
    }

    if (b <= 0.0)
    {
        normal = computeGradient<D, T>(newQ, cone);
        contactPoint = newQ;
        penetration = b;
        return true;
    }

    if (c <= 0.0)
    {
        normal = computeGradient<D, T>(newR, cone);
        contactPoint = newR;
        penetration = c;
        return true;
    }
    return false;

/*    return optimizeTriangleFW<D, T>(newP, newQ, newR, cone, contactPoint,
                                    normal, penetration, maxIterations);*/
}
} // namespace grid

#endif // GRID_LOCAL_OPTIMIZATION_HPP
