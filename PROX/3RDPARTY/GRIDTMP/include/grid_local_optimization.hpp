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
    T hx = cell_size.x() * 5.0;
    T hy = cell_size.y() * 5.0;
    T hz = cell_size.z() * 5.0;

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
Eigen::Matrix<T, 3, 1> computeGradient_Robust(const Eigen::Matrix<T, 3, 1>& p,
                                              const grid::Grid<D, T>& grid)
{
    Eigen::Matrix<T, 3, 1> diff = (grid.m_max - grid.m_min);
    Eigen::Matrix<T, 3, 1> cell_size = Eigen::Matrix<T, 3, 1>(
        diff.x() / (grid.m_nodes.x() - 1), diff.y() / (grid.m_nodes.y() - 1),
        diff.z() / (grid.m_nodes.z() - 1));

    T hx = cell_size.x();
    T hy = cell_size.y();
    T hz = cell_size.z();

    T dx, dy, dz;

    if (p.x() - hx < grid.m_min.x())
    {
        //use forward difference
        dx = (grid::value_at_2(grid,
                               (p + Eigen::Matrix<T, 3, 1>(hx, 0, 0)).eval())
              - grid::value_at_2(grid, p))
           / hx;
    }
    else if (p.x() + hx > grid.m_max.x())
    {
        //use backward difference
        dx = (grid::value_at_2(grid, p)
              - grid::value_at_2(grid,
                                 (p - Eigen::Matrix<T, 3, 1>(hx, 0, 0)).eval()))
           / hx;
    }
    else
    {
        //Use central difference
        dx = (grid::value_at_2(grid,
                               (p + Eigen::Matrix<T, 3, 1>(hx, 0, 0)).eval())
              - grid::value_at_2(grid,
                                 (p - Eigen::Matrix<T, 3, 1>(hx, 0, 0)).eval()))
           / (2 * hx);
    }

    if (p.y() - hy < grid.m_min.y())
    {
        //use forward difference
        dy = (grid::value_at_2(grid,
                               (p + Eigen::Matrix<T, 3, 1>(hy, 0, 0)).eval())
              - grid::value_at_2(grid, p))
           / hy;
    }
    else if (p.y() + hy > grid.m_max.y())
    {
        //use backward difference
        dy = (grid::value_at_2(grid, p)
              - grid::value_at_2(grid,
                                 (p - Eigen::Matrix<T, 3, 1>(hy, 0, 0)).eval()))
           / hy;
    }
    else
    {
        //Use central difference
        dy = (grid::value_at_2(grid,
                               ((p + Eigen::Matrix<T, 3, 1>(hy, 0, 0))).eval())
              - grid::value_at_2(grid,
                                 (p - Eigen::Matrix<T, 3, 1>(hy, 0, 0)).eval()))
           / (2 * hy);
    }
    // Z component
    if (p.z() - hz < grid.m_min.z())
    {
        //use forward difference
        dz = (grid::value_at_2(grid,
                               (p + Eigen::Matrix<T, 3, 1>(hz, 0, 0)).eval())
              - grid::value_at_2(grid, p))
           / hz;
    }
    else if (p.z() + hz > grid.m_max.z())
    {
        //use backward difference
        dz = (grid::value_at_2(grid, p)
              - grid::value_at_2(grid,
                                 (p - Eigen::Matrix<T, 3, 1>(hz, 0, 0)).eval()))
           / hz;
    }
    else
    {
        //Use central difference
        dz = (grid::value_at_2(grid,
                               (p + Eigen::Matrix<T, 3, 1>(hz, 0, 0)).eval())
              - grid::value_at_2(grid,
                                 (p - Eigen::Matrix<T, 3, 1>(hz, 0, 0)).eval()))
           / (2 * hz);
    }

    return Eigen::Matrix<T, 3, 1>(dx, dy, dz);
}

template <typename T>
T computeTriangleNormalCone(const Eigen::Matrix<T, 3, 1>& p,
                            const Eigen::Matrix<T, 3, 1>& q,
                            const Eigen::Matrix<T, 3, 1>& r)
{
    // For a flat triangle, the normal cone is narrow (just the face normal)
    // Base cone angle on triangle size and curvature
    Eigen::Matrix<T, 3, 1> edges[3] = {q - p, r - q, p - r};
    T max_edge_length = T(0);
    for (int i = 0; i < 3; ++i)
    {
        max_edge_length = std::max(max_edge_length, edges[i].norm());
    }

    // Larger triangles have wider normal cones
    // You can make this more sophisticated based on your mesh curvature
    T base_cone_angle = T(5) * M_PI / T(180); // 5 degrees base
    T size_factor = std::min(
        T(1), max_edge_length / T(10)); // Scale with size up to 10 units

    return base_cone_angle * (T(1) + size_factor * T(2)); // 5-15 degrees
}

template <typename T>
Eigen::Matrix<T, 3, 1>
projectToNormalCone(const Eigen::Matrix<T, 3, 1>& input_normal,
                    const Eigen::Matrix<T, 3, 1>& cone_axis, T cone_angle)
{
    T dot_product = input_normal.dot(cone_axis);
    T angle = std::acos(std::clamp(dot_product, T(-1), T(1)));

    // If input normal is within the cone, use it as-is
    if (angle <= cone_angle) { return input_normal; }

    // Otherwise, project onto the cone boundary
    // Find the closest point on the cone boundary to the input normal

    // Project input_normal onto the plane perpendicular to cone_axis
    Eigen::Matrix<T, 3, 1> perpendicular
        = input_normal - cone_axis * dot_product;
    T perpendicular_length = perpendicular.norm();

    if (perpendicular_length < T(1e-10))
    {
        // Input normal is parallel to cone_axis but beyond the cone
        // This shouldn't happen if angle > cone_angle, but handle it
        return cone_axis;
    }

    perpendicular.normalize();

    // The projected normal is on the cone boundary
    Eigen::Matrix<T, 3, 1> projected_normal
        = cone_axis * std::cos(cone_angle)
        + perpendicular * std::sin(cone_angle);

    return projected_normal.normalized();
}

template <typename D, typename T>
void projectToIsosurface(Eigen::Matrix<T, 3, 1>& point,
                         const grid::Grid<D, T>& sdf, T target_value = T(0),
                         size_t max_projection_iters = 10, T tolerance = 1e-6)
{
    T current_value = grid::value_at_2(sdf, point);

    for (size_t i = 0; i < max_projection_iters
                       && std::abs(current_value - target_value) > tolerance;
         ++i)
    {
        // Compute gradient at current point
        Eigen::Matrix<T, 3, 1> grad = computeGradient_Working(point, sdf);

        // Avoid division by zero
        T grad_norm_sq = grad.squaredNorm();
        if (grad_norm_sq < 1e-12) { break; }

        // Move in the direction of the gradient to reach the isosurface
        // The step is proportional to the current SDF value
        T step = (current_value - target_value) / grad_norm_sq;
        point = point - step * grad;

        current_value = grid::value_at_2(sdf, point);
    }
}

template <typename T>
Eigen::Matrix<T, 3, 1> projectToNormalConeIntersection(
    const Eigen::Matrix<T, 3, 1>& initial_normal,
    const Eigen::Matrix<T, 3, 1>& triangle_normal,
    const Eigen::Matrix<T, 3, 1>& sdf_gradient,
    T triangle_cone_angle = T(M_PI / 6), // 30 degrees default
    T sdf_cone_angle = T(M_PI / 4)) // 45 degrees default
{
    // Normalize inputs
    Eigen::Matrix<T, 3, 1> n_tri = triangle_normal.normalized();
    Eigen::Matrix<T, 3, 1> n_sdf = sdf_gradient.normalized();
    Eigen::Matrix<T, 3, 1> n_initial = initial_normal.normalized();

    // Check if initial normal is already within both cones
    T dot_tri = n_initial.dot(n_tri);
    T dot_sdf = n_initial.dot(n_sdf);

    T cos_tri_angle = std::cos(triangle_cone_angle);
    T cos_sdf_angle = std::cos(sdf_cone_angle);

    if (dot_tri >= cos_tri_angle && dot_sdf >= cos_sdf_angle)
    {
        return n_initial; // Already in intersection
    }

    // Compute the plane of intersection between the two cones
    Eigen::Matrix<T, 3, 1> axis1 = n_tri;
    Eigen::Matrix<T, 3, 1> axis2
        = (n_sdf - n_sdf.dot(n_tri) * n_tri).normalized();

    // Parameterize the search in 2D
    auto compute_normal = [&](T angle) -> Eigen::Matrix<T, 3, 1> {
        return (std::cos(angle) * axis1 + std::sin(angle) * axis2).normalized();
    };

    // Find the angle range that satisfies both cone constraints
    T best_angle = 0;
    T min_distance = std::numeric_limits<T>::max();

    // Sample multiple angles to find the best projection
    const int samples = 32;
    for (int i = 0; i <= samples; ++i)
    {
        T angle = (T(i) / T(samples)) * T(2 * M_PI);
        Eigen::Matrix<T, 3, 1> candidate = compute_normal(angle);

        T dot_with_tri = candidate.dot(n_tri);
        T dot_with_sdf = candidate.dot(n_sdf);

        // Check if within both cones
        bool in_tri_cone = (dot_with_tri >= cos_tri_angle);
        bool in_sdf_cone = (dot_with_sdf >= cos_sdf_angle);

        if (in_tri_cone && in_sdf_cone)
        {
            // This candidate is in the intersection, use it directly
            return candidate;
        }

        // Compute distance to initial normal for candidates outside intersection
        T distance = (candidate - n_initial).norm();
        if (distance < min_distance)
        {
            min_distance = distance;
            best_angle = angle;
        }
    }

    // If no point found in intersection, find closest point to intersection boundary
    // This is more sophisticated - we want the point that minimizes distance to both cones
    Eigen::Matrix<T, 3, 1> projected = compute_normal(best_angle);

    // If still not in intersection, project to the nearest cone boundary
    if (projected.dot(n_tri) < cos_tri_angle)
    {
        // Project to triangle cone boundary
        Eigen::Matrix<T, 3, 1> rejection
            = projected - projected.dot(n_tri) * n_tri;
        if (rejection.norm() > 1e-12) { rejection.normalize(); }
        projected = cos_tri_angle * n_tri
                  + std::sin(std::acos(cos_tri_angle)) * rejection;
    }

    if (projected.dot(n_sdf) < cos_sdf_angle)
    {
        // Project to SDF cone boundary
        Eigen::Matrix<T, 3, 1> rejection
            = projected - projected.dot(n_sdf) * n_sdf;
        if (rejection.norm() > 1e-12) { rejection.normalize(); }
        projected = cos_sdf_angle * n_sdf
                  + std::sin(std::acos(cos_sdf_angle)) * rejection;
    }

    return projected.normalized();
}

template <typename T>
T estimateTriangleNormalConeAngle(const Eigen::Matrix<T, 3, 1>& p,
                                  const Eigen::Matrix<T, 3, 1>& q,
                                  const Eigen::Matrix<T, 3, 1>& r)
{
    // Compute triangle normal
    Eigen::Matrix<T, 3, 1> tri_normal = (q - p).cross(r - p).normalized();

    // Estimate cone angle based on triangle shape
    // For a perfectly flat triangle, angle is small
    // For highly curved regions, angle is larger

    // Compute edge lengths
    T a = (q - p).norm();
    T b = (r - q).norm();
    T c = (p - r).norm();

    // Use triangle aspect ratio to estimate curvature
    T max_edge = std::max({a, b, c});
    T min_edge = std::min({a, b, c});
    T aspect_ratio = max_edge / (min_edge + 1e-12);

    // Base angle + aspect ratio contribution
    T base_angle = T(M_PI / 12); // 15 degrees base
    T aspect_contribution = std::min(T(0.2), T(0.05) * (aspect_ratio - T(1)));

    return base_angle + aspect_contribution;
}

template <typename D, typename T>
T estimateSDFNormalConeAngle(const Eigen::Matrix<T, 3, 1>& point,
                             const grid::Grid<D, T>& sdf,
                             T search_radius = T(0.01))
{
    // Sample SDF normals in neighborhood to estimate local curvature
    std::vector<Eigen::Matrix<T, 3, 1>> neighbor_normals;

    // Sample in a small neighborhood
    const int samples = 8;
    for (int i = 0; i < samples; ++i)
    {
        T angle = (T(i) / T(samples)) * T(2 * M_PI);
        Eigen::Matrix<T, 3, 1> offset(search_radius * std::cos(angle),
                                      search_radius * std::sin(angle), T(0));

        // Apply random rotation to avoid bias
        Eigen::Matrix<T, 3, 1> sample_point = point + offset;
        Eigen::Matrix<T, 3, 1> grad
            = computeGradient_Working(sample_point, sdf);

        if (grad.norm() > 1e-12)
        {
            neighbor_normals.push_back(grad.normalized());
        }
    }

    if (neighbor_normals.empty())
    {
        return T(M_PI / 6); // Default 30 degrees
    }

    // Compute maximum angle between center normal and neighbors
    Eigen::Matrix<T, 3, 1> center_normal
        = computeGradient_Working(point, sdf).normalized();
    T max_angle = 0;

    for (const auto& neighbor : neighbor_normals)
    {
        T dot_product = center_normal.dot(neighbor);
        T angle = std::acos(std::clamp(dot_product, T(-1), T(1)));
        max_angle = std::max(max_angle, angle);
    }

    // Add safety margin
    return max_angle * T(1.5) + T(M_PI / 36); // + 5 degrees
}

template <typename D, typename T>
bool optimizeTriangleFW_Working_old(const Eigen::Matrix<T, 3, 1>& p,
                                    const Eigen::Matrix<T, 3, 1>& q,
                                    const Eigen::Matrix<T, 3, 1>& r,
                                    const grid::Grid<D, T>& sdf,
                                    Eigen::Matrix<T, 3, 1>& contactPoint,
                                    Eigen::Matrix<T, 3, 1>& normal,
                                    T& penetration, size_t maxIterations = 32)
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
    penetration = grid::value_at_2<D, T>(sdf, contactPoint);
    if (penetration <= T(0))
    {
        projectToIsosurface(contactPoint, sdf, T(0));

        // Compute initial SDF normal
        Eigen::Matrix<T, 3, 1> sdf_normal
            = computeGradient_Working(contactPoint, sdf).normalized();

        // Compute triangle normal
        Eigen::Matrix<T, 3, 1> triangle_normal
            = (q - p).cross(r - p).normalized();

        // Estimate cone angles
        T tri_cone_angle = estimateTriangleNormalConeAngle(p, q, r);
        T sdf_cone_angle = estimateSDFNormalConeAngle(contactPoint, sdf);

        // Project to cone intersection
        normal = projectToNormalConeIntersection(sdf_normal, triangle_normal,
                                                 sdf_normal, tri_cone_angle,
                                                 sdf_cone_angle);

        contactPoint = x;
    }

    return penetration <= T(0);
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
    penetration = grid::value_at_2<D, T>(sdf, contactPoint);
    normal = computeGradient_Working(contactPoint, sdf);
    return penetration <= T(0);
}

template <typename D, typename T>
bool optimizeTriangleFWTransform(
    const Eigen::Matrix<T, 3, 1>& p, const Eigen::Matrix<T, 3, 1>& q,
    const Eigen::Matrix<T, 3, 1>& r, const Eigen::Matrix<T, 3, 1>& translation,
    const Eigen::Quaternion<T> rotation, const grid::Grid<D, T>& cone,
    Eigen::Matrix<T, 3, 1>& contactPoint, Eigen::Matrix<T, 3, 1>& normal,
    T& penetration, int maxIterations = 32)
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
