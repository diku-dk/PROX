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
inline Eigen::Matrix<T, 3, 1> gradient_at_2(Grid<D, T> const& grid,
                                            const Eigen::Matrix<T, 3, 1>& pos)
{
    // Outside grid: mimic value_at_2 behaviour and return radial gradient
    if (pos.x() < grid.min().x() || pos.x() > grid.max().x()
        || pos.y() < grid.min().y() || pos.y() > grid.max().y()
        || pos.z() < grid.min().z() || pos.z() > grid.max().z())
    {
        Eigen::Matrix<T, 3, 1> clamped;
        clamped.x() = std::clamp(pos.x(), grid.min().x(), grid.max().x());
        clamped.y() = std::clamp(pos.y(), grid.min().y(), grid.max().y());
        clamped.z() = std::clamp(pos.z(), grid.min().z(), grid.max().z());

        D d_val = value_at_2(grid, clamped);

        Eigen::Matrix<T, 3, 1> diff = pos - clamped;
        T dist = diff.norm();

        if (dist <= std::numeric_limits<T>::epsilon())
        {
            // Exactly on boundary — fall back to inside-grid derivative if possible
            // compute gradient at clamped (which is on boundary) by proceeding below
            // (we will continue to the inside-grid code path after clamping nodes).
        }
        else
        {

            T sign = (d_val >= D(0)) ? T(1) : T(-1);
            return sign * diff / dist;
        }
    }

    // --- inside-grid derivative (or boundary point that fell through) ---

    Eigen::Matrix<size_t, 3, 1> nodes0;
    Eigen::Matrix<size_t, 3, 1> nodes1;
    enclosing_indices(grid, pos, nodes0, nodes1);

    nodes0 = Eigen::Matrix<size_t, 3, 1>(
        std::clamp<size_t>(nodes0.x(), 0, grid.I() - 1),
        std::clamp<size_t>(nodes0.y(), 0, grid.J() - 1),
        std::clamp<size_t>(nodes0.z(), 0, grid.K() - 1));
    nodes1 = Eigen::Matrix<size_t, 3, 1>(
        std::clamp<size_t>(nodes1.x(), 0, grid.I() - 1),
        std::clamp<size_t>(nodes1.y(), 0, grid.J() - 1),
        std::clamp<size_t>(nodes1.z(), 0, grid.K() - 1));

    D const d000 = grid(nodes0);
    D const d001 = grid(
        Eigen::Matrix<size_t, 3, 1>(nodes1.x(), nodes0.y(), nodes0.z())); // x+
    D const d010 = grid(
        Eigen::Matrix<size_t, 3, 1>(nodes0.x(), nodes1.y(), nodes0.z())); // y+
    D const d011 = grid(Eigen::Matrix<size_t, 3, 1>(nodes1.x(), nodes1.y(),
                                                    nodes0.z())); // x+ y+
    D const d100 = grid(
        Eigen::Matrix<size_t, 3, 1>(nodes0.x(), nodes0.y(), nodes1.z())); // z+
    D const d101 = grid(Eigen::Matrix<size_t, 3, 1>(nodes1.x(), nodes0.y(),
                                                    nodes1.z())); // x+ z+
    D const d110 = grid(Eigen::Matrix<size_t, 3, 1>(nodes0.x(), nodes1.y(),
                                                    nodes1.z())); // y+ z+
    D const d111 = grid(Eigen::Matrix<size_t, 3, 1>(nodes1.x(), nodes1.y(),
                                                    nodes1.z())); // x+ y+ z+

    // local coordinates in cell [0,1]
    const Eigen::Matrix<T, 3, 1> stu
        = (pos
           - (nodes0.template cast<T>().cwiseProduct(grid.dir()) + grid.min()))
              .cwiseQuotient(grid.dir());
    const T u = stu.x(); // x local
    const T v = stu.y(); // y local
    const T w = stu.z(); // z local

    // Map names to clearer f_{xyz} where indices are (x_bit,y_bit,z_bit):
    // f000 = d000
    // f100 = d001  (x+)
    // f010 = d010  (y+)
    // f110 = d011  (x+ y+)
    // f001 = d100  (z+)
    // f101 = d101  (x+ z+)
    // f011 = d110  (y+ z+)
    // f111 = d111  (x+ y+ z+)

    // derivative wrt local u (x)
    T df_du = (T(1) - v) * (T(1) - w) * (T(d001) - T(d000))
            + v * (T(1) - w) * (T(d011) - T(d010))
            + (T(1) - v) * w * (T(d101) - T(d100))
            + v * w * (T(d111) - T(d110));

    // derivative wrt local v (y)
    T df_dv = (T(1) - u) * (T(1) - w) * (T(d010) - T(d000))
            + u * (T(1) - w) * (T(d011) - T(d001))
            + (T(1) - u) * w * (T(d110) - T(d100))
            + u * w * (T(d111) - T(d101));

    // derivative wrt local w (z)
    T df_dw = (T(1) - u) * (T(1) - v) * (T(d100) - T(d000))
            + u * (T(1) - v) * (T(d101) - T(d001))
            + (T(1) - u) * v * (T(d110) - T(d010))
            + u * v * (T(d111) - T(d011));

    // convert to spatial derivatives by dividing by spacing (grid.dir())
    Eigen::Matrix<T, 3, 1> g;
    g.x() = df_du / grid.dir().x();
    g.y() = df_dv / grid.dir().y();
    g.z() = df_dw / grid.dir().z();

    return g;
}

// Helper to get normalized normal with safe fallback
template <typename D, typename T>
inline Eigen::Matrix<T, 3, 1> normal_at_2(Grid<D, T> const& grid,
                                          const Eigen::Matrix<T, 3, 1>& pos,
                                          T eps = T(1e-12))
{
    auto g = gradient_at_2(grid, pos);
    T nrm = g.norm();
    if (nrm <= eps)
    {
        // fallback: try small finite-diff or return a default normal

        return Eigen::Matrix<T, 3, 1>(T(0), T(0), T(1));
    }
    return g / nrm;
}

template <typename T>
inline T cubicCR(const T t, const T P0, const T P1, const T P2, const T P3)
{
    // Coeffs for Catmull-Rom (tension=0.5) as 0.5 * (a t^3 + b t^2 + c t + d)
    const T a = -P0 + 3 * P1 - 3 * P2 + P3;
    const T b = 2 * P0 - 5 * P1 + 4 * P2 - P3;
    const T c = -P0 + P2;
    const T d = 2 * P1;
    return T(0.5) * ((a * t + b) * t * t + (c * t + d));
}

// Derivative wrt t
template <typename T>
inline T cubicCR_deriv(const T t, const T P0, const T P1, const T P2,
                       const T P3)
{
    const T a = -P0 + 3 * P1 - 3 * P2 + P3;
    const T b = 2 * P0 - 5 * P1 + 4 * P2 - P3;
    const T c = -P0 + P2;
    // derivative of 0.5*(a t^3 + b t^2 + c t + d) -> 0.5*(3 a t^2 + 2 b t + c)
    return T(0.5) * ((T(3) * a * t + T(2) * b) * t + c);
}

// --- tricubic Catmull-Rom: compute value and gradient (analytic) ---
// returns pair: <value, gradient vector>
template <typename D, typename T>
inline Eigen::Matrix<T, 3, 1>
tricubicCR_value_and_gradient(Grid<D, T> const& grid,
                              const Eigen::Matrix<T, 3, 1>& pos)
{
    if (pos.x() < grid.min().x() || pos.x() > grid.max().x()
        || pos.y() < grid.min().y() || pos.y() > grid.max().y()
        || pos.z() < grid.min().z() || pos.z() > grid.max().z())
    {
        Eigen::Matrix<T, 3, 1> clamped;
        clamped.x() = std::clamp(pos.x(), grid.min().x(), grid.max().x());
        clamped.y() = std::clamp(pos.y(), grid.min().y(), grid.max().y());
        clamped.z() = std::clamp(pos.z(), grid.min().z(), grid.max().z());

        D d_val = value_at_2(grid, clamped);
        Eigen::Matrix<T, 3, 1> diff = pos - clamped;
        T dist = diff.norm();
        if (dist <= std::numeric_limits<T>::epsilon())
        {
            // fall through to interior computation at clamped point
            // set pos = clamped for subsequent computation
        }
        else
        {
            T sign = (d_val >= D(0)) ? T(1) : T(-1);
            Eigen::Matrix<T, 3, 1> g = sign * diff / dist;
            return g; // value and gradient
        }
    }

    // compute continuous grid coordinates (in grid-sample units)
    Eigen::Matrix<T, 3, 1> local = (pos - grid.min()).cwiseQuotient(grid.dir());
    // integer "central" index: floor(local)
    int ix = static_cast<int>(std::floor(local.x()));
    int iy = static_cast<int>(std::floor(local.y()));
    int iz = static_cast<int>(std::floor(local.z()));

    // local param t in [0,1]
    T tx = local.x() - T(ix);
    T ty = local.y() - T(iy);
    T tz = local.z() - T(iz);

    // indices for cubic kernel: i0 = ix-1, i1 = ix, i2 = ix+1, i3 = ix+2
    // clamp to valid grid ranges
    auto clamp_idx = [&](int v, int maxv) -> int
    {
        if (v < 0) return 0;
        if (v > maxv) return maxv;
        return v;
    };
    const int imax = static_cast<int>(grid.I()) - 1;
    const int jmax = static_cast<int>(grid.J()) - 1;
    const int kmax = static_cast<int>(grid.K()) - 1;

    int ix0 = clamp_idx(ix - 1, imax);
    int ix1 = clamp_idx(ix, imax);
    int ix2 = clamp_idx(ix + 1, imax);
    int ix3 = clamp_idx(ix + 2, imax);

    int iy0 = clamp_idx(iy - 1, jmax);
    int iy1 = clamp_idx(iy, jmax);
    int iy2 = clamp_idx(iy + 1, jmax);
    int iy3 = clamp_idx(iy + 2, jmax);

    int iz0 = clamp_idx(iz - 1, kmax);
    int iz1 = clamp_idx(iz, kmax);
    int iz2 = clamp_idx(iz + 1, kmax);
    int iz3 = clamp_idx(iz + 2, kmax);

    // Load 4x4x4 samples (64) into a small array for convenience
    // sample(xi, yj, zk) naming
    D p[4][4][4];
    int xs[4] = {ix0, ix1, ix2, ix3};
    int ys[4] = {iy0, iy1, iy2, iy3};
    int zs[4] = {iz0, iz1, iz2, iz3};
    for (int kz = 0; kz < 4; ++kz)
        for (int jy = 0; jy < 4; ++jy)
            for (int ix_ = 0; ix_ < 4; ++ix_)
                p[ix_][jy][kz] = grid(Eigen::Matrix<size_t, 3, 1>(
                    (size_t)xs[ix_], (size_t)ys[jy], (size_t)zs[kz]));

    // Stage 1: along X for each (y,z) -> compute 4x4 intermediate values: vx[jy][kz]
    T vx[4][4]; // values from cubic along x
    T vx_dx[4][4]; // derivative wrt tx of the cubic along x
    for (int kz = 0; kz < 4; ++kz)
    {
        for (int jy = 0; jy < 4; ++jy)
        {
            // samples along x: p[0..3][jy][kz]
            T P0 = T(p[0][jy][kz]);
            T P1 = T(p[1][jy][kz]);
            T P2 = T(p[2][jy][kz]);
            T P3 = T(p[3][jy][kz]);
            vx[jy][kz] = cubicCR(tx, P0, P1, P2, P3);
            vx_dx[jy][kz] = cubicCR_deriv(tx, P0, P1, P2, P3);
        }
    }

    // Stage 2: along Y for each z -> compute 4 intermediates: vy[kz] and vy_dy[kz] (dy derivative)
    T vy[4];
    T vy_dy[4];
    T vy_dx[4]; // derivatives wrt tx propagated through y (for ∂/∂x later)
    for (int kz = 0; kz < 4; ++kz)
    {
        // samples along y are vx[0..3][kz]
        vy[kz] = cubicCR(ty, vx[0][kz], vx[1][kz], vx[2][kz], vx[3][kz]);
        vy_dy[kz]
            = cubicCR_deriv(ty, vx[0][kz], vx[1][kz], vx[2][kz], vx[3][kz]);
        // propagate dx derivatives: treat vx_dx[*][kz] as samples and cubic-interp them (not deriv)
        vy_dx[kz] = cubicCR(ty, vx_dx[0][kz], vx_dx[1][kz], vx_dx[2][kz],
                            vx_dx[3][kz]);
    }

    // Stage 3: along Z final -> value and derivative in z; also propagate derivatives for x and y
    T value = cubicCR(tz, vy[0], vy[1], vy[2], vy[3]);
    T dval_dz = cubicCR_deriv(tz, vy[0], vy[1], vy[2], vy[3]);
    // derivative wrt tx (x direction): cubic in z of vy_dx
    T dval_dtx = cubicCR(tz, vy_dx[0], vy_dx[1], vy_dx[2], vy_dx[3]);
    // derivative wrt ty (y direction): we propagated vy_dy, now cubic in z
    T dval_dty = cubicCR(tz, vy_dy[0], vy_dy[1], vy_dy[2], vy_dy[3]);

    // convert param derivatives to spatial derivatives by dividing by spacing
    Eigen::Matrix<T, 3, 1> g;
    // ∂/∂x = (1/hx) * ∂/∂tx  where tx is local param along x grid samples
    g.x() = dval_dtx / grid.dir().x();
    g.y() = dval_dty / grid.dir().y();
    g.z() = dval_dz / grid.dir().z();

    return g;
}

template <typename D, typename T>
Eigen::Matrix<T, 3, 1>
computeGradient_TrilinearAnalytic(const Eigen::Matrix<T, 3, 1>& p,
                                  const grid::Grid<D, T>& grid)
{
    using Vec3T = Eigen::Matrix<T, 3, 1>;
    using Idx3 = Eigen::Matrix<size_t, 3, 1>;

    // clamp query into grid bounds (so the cell computation below is safe).
    Vec3T p_clamped;
    p_clamped.x() = std::clamp(p.x(), grid.m_min.x(), grid.m_max.x());
    p_clamped.y() = std::clamp(p.y(), grid.m_min.y(), grid.m_max.y());
    p_clamped.z() = std::clamp(p.z(), grid.m_min.z(), grid.m_max.z());

    // continuous grid coordinates in voxel units (0..N-1)
    const Vec3T diff
        = (p_clamped - grid.m_min)
              .cwiseQuotient(grid.m_max - grid.m_min)
              .cwiseProduct(Vec3T(static_cast<T>(grid.m_nodes.x() - 1),
                                  static_cast<T>(grid.m_nodes.y() - 1),
                                  static_cast<T>(grid.m_nodes.z() - 1)));
    // Alternative (if grid.dir() exists): diff = (p_clamped - grid.min()).cwiseQuotient(grid.dir());

    // compute lower cell indices (i0,j0,k0) = floor(diff) clamped to [0, N-2]
    long ix = static_cast<long>(std::floor(diff.x()));
    long iy = static_cast<long>(std::floor(diff.y()));
    long iz = static_cast<long>(std::floor(diff.z()));

    ix = std::clamp<long>(ix, 0L, static_cast<long>(grid.m_nodes.x()) - 2L);
    iy = std::clamp<long>(iy, 0L, static_cast<long>(grid.m_nodes.y()) - 2L);
    iz = std::clamp<long>(iz, 0L, static_cast<long>(grid.m_nodes.z()) - 2L);

    const size_t i0 = static_cast<size_t>(ix), i1 = i0 + 1;
    const size_t j0 = static_cast<size_t>(iy), j1 = j0 + 1;
    const size_t k0 = static_cast<size_t>(iz), k1 = k0 + 1;

    // read the 8 corner samples (reuse a single idx to avoid temporaries)
    Idx3 idx;
    D d000, d100, d010, d110, d001, d101, d011, d111;
    idx.x() = i0;
    idx.y() = j0;
    idx.z() = k0;
    d000 = grid(idx);
    idx.x() = i1;
    idx.y() = j0;
    idx.z() = k0;
    d100 = grid(idx);
    idx.x() = i0;
    idx.y() = j1;
    idx.z() = k0;
    d010 = grid(idx);
    idx.x() = i1;
    idx.y() = j1;
    idx.z() = k0;
    d110 = grid(idx);
    idx.x() = i0;
    idx.y() = j0;
    idx.z() = k1;
    d001 = grid(idx);
    idx.x() = i1;
    idx.y() = j0;
    idx.z() = k1;
    d101 = grid(idx);
    idx.x() = i0;
    idx.y() = j1;
    idx.z() = k1;
    d011 = grid(idx);
    idx.x() = i1;
    idx.y() = j1;
    idx.z() = k1;
    d111 = grid(idx);

    // compute local (u,v,w) in [0,1] relative to node (i0,j0,k0)
    // We need the world coords of node (i0,j0,k0). Use grid spacing:
    Vec3T grid_extent = (grid.m_max - grid.m_min);
    Vec3T cell_size
        = Vec3T(grid_extent.x() / static_cast<T>(grid.m_nodes.x() - 1),
                grid_extent.y() / static_cast<T>(grid.m_nodes.y() - 1),
                grid_extent.z() / static_cast<T>(grid.m_nodes.z() - 1));
    Vec3T node0_world
        = Vec3T(static_cast<T>(i0), static_cast<T>(j0), static_cast<T>(k0))
              .cwiseProduct(cell_size)
        + grid.m_min;

    Vec3T stu = (p_clamped - node0_world).cwiseQuotient(cell_size);
    // numeric safety clamp
    stu.x() = std::clamp<T>(stu.x(), (T)0, (T)1);
    stu.y() = std::clamp<T>(stu.y(), (T)0, (T)1);
    stu.z() = std::clamp<T>(stu.z(), (T)0, (T)1);

    const T u = stu.x(), v = stu.y(), w = stu.z();

    // analytic derivatives of trilinear interpolant wrt (u,v,w)
    // (these are EXACT for the trilinear interpolant)
    // ∂f/∂u (cell coords)
    T df_du = (1 - v) * (1 - w) * (static_cast<T>(d100) - static_cast<T>(d000))
            + v * (1 - w) * (static_cast<T>(d110) - static_cast<T>(d010))
            + (1 - v) * w * (static_cast<T>(d101) - static_cast<T>(d001))
            + v * w * (static_cast<T>(d111) - static_cast<T>(d011));

    // ∂f/∂v
    T df_dv = (1 - u) * (1 - w) * (static_cast<T>(d010) - static_cast<T>(d000))
            + u * (1 - w) * (static_cast<T>(d110) - static_cast<T>(d100))
            + (1 - u) * w * (static_cast<T>(d011) - static_cast<T>(d001))
            + u * w * (static_cast<T>(d111) - static_cast<T>(d101));

    // ∂f/∂w
    T df_dw = (1 - u) * (1 - v) * (static_cast<T>(d001) - static_cast<T>(d000))
            + u * (1 - v) * (static_cast<T>(d101) - static_cast<T>(d100))
            + (1 - u) * v * (static_cast<T>(d011) - static_cast<T>(d010))
            + u * v * (static_cast<T>(d111) - static_cast<T>(d110));

    // convert to spatial derivatives dividing by cell sizes
    T dfdx = df_du / cell_size.x();
    T dfdy = df_dv / cell_size.y();
    T dfdz = df_dw / cell_size.z();

    return Eigen::Matrix<T, 3, 1>(dfdx, dfdy, dfdz);
}

template <typename D, typename T>
Eigen::Matrix<T, 3, 1> computeGradient_Working(const Eigen::Matrix<T, 3, 1>& p,
                                               const grid::Grid<D, T>& grid)
{
    return computeGradient_TrilinearAnalytic(p, grid);
    //return tricubicCR_value_and_gradient(grid, p);
    // Calculate grid cell size based on resolution and bounds
    Eigen::Matrix<T, 3, 1> diff = (grid.m_max - grid.m_min);
    Eigen::Matrix<T, 3, 1> cell_size = Eigen::Matrix<T, 3, 1>(
        diff.x() / (grid.m_nodes.x() - 1), diff.y() / (grid.m_nodes.y() - 1),
        diff.z() / (grid.m_nodes.z() - 1));

    // Use cell size for finite differences
    T hx = cell_size.x() * 1.0;
    T hy = cell_size.y() * 1.0;
    T hz = cell_size.z() * 1.0;

    // Central difference for gradient approximation
    T dx = (grid::value_at_2(grid,
                             Eigen::Matrix<T, 3, 1>(p.x() + hx, p.y(), p.z()))
            - grid::value_at_2(
                grid, Eigen::Matrix<T, 3, 1>(p.x() - hx, p.y(), p.z())))
         / (2.0 * hx);

    T dy = (grid::value_at_2(grid,
                             Eigen::Matrix<T, 3, 1>(p.x(), p.y() + hy, p.z()))
            - grid::value_at_2(
                grid, Eigen::Matrix<T, 3, 1>(p.x(), p.y() - hy, p.z())))
         / (2.0 * hy);

    T dz = (grid::value_at_2(grid,
                             Eigen::Matrix<T, 3, 1>(p.x(), p.y(), p.z() + hz))
            - grid::value_at_2(
                grid, Eigen::Matrix<T, 3, 1>(p.x(), p.y(), p.z() - hz)))
         / (2.0 * hz);

    return (Eigen::Matrix<T, 3, 1>(dx, dy, dz));
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

template <typename T> struct PointPenetrations
{
    EigenVector3<T> point;
    T distToSDF;
};

template <typename T>
std::vector<PointPenetrations<T>>
getPenetrationForTris(const EigenVector3<T>& A, const EigenVector3<T>& B,
                      const EigenVector3<T>& C, const Grid<T, T>& sdf,
                      const EigenVector3<T>& sdfTrans,
                      const EigenQuaternion<T>& sdfRot, int N = 100)
{

    if (N <= 0) throw std::invalid_argument("N must be positive");

    size_t M = static_cast<size_t>(N + 1) * static_cast<size_t>(N + 2) / 2;
    std::vector<PointPenetrations<T>> pts;
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
            PointPenetrations<T> pointPenetration;
            pointPenetration.point = u * A + v * B + w * C;
            pointPenetration.distToSDF = valueAtProjection<T>(
                sdf, pointPenetration.point, sdfTrans, sdfRot);
            pts.push_back(pointPenetration);
        }
    }
    return pts;
}

template <typename D, typename T>
bool optimizeTriangleFW_NOFWA(
    const Eigen::Matrix<T, 3, 1>& p, const Eigen::Matrix<T, 3, 1>& q,
    const Eigen::Matrix<T, 3, 1>& r, const Grid<T, T>& sdf,
    const EigenVector3<T>& sdfTrans, const EigenQuaternion<T>& sdfRot,
    Eigen::Matrix<T, 3, 1>& contactPoint, Eigen::Matrix<T, 3, 1>& normal,
    T& penetration, size_t maxIterations = 32)
{
    std::vector<PointPenetrations<T>> pens
        = getPenetrationForTris(p, q, r, sdf, sdfTrans, sdfRot, 15);
    T minPenetration = std::numeric_limits<T>::max();
    EigenVector3<T> minPoint;
    for (size_t i = 0; i < pens.size(); ++i)
    {
        PointPenetrations<T> currPen = pens[i];
        if (currPen.distToSDF < minPenetration)
        {
            minPoint = currPen.point;
            minPenetration = currPen.distToSDF;
        }
    }
    contactPoint = minPoint;
    penetration = minPenetration;
    normal = gradientAtProjection(contactPoint, sdf, sdfTrans, sdfRot);
    return penetration <= T(0.0001);
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

    T max_sq = d0 > d1 ? d0 : d1;
    max_sq = d2 > max_sq ? d2 : max_sq;

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
        Eigen::Matrix<T, 3, 1> gradient
            = computeGradient_Working(x, sdf).normalized();
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
        penetration = grid::value_at_2<D, T>(sdf, x);
        if (penetration <= 0.0 + 1e-8) { break; }

        //Eigen::Matrix<T, 3, 1> sitmp = (si.transpose().eval()).dot(gradient);
        T alpha = T(2) / (T(i) + T(2));
        //xi+1 = xi...
        x = x + alpha * (si - x);
    }

    contactPoint = x;
    penetration = grid::value_at_2<D, T>(sdf, contactPoint);
    normal = computeGradient_Working(contactPoint, sdf);
    return penetration <= T(0.0001);
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
