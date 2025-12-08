#ifndef GRID_VALUE_AT_H
#define GRID_VALUE_AT_H

#include <grid_enclosing_indices.h>

#include <boost/numeric/conversion/bounds.hpp>
#include <boost/cast.hpp>

#include <cmath>

namespace grid
{
    template<typename D,typename T>
    inline D trilinear(
        const D d000
        , const D d001
        , const D d010
        , const D d011
        , const D d100
        , const D d101
        , const D d110
        , const D d111
        , const Eigen::Matrix<T, 3, 1>& stu
        )
    {
        D const x00 = static_cast<D>( (d001 - d000)*stu.x() + d000 );
        D const x01 = static_cast<D>( (d011 - d010)*stu.x() + d010 );
        D const x10 = static_cast<D>( (d101 - d100)*stu.x() + d100 );
        D const x11 = static_cast<D>( (d111 - d110)*stu.x() + d110 );
        D const y0  = static_cast<D>( ( x01 -  x00)*stu.y() +  x00 );
        D const y1  = static_cast<D>( ( x11 -  x10)*stu.y() +  x10 );
        return static_cast<D>( (y1-y0)*stu.z() + y0 );
    }

    /*template<typename D,typename T>
    inline typename std::enable_if<!std::is_same<D,T>::value, D>::type value_at(Grid<D,T> const & grid, const Eigen::Matrix<T, 3, 1>& pos)
    {
        (void)0;
        return value_at<D,T>(grid,pos);
    }*/

    template<typename D,typename T>
    inline D value_at(Grid<D,T> const & grid, const Eigen::Matrix<T, 3, 1>& pos)
    {
        Eigen::Matrix<size_t, 3, 1> nodes0;
        Eigen::Matrix<size_t, 3, 1> nodes1;
        enclosing_indices( grid, pos, nodes0, nodes1);

        nodes0 = Eigen::Matrix<size_t, 3, 1>(
            std::max<size_t>(std::min<T>(nodes0.x(), grid.I() - 1), size_t(0)),
            std::max<size_t>(std::min<T>(nodes0.y(), grid.J() - 1), size_t(0)),
            std::max<size_t>(std::min<T>(nodes0.z(), grid.K() - 1), size_t(0)));
        nodes1 = Eigen::Matrix<size_t, 3, 1>(
            std::max<size_t>(std::min<size_t>(nodes1.x(), grid.I() - 1),
                             size_t(0)),
            std::max<size_t>(std::min<size_t>(nodes1.y(), grid.J() - 1),
                             size_t(0)),
            std::max<size_t>(std::min<size_t>(nodes1.z(), grid.K() - 1),
                             size_t(0)));
        D const d000 = grid(nodes0);
        D const d001 = grid( Eigen::Matrix<size_t, 3, 1>(nodes1.x(), nodes0.y(), nodes0.z()) );
        D const d010 = grid( Eigen::Matrix<size_t, 3, 1>(nodes0.x(), nodes1.y(), nodes0.z()) );
        D const d011 = grid( Eigen::Matrix<size_t, 3, 1>(nodes1.x(), nodes1.y(), nodes0.z()) );
        D const d100 = grid( Eigen::Matrix<size_t, 3, 1>(nodes0.x(), nodes0.y(), nodes1.z()) );
        D const d101 = grid( Eigen::Matrix<size_t, 3, 1>(nodes1.x(), nodes0.y(), nodes1.z()) );
        D const d110 = grid( Eigen::Matrix<size_t, 3, 1>(nodes0.x(), nodes1.y(), nodes1.z()) );
        D const d111 = grid( Eigen::Matrix<size_t, 3, 1>(nodes1.x(), nodes1.y(), nodes1.z()) );

        const Eigen::Matrix<T, 3, 1> stu = (pos - (nodes0.template cast<T>().cwiseProduct(grid.dir()) + grid.min())) .cwiseQuotient(grid.dir());

        return trilinear<D,T>( d000, d001, d010, d011, d100, d101, d110, d111, stu );
    }

    template <typename D, typename T>
    inline void enclosing_cell_indices2(
        Grid<D, T> const& grid, const Eigen::Matrix<T, 3, 1>& pos,
        Eigen::Matrix<size_t, 3, 1>& node0, // lower corner indices
        Eigen::Matrix<size_t, 3, 1>& node1 // upper corner indices = node0 + 1
    )
    {
        using std::floor;

        // compute continuous grid coordinates (in grid spacing units)
        Eigen::Matrix<T, 3, 1> diff
            = (pos - grid.min()).cwiseQuotient(grid.dir());

        // floor to get lower vertex indices (may be negative)
        long ix = static_cast<long>(std::floor(diff.x()));
        long iy = static_cast<long>(std::floor(diff.y()));
        long iz = static_cast<long>(std::floor(diff.z()));

        // Ensure grid has at least 2 samples along each axis
        assert(grid.I() >= 2 && grid.J() >= 2 && grid.K() >= 2);

        // clamp to valid lower-corner index range [0, size-2] so node1 = node0 + 1 is valid
        long max_ix = static_cast<long>(grid.I()) - 2;
        long max_iy = static_cast<long>(grid.J()) - 2;
        long max_iz = static_cast<long>(grid.K()) - 2;

        ix = std::clamp(ix, 0l, max_ix);
        iy = std::clamp(iy, 0l, max_iy);
        iz = std::clamp(iz, 0l, max_iz);

        node0.x() = static_cast<size_t>(ix);
        node0.y() = static_cast<size_t>(iy);
        node0.z() = static_cast<size_t>(iz);

        node1 = node0 + Eigen::Matrix<size_t, 3, 1>(1, 1, 1);
    }

    template <typename D, typename T> inline D lerp(const D& a, const D& b, T t)
    {
        return static_cast<D>(a + (b - a) * t);
    }

    // Trilinear with explicit corner ordering: (i,j,k) where i=x, j=y, k=z.
    // Arguments ordered so it's straightforward to match grid indices.
    template <typename D, typename T>
    inline D
    trilinear_sample(const D d000, // (i0, j0, k0)
                     const D d100, // (i1, j0, k0)
                     const D d010, // (i0, j1, k0)
                     const D d110, // (i1, j1, k0)
                     const D d001, // (i0, j0, k1)
                     const D d101, // (i1, j0, k1)
                     const D d011, // (i0, j1, k1)
                     const D d111, // (i1, j1, k1)
                     const Eigen::Matrix<T, 3, 1>& stu // (u,v,w) in [0,1]
    )
    {
        // use double-like type for intermediate weights for numeric stability
        T u = stu.x();
        T v = stu.y();
        T w = stu.z();

        // Interpolate along x for k0 / k1 slices
        D c00 = lerp(d000, d100, u); // z=k0, y=j0..j0
        D c10 = lerp(d010, d110, u); // z=k0, y=j1..j1
        D c01 = lerp(d001, d101, u); // z=k1, y=j0..j0
        D c11 = lerp(d011, d111, u); // z=k1, y=j1..j1

        // Interpolate along y
        D c0 = lerp(c00, c10, v); // z=k0
        D c1 = lerp(c01, c11, v); // z=k1

        // Interpolate along z
        return lerp(c0, c1, w);
    }

    // Internal fast trilinear sample (no bounds recursion). Uses new enclosing_cell_indices2.
    template <typename D, typename T>
    inline D sample_trilinear_local(Grid<D, T> const& grid,
                                    const Eigen::Matrix<T, 3, 1>& pos)
    {
        Eigen::Matrix<size_t, 3, 1> n0, n1;
        enclosing_cell_indices2(grid, pos, n0, n1);

        const size_t i0 = n0.x(), i1 = n1.x();
        const size_t j0 = n0.y(), j1 = n1.y();
        const size_t k0 = n0.z(), k1 = n1.z();

        // read corners with explicit mapping (i=x, j=y, k=z)
        D d000 = grid(Eigen::Matrix<size_t, 3, 1>(i0, j0, k0));
        D d100 = grid(Eigen::Matrix<size_t, 3, 1>(i1, j0, k0));
        D d010 = grid(Eigen::Matrix<size_t, 3, 1>(i0, j1, k0));
        D d110 = grid(Eigen::Matrix<size_t, 3, 1>(i1, j1, k0));
        D d001 = grid(Eigen::Matrix<size_t, 3, 1>(i0, j0, k1));
        D d101 = grid(Eigen::Matrix<size_t, 3, 1>(i1, j0, k1));
        D d011 = grid(Eigen::Matrix<size_t, 3, 1>(i0, j1, k1));
        D d111 = grid(Eigen::Matrix<size_t, 3, 1>(i1, j1, k1));

        // compute local (u,v,w) in [0,1] using the lower node position in world coords
        Eigen::Matrix<T, 3, 1> node0_world
            = n0.template cast<T>().cwiseProduct(grid.dir()) + grid.min();
        Eigen::Matrix<T, 3, 1> stu
            = (pos - node0_world).cwiseQuotient(grid.dir());

        // clamp numeric round-off (in case pos is exactly on max) to [0,1]
        stu.x() = std::clamp<T>(stu.x(), (T)0, (T)1);
        stu.y() = std::clamp<T>(stu.y(), (T)0, (T)1);
        stu.z() = std::clamp<T>(stu.z(), (T)0, (T)1);

        return trilinear_sample<D, T>(d000, d100, d010, d110, d001, d101, d011,
                                      d111, stu);
    }

    template <typename T>
    inline void quadratic_basis(const T u, T& Lm1, T& L0, T& L1)
    {
        Lm1 = (u * (u - (T)1)) * (T)0.5; // u*(u-1)/2
        L0 = (T)1 - u * u; // 1 - u^2
        L1 = (u * (u + (T)1)) * (T)0.5; // u*(u+1)/2
    }

    // local triquadratic sampler: uses 3x3x3 samples around nearest node (m-1,m,m+1)
    template <typename D, typename T>
    inline D sample_triquadratic_local(Grid<D, T> const& grid,
                                       const Eigen::Matrix<T, 3, 1>& pos)
    {
        using Vec3T = Eigen::Matrix<T, 3, 1>;
        using Idx3 = Eigen::Matrix<size_t, 3, 1>;

        // compute continuous grid coords
        Vec3T diff = (pos - grid.min()).cwiseQuotient(grid.dir());

        // if grid is too small to form a 3-wide stencil, fallback
        if (grid.I() < 2 || grid.J() < 2 || grid.K() < 2)
        {
            // fallback to nearest-vertex sample
            Idx3 idx;
            idx.x() = static_cast<size_t>(std::clamp<long>(
                std::lround(diff.x()), 0L,
                static_cast<long>(std::max<size_t>(1, grid.I()) - 1)));
            idx.y() = static_cast<size_t>(std::clamp<long>(
                std::lround(diff.y()), 0L,
                static_cast<long>(std::max<size_t>(1, grid.J()) - 1)));
            idx.z() = static_cast<size_t>(std::clamp<long>(
                std::lround(diff.z()), 0L,
                static_cast<long>(std::max<size_t>(1, grid.K()) - 1)));
            return grid(idx);
        }

        // choose central vertex index m = round(diff)
        long mx = static_cast<long>(std::lround(diff.x()));
        long my = static_cast<long>(std::lround(diff.y()));
        long mz = static_cast<long>(std::lround(diff.z()));

        // clamp m so that m-1 .. m+1 remain within [0, N-1]
        long min_mx = 1, max_mx = static_cast<long>(grid.I()) - 2;
        long min_my = 1, max_my = static_cast<long>(grid.J()) - 2;
        long min_mz = 1, max_mz = static_cast<long>(grid.K()) - 2;

        // if the grid dimension is exactly 2, adjust clamps (will fallback below to trilinear)
        if (grid.I() < 3) { min_mx = max_mx = 0; }
        if (grid.J() < 3) { min_my = max_my = 0; }
        if (grid.K() < 3) { min_mz = max_mz = 0; }

        mx = std::clamp(mx, min_mx, max_mx);
        my = std::clamp(my, min_my, max_my);
        mz = std::clamp(mz, min_mz, max_mz);

        // compute local u = (x - x_m)/dx where x_m is world coord of vertex m
        Vec3T m_world
            = Eigen::Matrix<T, 3, 1>(static_cast<T>(mx), static_cast<T>(my),
                                     static_cast<T>(mz))
                  .cwiseProduct(grid.dir())
            + grid.min();
        T u = (pos.x() - m_world.x()) / grid.dir().x();
        T v = (pos.y() - m_world.y()) / grid.dir().y();
        T w = (pos.z() - m_world.z()) / grid.dir().z();

        // basis functions (degree-2 Lagrange on nodes -1,0,1)
        T Lxu_m1, Lxu_0, Lxu_1;
        T Lyv_m1, Lyv_0, Lyv_1;
        T Lzw_m1, Lzw_0, Lzw_1;
        quadratic_basis(u, Lxu_m1, Lxu_0, Lxu_1);
        quadratic_basis(v, Lyv_m1, Lyv_0, Lyv_1);
        quadratic_basis(w, Lzw_m1, Lzw_0, Lzw_1);

        // sample 3x3x3 values centered at (mx,my,mz) using offsets -1,0,1
        // accumulate separably: first along x for each (jy,kz), then y, then z.
        D tmp_yz[3][3]; // after x interpolation -> tmp_yz[j+1][k+1]
        for (int dj = -1; dj <= 1; ++dj)
        {
            for (int dk = -1; dk <= 1; ++dk)
            {
                // x-interpolation for fixed (j,k)
                long jidx = static_cast<long>(my + dj);
                long kidx = static_cast<long>(mz + dk);

                // clamp vertex indices to valid range [0, N-1]
                size_t iA = static_cast<size_t>(std::clamp<long>(
                    mx - 1, 0L, static_cast<long>(grid.I() - 1)));
                size_t iB = static_cast<size_t>(std::clamp<long>(
                    mx + 0, 0L, static_cast<long>(grid.I() - 1)));
                size_t iC = static_cast<size_t>(std::clamp<long>(
                    mx + 1, 0L, static_cast<long>(grid.I() - 1)));

                size_t jS = static_cast<size_t>(std::clamp<long>(
                    jidx, 0L, static_cast<long>(grid.J() - 1)));
                size_t kS = static_cast<size_t>(std::clamp<long>(
                    kidx, 0L, static_cast<long>(grid.K() - 1)));

                D fA = grid(Eigen::Matrix<size_t, 3, 1>(iA, jS, kS));
                D fB = grid(Eigen::Matrix<size_t, 3, 1>(iB, jS, kS));
                D fC = grid(Eigen::Matrix<size_t, 3, 1>(iC, jS, kS));

                D interp_x
                    = static_cast<D>(Lxu_m1 * fA + Lxu_0 * fB + Lxu_1 * fC);
                tmp_yz[dj + 1][dk + 1] = interp_x;
            }
        }

        // now interpolate along y for each k
        D tmp_z[3];
        for (int dk = 0; dk < 3; ++dk)
        {
            D c0 = tmp_yz[0][dk];
            D c1 = tmp_yz[1][dk];
            D c2 = tmp_yz[2][dk];
            tmp_z[dk] = static_cast<D>(Lyv_m1 * c0 + Lyv_0 * c1 + Lyv_1 * c2);
        }

        // finally along z
        D cfinal = static_cast<D>(Lzw_m1 * tmp_z[0] + Lzw_0 * tmp_z[1]
                                  + Lzw_1 * tmp_z[2]);
        return cfinal;
    }

    // fallback trilinear sampler (compact, similar to your previous trilinear)
    template <typename D, typename T>
    inline D sample_trilinear_fallback(Grid<D, T> const& grid,
                                       const Eigen::Matrix<T, 3, 1>& pos)
    {
        // (re-use your earlier code's logic here or call your existing value_at)
        // For brevity, assume you have `value_at(grid,pos)` implemented as earlier.
        return value_at<D, T>(grid,
                              pos); // <-- your existing trilinear function
    }

    // The improved value_at_2: outside handling + fast 27-point quadratic sampler
    template <typename D, typename T>
    inline D value_at_2_quadratic(Grid<D, T> const& grid,
                                  const Eigen::Matrix<T, 3, 1>& pos)
    {
        // preserve your outside-grid behavior exactly, but call local sampler instead of recursion
        if (pos.x() < grid.min().x() || pos.x() > grid.max().x()
            || pos.y() < grid.min().y() || pos.y() > grid.max().y()
            || pos.z() < grid.min().z() || pos.z() > grid.max().z())
        {
            // clamp to nearest point on grid boundary
            Eigen::Matrix<T, 3, 1> clamped;
            clamped.x() = std::clamp(pos.x(), grid.min().x(), grid.max().x());
            clamped.y() = std::clamp(pos.y(), grid.min().y(), grid.max().y());
            clamped.z() = std::clamp(pos.z(), grid.min().z(), grid.max().z());

            // sample at clamped (use local sampler; avoids recursion)
            D d_val;
            // if grid too small for quadratic, fallback to trilinear sample
            if (grid.I() >= 3 && grid.J() >= 3 && grid.K() >= 3)
                d_val = sample_triquadratic_local<D, T>(grid, clamped);
            else
                d_val = sample_trilinear_fallback<D, T>(grid, clamped);

            T euclidean_dist = (pos - clamped).norm();
            return d_val >= static_cast<D>(0)
                     ? d_val + static_cast<D>(euclidean_dist)
                     : d_val - static_cast<D>(euclidean_dist);
        }

        // inside grid: use triquadratic if possible, otherwise fallback
        if (grid.I() >= 3 && grid.J() >= 3 && grid.K() >= 3)
            return sample_triquadratic_local<D, T>(grid, pos);
        else
            return sample_trilinear_fallback<D, T>(grid, pos);
    }

    template <typename D, typename T>
    inline D value_at_2(Grid<D, T> const& grid,
                        const Eigen::Matrix<T, 3, 1>& pos)
    {
        return value_at_2_quadratic(grid, pos);
        // Check if point is outside grid bounds
        if (pos.x() < grid.min().x() || pos.x() > grid.max().x()
            || pos.y() < grid.min().y() || pos.y() > grid.max().y()
            || pos.z() < grid.min().z() || pos.z() > grid.max().z())
        {

            // Find closest point on grid boundary
            Eigen::Matrix<T, 3, 1> clamped;
            clamped.x() = std::clamp(pos.x(), grid.min().x(), grid.max().x());
            clamped.y() = std::clamp(pos.y(), grid.min().y(), grid.max().y());
            clamped.z() = std::clamp(pos.z(), grid.min().z(), grid.max().z());

            D d_val = value_at_2(grid, clamped);

            T euclidean_dist = (pos - clamped).norm();

            return d_val >= 0 ? d_val + euclidean_dist : d_val - euclidean_dist;
        }

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
            Eigen::Matrix<size_t, 3, 1>(nodes1.x(), nodes0.y(), nodes0.z()));
        D const d010 = grid(
            Eigen::Matrix<size_t, 3, 1>(nodes0.x(), nodes1.y(), nodes0.z()));
        D const d011 = grid(
            Eigen::Matrix<size_t, 3, 1>(nodes1.x(), nodes1.y(), nodes0.z()));
        D const d100 = grid(
            Eigen::Matrix<size_t, 3, 1>(nodes0.x(), nodes0.y(), nodes1.z()));
        D const d101 = grid(
            Eigen::Matrix<size_t, 3, 1>(nodes1.x(), nodes0.y(), nodes1.z()));
        D const d110 = grid(
            Eigen::Matrix<size_t, 3, 1>(nodes0.x(), nodes1.y(), nodes1.z()));
        D const d111 = grid(
            Eigen::Matrix<size_t, 3, 1>(nodes1.x(), nodes1.y(), nodes1.z()));

        const Eigen::Matrix<T, 3, 1> stu
            = (pos
               - (nodes0.template cast<T>().cwiseProduct(grid.dir())
                  + grid.min()))
                  .cwiseQuotient(grid.dir());

        return trilinear<D, T>(d000, d001, d010, d011, d100, d101, d110, d111,
                               stu);
    }

    template <typename D, typename T>
    inline D value_at_2_get_closest_point(Grid<D, T> const& grid,
                                          const Eigen::Matrix<T, 3, 1>& pos,
                                          Eigen::Matrix<T, 3, 1>& closestPoint)
    {
        // Check if point is outside grid bounds
        if (pos.x() < grid.min().x() || pos.x() > grid.max().x()
            || pos.y() < grid.min().y() || pos.y() > grid.max().y()
            || pos.z() < grid.min().z() || pos.z() > grid.max().z())
        {

            // Find closest point on grid boundary
            Eigen::Matrix<T, 3, 1> clamped;
            clamped.x() = std::clamp(pos.x(), grid.min().x(), grid.max().x());
            clamped.y() = std::clamp(pos.y(), grid.min().y(), grid.max().y());
            clamped.z() = std::clamp(pos.z(), grid.min().z(), grid.max().z());

            // Get distance value at boundary point
            D d_val = value_at_2(grid, clamped);

            // Calculate Euclidean distance from boundary
            T euclidean_dist = (pos - clamped).norm();

            // Adjust signed distance based on boundary value
            return d_val >= 0 ? d_val + euclidean_dist : d_val - euclidean_dist;
        }

        // Original interpolation code for points inside grid
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
            Eigen::Matrix<size_t, 3, 1>(nodes1.x(), nodes0.y(), nodes0.z()));
        D const d010 = grid(
            Eigen::Matrix<size_t, 3, 1>(nodes0.x(), nodes1.y(), nodes0.z()));
        D const d011 = grid(
            Eigen::Matrix<size_t, 3, 1>(nodes1.x(), nodes1.y(), nodes0.z()));
        D const d100 = grid(
            Eigen::Matrix<size_t, 3, 1>(nodes0.x(), nodes0.y(), nodes1.z()));
        D const d101 = grid(
            Eigen::Matrix<size_t, 3, 1>(nodes1.x(), nodes0.y(), nodes1.z()));
        D const d110 = grid(
            Eigen::Matrix<size_t, 3, 1>(nodes0.x(), nodes1.y(), nodes1.z()));
        D const d111 = grid(
            Eigen::Matrix<size_t, 3, 1>(nodes1.x(), nodes1.y(), nodes1.z()));

        const Eigen::Matrix<T, 3, 1> stu
            = (pos
               - (nodes0.template cast<T>().cwiseProduct(grid.dir())
                  + grid.min()))
                  .cwiseQuotient(grid.dir());
        closestPoint = stu;

        return trilinear<D, T>(d000, d001, d010, d011, d100, d101, d110, d111,
                               stu);
    }
}

// GRID_VALUE_AT_H
#endif
