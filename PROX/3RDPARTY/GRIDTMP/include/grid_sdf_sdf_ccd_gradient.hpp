#ifndef GRID_SDF_SDF_CCD_GRADIENT_HPP
#define GRID_SDF_SDF_CCD_GRADIENT_HPP

#include <grid_grid.h>
#include <grid_sdf_sdf_voxelize.hpp>
#include <grid_local_optimization.hpp>

namespace SDFSDFContact
{
template <typename T>
EigenVector3<T> gridIndexToPos(int ix, int iy, int iz,
                               const EigenVector3<T>& global_mmin, T spacing)
{
    return {global_mmin.x() + T(ix) * spacing,
            global_mmin.y() + T(iy) * spacing,
            global_mmin.z() + T(iz) * spacing};
}

bool inBounds(int i, int n) { return i >= 0 && i < n; }

template <typename T> T length(const EigenVector3<T>& v)
{
    return sqrt(v.x() * v.x() + v.y() * v.y() + v.z() * v.z());
}

template <typename T> EigenVector3<T> normalize_safe(const EigenVector3<T>& v)
{
    T L = length(v);
    if (L <= 1e-15) return EigenVector3<T>(0, 0, 0);
    return v * (1.0 / L);
}

template <typename D, typename T>
void selectFeaturePointsPerVoxelOld(std::vector<SDFVoxel<T>>& voxels,
                                    const EigenVector3<T>& global_mmin,
                                    T spacing, int nx, int ny, int nz,
                                    grid::Grid<D, T> sampler)
{
    // Parameters / heuristics (tweakable)
    const T extremum_eps = 1e-8; // tolerance for strict extremum test
    const T min_discontinuity_threshold
        = spacing; // minimal absolute jump considered disruptive
    const T discontinuity_multiplier
        = 3.0; // how many times above mean diff is "large"
    const int global_max_per_voxel
        = 32; // hard cap on selected points for a voxel
    const int min_per_voxel = 1; // ensure at least this many points selected

    // offsets for 26-neighborhood
    std::vector<std::array<int, 3>> neighborOffsets;
    for (int dz = -1; dz <= 1; ++dz)
        for (int dy = -1; dy <= 1; ++dy)
            for (int dx = -1; dx <= 1; ++dx)
                if (!(dx == 0 && dy == 0 && dz == 0))
                    neighborOffsets.push_back({dx, dy, dz});

    // For convenience, pre-allocate a small cache for values inside a voxel:
    // We'll need to frequently re-query neighbor values possibly outside the voxel,
    // so we do not attempt any global cache - we simply call the sampler/grid-sample
    // on demand. But we will compute per-voxel data arrays for the points inside it.
    for (SDFVoxel<T>& v : voxels)
    {
        v.selected.clear();
        if (v.pointCount == 0) continue; // skip

        // First pass: iterate all points in voxel and gather sdf values and neighbor diffs
        struct GridPoint
        {
            int ix, iy, iz;
            EigenVector3<T> pos;
            T val;
        };

        std::vector<GridPoint> points;
        points.reserve(std::max(1u, v.pointCount));

        for (int iz = v.izMin; iz <= v.izMax; ++iz)
        {
            if (!inBounds(iz, nz)) continue;
            for (int iy = v.iyMin; iy <= v.iyMax; ++iy)
            {
                if (!inBounds(iy, ny)) continue;
                for (int ix = v.ixMin; ix <= v.ixMax; ++ix)
                {
                    if (!inBounds(ix, nx)) continue;
                    EigenVector3<T> p
                        = gridIndexToPos(ix, iy, iz, global_mmin, spacing);
                    T val = grid::value_at_2(sampler, p);
                    if (!isfinite(val)) continue; // skip NaNs/Infs
                    points.push_back({ix, iy, iz, p, val});
                }
            }
        }

        if (points.empty()) continue;

        // Precompute per-point neighbor statistics (max abs diff, mean abs diff, extreum score)
        struct PStats
        {
            T max_abs_diff = 0.0;
            T mean_abs_diff = 0.0;
            T extremum_score
                = 0.0; // how strong a peak/valley: min difference to neighbors
            int neighbor_count = 0;
            bool is_peak = false;
            bool is_valley = false;
        };

        std::vector<PStats> stats(points.size());

        // Build a quick lambda to get neighbor value (if neighbor index out-of-bounds return NaN)
        auto gridValue = [&](int gx, int gy, int gz) -> std::pair<bool, T>
        {
            if (!inBounds(gx, nx) || !inBounds(gy, ny) || !inBounds(gz, nz))
                return {false, 0.0};
            EigenVector3<T> p
                = gridIndexToPos(gx, gy, gz, global_mmin, spacing);
            T v = grid::value_at_2(sampler, p);
            if (!isfinite(v)) return {false, 0.0};
            return {true, v};
        };

        for (size_t pi = 0; pi < points.size(); ++pi)
        {
            const auto& gp = points[pi];
            PStats s;
            T acc = 0.0;
            int cnt = 0;
            T local_min_diff = std::numeric_limits<T>::infinity();
            bool all_lower
                = true; // candidate for peak: all neighbors lower (by eps)
            bool all_higher
                = true; // candidate for valley: all neighbors higher (by eps)

            for (auto& off : neighborOffsets)
            {
                int nxg = gp.ix + off[0];
                int nyg = gp.iy + off[1];
                int nzg = gp.iz + off[2];
                auto [ok, nv] = gridValue(nxg, nyg, nzg);
                if (!ok) continue;
                T diff = gp.val - nv;
                T ad = fabs(diff);
                acc += ad;
                ++cnt;
                local_min_diff = std::min(local_min_diff, ad);
                if (!(diff > extremum_eps))
                    all_lower = false; // not strictly greater than neighbor
                if (!(diff < -extremum_eps))
                    all_higher = false; // not strictly lower than neighbor
            }
            s.neighbor_count = cnt;
            s.mean_abs_diff = (cnt > 0) ? (acc / T(cnt)) : 0.0;
            s.max_abs_diff = 0.0;
            // recompute max_abs_diff explicitly (we want max not mean)
            T maxad = 0.0;
            for (auto& off : neighborOffsets)
            {
                int nxg = gp.ix + off[0];
                int nyg = gp.iy + off[1];
                int nzg = gp.iz + off[2];
                auto [ok, nv] = gridValue(nxg, nyg, nzg);
                if (!ok) continue;
                maxad = std::max(maxad, fabs(gp.val - nv));
            }
            s.max_abs_diff = maxad;
            s.extremum_score
                = (local_min_diff == std::numeric_limits<T>::infinity())
                    ? 0.0
                    : local_min_diff;
            s.is_peak = all_lower && (s.neighbor_count > 0);
            s.is_valley = all_higher && (s.neighbor_count > 0);
            stats[pi] = s;
        }

        // Stage 1: collect extrema (peaks and valleys)
        std::vector<std::pair<T, int>> extrema_candidates; // (score, index)
        for (size_t i = 0; i < points.size(); ++i)
        {
            if (stats[i].is_peak)
            {
                // score: how strongly above neighbors: use extremum_score
                extrema_candidates.emplace_back(stats[i].extremum_score,
                                                (int)i);
            }
            else if (stats[i].is_valley)
            {
                // valley score also positive (how strongly below neighbors)
                extrema_candidates.emplace_back(stats[i].extremum_score,
                                                (int)i);
            }
        }

        // sort descending by score (strongest extrema first)
        sort(extrema_candidates.begin(), extrema_candidates.end(),
             [](const std::pair<T, int>& a, const std::pair<T, int>& b)
             { return a.first > b.first; });

        int cap = std::min(global_max_per_voxel,
                           std::max(min_per_voxel, int(v.pointCount / 8)));
        // Accept extrema up to cap (we'll later possibly add discontinuities if room or do fallback)
        for (size_t ei = 0;
             ei < extrema_candidates.size() && (int)v.selected.size() < cap;
             ++ei)
        {
            int idx = extrema_candidates[ei].second;
            SelectedPoint<T> sp;
            sp.pos = points[idx].pos;
            sp.sdf = points[idx].val;
            sp.gx = points[idx].ix;
            sp.gy = points[idx].iy;
            sp.gz = points[idx].iz;
            if (stats[idx].is_peak)
                sp.type = SelectedPoint<T>::PEAK;
            else
                sp.type = SelectedPoint<T>::VALLEY;
            v.selected.push_back(sp);
        }

        // Stage 2: if we did not get "enough" (or want additional important points),
        // find discontinuities. We compute a voxel-typical scale and pick points with
        // max_abs_diff significantly larger than that.
        if ((int)v.selected.size() < cap)
        {
            // compute typical mean over points of mean_abs_diff
            T mean_mean_diff = 0.0;
            for (size_t i = 0; i < points.size(); ++i)
                mean_mean_diff += stats[i].mean_abs_diff;
            mean_mean_diff /= T(points.size());
            T threshold = std::max(min_discontinuity_threshold,
                                   discontinuity_multiplier * mean_mean_diff);

            // gather candidates with their jump severity
            std::vector<std::pair<T, int>> disc_cands; // (max_abs_diff, index)
            for (size_t i = 0; i < points.size(); ++i)
            {
                // skip if already selected (same grid index)
                bool already = false;
                for (auto& sp : v.selected)
                {
                    if (sp.gx == points[i].ix && sp.gy == points[i].iy
                        && sp.gz == points[i].iz)
                    {
                        already = true;
                        break;
                    }
                }
                if (already) continue;
                if (stats[i].max_abs_diff >= threshold)
                {
                    disc_cands.emplace_back(stats[i].max_abs_diff, (int)i);
                }
            }
            // sort descending by jump magnitude
            std::sort(disc_cands.begin(), disc_cands.end(),
                      [](const std::pair<T, int>& a, const std::pair<T, int>& b)
                      { return a.first > b.first; });

            for (size_t di = 0;
                 di < disc_cands.size() && (int)v.selected.size() < cap; ++di)
            {
                int idx = disc_cands[di].second;
                SelectedPoint<T> sp;
                sp.pos = points[idx].pos;
                sp.sdf = points[idx].val;
                sp.gx = points[idx].ix;
                sp.gy = points[idx].iy;
                sp.gz = points[idx].iz;
                sp.type = SelectedPoint<T>::DISCONTINUITY;
                v.selected.push_back(sp);
            }
        }

        // Stage 3: if still nothing selected, create a representative isosurface point
        if (v.selected.empty())
        {
            // Compute center of voxel (world coords)
            EigenVector3<T> center
                = EigenVector3<T>((v.vmin.x() + v.vmax.x()) * 0.5,
                                  (v.vmin.y() + v.vmax.y()) * 0.5,
                                  (v.vmin.z() + v.vmax.z()) * 0.5);
            T sdf_center = grid::value_at_2(sampler, center);
            if (!isfinite(sdf_center)) sdf_center = 1e9; // avoid NaN

            // If center already on isosurface (within small tol) take it
            const T iso_eps = 1e-8;
            if (fabs(sdf_center) <= iso_eps)
            {
                SelectedPoint<T> sp;
                sp.pos = center;
                sp.sdf = sdf_center;
                sp.type = SelectedPoint<T>::CENTER_ISOSURFACE;
                v.selected.push_back(sp);
            }
            else
            {
                // Estimate gradient and project to isosurface: p_iso = center - sdf_center * grad/|grad|
                T h = spacing * 0.5;

                EigenVector3<T> grad
                    = grid::computeGradient_Working(center, sampler);
                T gnorm = length(grad);
                if (gnorm > 1e-12)
                {
                    EigenVector3<T> n = normalize_safe(grad);
                    EigenVector3<T> p_iso = center - n * sdf_center;
                    SelectedPoint<T> sp;
                    sp.pos = p_iso;
                    sp.sdf = grid::value_at_2(
                        sampler,
                        p_iso); //sampler(p_iso); // should be close to zero
                    sp.type = SelectedPoint<T>::PROJECTED_ISOSURFACE;
                    v.selected.push_back(sp);
                }
                else
                {
                    // If gradient too small, fallback to selecting the closest grid point in this voxel
                    T best_abs = std::numeric_limits<T>::infinity();
                    size_t best_idx = 0;
                    for (size_t i = 0; i < points.size(); ++i)
                    {
                        T ab = fabs(points[i].val);
                        if (ab < best_abs)
                        {
                            best_abs = ab;
                            best_idx = i;
                        }
                    }
                    SelectedPoint<T> sp;
                    sp.pos = points[best_idx].pos;
                    sp.sdf = points[best_idx].val;
                    sp.gx = points[best_idx].ix;
                    sp.gy = points[best_idx].iy;
                    sp.gz = points[best_idx].iz;
                    sp.type = SelectedPoint<T>::NEAREST_GRID;
                    v.selected.push_back(sp);
                }
            }
        } // end stage 3

        // Final cap enforcement: if we somehow still exceed cap, keep strongest by priority:
        // Priority order: PEAK/VALLEY (extrema) > DISCONTINUITY > PROJECTED/CENTER/NEAREST
        if ((int)v.selected.size() > cap)
        {
            auto priority = [](const SelectedPoint<T>& sp) -> int
            {
                switch (sp.type)
                {
                case SelectedPoint<T>::PEAK:
                case SelectedPoint<T>::VALLEY:        return 3;
                case SelectedPoint<T>::DISCONTINUITY: return 2;
                default:                              return 1;
                }
            };
            std::sort(v.selected.begin(), v.selected.end(),
                      [&](const SelectedPoint<T>& a, const SelectedPoint<T>& b)
                      {
                          int pa = priority(a), pb = priority(b);
                          if (pa != pb) return pa > pb;
                          // tie-breaker: prefer smaller absolute sdf (closer to iso) for projection/nearest
                          return std::fabs(a.sdf) < std::fabs(b.sdf);
                      });
            v.selected.resize(cap);
        }

    } // end loop voxels
}

template <typename T>
std::vector<SelectedPoint<T>>
filterClosePointsStable(const std::vector<SelectedPoint<T>>& pts, T minDist)
{
    std::vector<SelectedPoint<T>> kept;
    T minDist2 = minDist * minDist;

    for (const auto& p : pts)
    {
        bool tooClose = false;
        for (const auto& k : kept)
        {
            if ((p.pos - k.pos).squaredNorm() < minDist2)
            { // squaredNorm -> no sqrt
                tooClose = true;
                break;
            }
        }
        if (!tooClose) kept.push_back(p);
    }
    return kept;
}

template <typename D, typename T>
void selectFeaturePointsPerVoxel(std::vector<SDFVoxel<T>>& voxels,
                                 const EigenVector3<T>& global_mmin, T spacing,
                                 int nx, int ny, int nz, grid::Grid<D, T> grid)
{
    // Heuristics / parameters (tweak to taste)
    const T extremum_eps = 1e-8; // tolerance for strict comparisons
    const T discontinuity_multiplier
        = 3.0; // how many times above mean neighbor-dist
    const T min_discontinuity_threshold = spacing * 0.5;
    const int global_max_per_voxel = 32;
    const int min_per_voxel = 1;

    // neighbor offsets (26-neighborhood)
    std::vector<std::array<int, 3>> neighborOffsets;
    for (int dz = -1; dz <= 1; ++dz)
        for (int dy = -1; dy <= 1; ++dy)
            for (int dx = -1; dx <= 1; ++dx)
                if (!(dx == 0 && dy == 0 && dz == 0))
                    neighborOffsets.push_back({dx, dy, dz});

    auto inBounds = [](int i, int n) { return (i >= 0 && i < n); };

    // world position from grid index
    auto gridIndexToPos = [&](int ix, int iy, int iz) -> EigenVector3<T>
    {
        return EigenVector3<T>{global_mmin.x() + ix * spacing,
                               global_mmin.y() + iy * spacing,
                               global_mmin.z() + iz * spacing};
    };

    auto vecLen = [&](const EigenVector3<T>& a) -> T
    { return std::sqrt(a.x() * a.x() + a.y() * a.y() + a.z() * a.z()); };
    auto normalize = [&](const EigenVector3<T>& a) -> EigenVector3<T>
    {
        T L = vecLen(a);
        if (L <= 1e-15) return EigenVector3<T>{0.0, 0.0, 0.0};
        return EigenVector3<T>{a.x() / L, a.y() / L, a.z() / L};
    };

    // Helper: attempt to project point p (world coords) to isosurface using local gradient.
    // If gradient is unreliable, attempt simple sign-change linear interpolation along axis-neighbors.
    // Returns pair<ok, p_iso>. Also returns sdf at original point via out_sdf.
    auto projectToIso = [&](const EigenVector3<T>& p, T out_sdf, int gx, int gy,
                            int gz) -> std::pair<bool, EigenVector3<T>>
    {
        // First try gradient projection
        T h = spacing * 0.5;
        EigenVector3<T> grad = grid::computeGradient_Working(p, grid);
        T gnorm = vecLen(grad);
        if (gnorm > 1e-12)
        {
            EigenVector3<T> n = normalize(grad);
            EigenVector3<T> p_iso = EigenVector3<T>{p.x() - out_sdf * n.x(),
                                                    p.y() - out_sdf * n.y(),
                                                    p.z() - out_sdf * n.z()};
            return {true, p_iso};
        }
        // Fallback: try to find neighbor with opposite sign and linear-interpolate to zero along that neighbor vector.
        // We use axis-aligned 6-neighborhood for sign-changes to keep it simple.
        std::array<int, 3> axes[6] = {
            std::array<int, 3>{ 1,  0,  0},
             std::array<int, 3>{-1,  0,  0},
            std::array<int, 3>{ 0,  1,  0},
             std::array<int, 3>{ 0, -1,  0},
            std::array<int, 3>{ 0,  0,  1},
             std::array<int, 3>{ 0,  0, -1}
        };
        for (int a = 0; a < 6; ++a)
        {
            int ni = gx + axes[a][0];
            int nj = gy + axes[a][1];
            int nk = gz + axes[a][2];
            if (!inBounds(ni, nx) || !inBounds(nj, ny) || !inBounds(nk, nz))
                continue;
            EigenVector3<T> qpos = gridIndexToPos(ni, nj, nk);
            T qval = grid::value_at_2(grid, qpos);
            if (!std::isfinite(qval)) continue;
            if (out_sdf * qval < 0.0)
            {
                // sign change: interpolate between p (value out_sdf) and qpos (value qval)
                T t = out_sdf / (out_sdf - qval); // t in [0,1]
                // p + t*(qpos - p) -> point where linear interpolation crosses zero
                EigenVector3<T> p_iso{p.x() + t * (qpos.x() - p.x()),
                                      p.y() + t * (qpos.y() - p.y()),
                                      p.z() + t * (qpos.z() - p.z())};
                return {true, p_iso};
            }
        }
        // failed
        return {
            false, EigenVector3<T>{0, 0, 0}
        };
    };

    // For each voxel, produce per-grid-point projected isosurface points, then analyze the projected points
    for (std::size_t vi = 0; vi < voxels.size(); ++vi)
    {
        SDFVoxel<T>& v = voxels[vi];
        v.selected.clear();
        if (v.pointCount == 0) continue;

        // gather grid points inside voxel
        struct GP
        {
            int ix, iy, iz;
            EigenVector3<T> p;
            T sdf;
        };

        std::vector<GP> pts;
        pts.reserve(std::max(1U, v.pointCount));
        for (int iz = v.izMin; iz <= v.izMax; ++iz)
        {
            if (!inBounds(iz, nz)) continue;
            for (int iy = v.iyMin; iy <= v.iyMax; ++iy)
            {
                if (!inBounds(iy, ny)) continue;
                for (int ix = v.ixMin; ix <= v.ixMax; ++ix)
                {
                    if (!inBounds(ix, nx)) continue;
                    EigenVector3<T> p = gridIndexToPos(ix, iy, iz);
                    T val = grid::value_at_2(grid, p);
                    if (!std::isfinite(val)) continue;
                    pts.push_back({ix, iy, iz, p, val});
                }
            }
        }
        if (pts.empty()) continue;

        // project each sample to isosurface (store only projected points that succeeded)
        struct Projection
        {
            bool ok;
            EigenVector3<T> p_iso;
            EigenVector3<T> normal;
            T orig_sdf;
        };

        std::vector<Projection> projs;
        projs.resize(pts.size());
        // also remember which projections succeeded for neighbor lookups
        for (std::size_t i = 0; i < pts.size(); ++i)
        {
            const GP& g = pts[i];
            // try to compute gradient and project
            T sdf_val = g.sdf;
            auto pr = projectToIso(g.p, sdf_val, g.ix, g.iy, g.iz);
            if (pr.first)
            {
                // compute normal at projected point if possible (for later use)
                EigenVector3<T> grad_at_orig
                    = grid::computeGradient_Working(g.p, grid);
                T gn = vecLen(grad_at_orig);
                EigenVector3<T> n = (gn > 1e-12) ? normalize(grad_at_orig)
                                                 : EigenVector3<T>{0, 0, 0};
                projs[i] = Projection{true, pr.second, n, sdf_val};
            }
            else
            {
                // projection failed; mark false
                projs[i] = Projection{
                    false, EigenVector3<T>{0, 0, 0},
                     EigenVector3<T>{0, 0, 0},
                    sdf_val
                };
            }
        }

        // compute neighbor statistics on projected points (distances and coordinate comparisons)
        // We'll compute for each valid projection:
        //  - axis-wise comparisons vs projected neighbors (for peak/valley detection)
        //  - max projected distance to any projected neighbor (for discontinuity)
        std::vector<bool> is_peak(pts.size(), false),
            is_valley(pts.size(), false);
        std::vector<T> max_proj_dist(pts.size(), 0.0),
            mean_proj_dist(pts.size(), 0.0);
        std::vector<int> proj_neighbor_count(pts.size(), 0);

        // build map from grid index tuple -> index in pts to look up neighbors quickly
        std::unordered_map<long long, std::size_t> idxMap;
        idxMap.reserve(pts.size() * 2);
        auto keyOf = [&](int ix, int iy, int iz) -> long long
        {
            // pack into 64-bit: ix (21 bits), iy (21), iz (21) — safe for reasonable grid sizes
            return ((long long)ix << 42) ^ ((long long)iy << 21)
                 ^ (long long)iz;
        };
        for (std::size_t i = 0; i < pts.size(); ++i)
        {
            idxMap[keyOf(pts[i].ix, pts[i].iy, pts[i].iz)] = i;
        }

        // For each projected point, examine neighbors' projected positions (if neighbors projected ok)
        for (std::size_t i = 0; i < pts.size(); ++i)
        {
            if (!projs[i].ok) continue;
            EigenVector3<T> P = projs[i].p_iso;
            T accd = 0.0;
            int cnt = 0;
            T maxd = 0.0;
            // gather neighbor projected positions
            std::vector<EigenVector3<T>> neighbor_proj_positions;
            neighbor_proj_positions.reserve(neighborOffsets.size());
            for (auto& off : neighborOffsets)
            {
                int nx_i = pts[i].ix + off[0];
                int ny_i = pts[i].iy + off[1];
                int nz_i = pts[i].iz + off[2];
                if (!inBounds(nx_i, nx) || !inBounds(ny_i, ny)
                    || !inBounds(nz_i, nz))
                    continue;
                auto it = idxMap.find(keyOf(nx_i, ny_i, nz_i));
                if (it == idxMap.end()) continue;
                std::size_t nidx = it->second;
                if (!projs[nidx].ok) continue;
                EigenVector3<T> Q = projs[nidx].p_iso;
                T d = vecLen(EigenVector3<T>{P.x() - Q.x(), P.y() - Q.y(),
                                             P.z() - Q.z()});
                accd += d;
                ++cnt;
                if (d > maxd) maxd = d;
                neighbor_proj_positions.push_back(Q);
            }
            proj_neighbor_count[i] = cnt;
            mean_proj_dist[i] = (cnt > 0) ? (accd / T(cnt)) : 0.0;
            max_proj_dist[i] = maxd;

            // Axis-wise extrema check (compare P.x vs neighbors' x, and similarly y,z)
            // We declare a peak if P.[axis] > all neighbors.[axis] + eps
            // and a valley if P.[axis] < all neighbors.[axis] - eps
            bool anyPeak = false, anyValley = false;
            if (cnt > 0)
            {
                // Check for each axis independently; if any axis yields strict extremum, mark it.
                // This is a simple, conservative criteria for "geometric" peak/valley on the surface.
                bool peakX = true, valleyX = true, peakY = true, valleyY = true,
                     peakZ = true, valleyZ = true;
                for (const EigenVector3<T>& Q : neighbor_proj_positions)
                {
                    if (!(P.x() > Q.x() + extremum_eps)) peakX = false;
                    if (!(P.x() < Q.x() - extremum_eps)) valleyX = false;
                    if (!(P.y() > Q.y() + extremum_eps)) peakY = false;
                    if (!(P.y() < Q.y() - extremum_eps)) valleyY = false;
                    if (!(P.z() > Q.z() + extremum_eps)) peakZ = false;
                    if (!(P.z() < Q.z() - extremum_eps)) valleyZ = false;
                }
                anyPeak = (peakX || peakY || peakZ);
                anyValley = (valleyX || valleyY || valleyZ);
            }
            is_peak[i] = anyPeak;
            is_valley[i] = anyValley;
        }

        // Stage 1: collect extrema candidates (based on projected positions)
        std::vector<std::pair<T, int>> extrema_candidates; // score, idx
        for (std::size_t i = 0; i < pts.size(); ++i)
        {
            if (projs[i].ok && (is_peak[i] || is_valley[i]))
            {
                // Use mean_proj_dist (how well-connected) as score: larger means stronger isolation
                T score = mean_proj_dist[i];
                extrema_candidates.emplace_back(score, (int)i);
            }
        }
        std::sort(extrema_candidates.begin(), extrema_candidates.end(),
                  [](const std::pair<T, int>& a, const std::pair<T, int>& b)
                  { return a.first > b.first; });

        // determine cap
        int cap = std::min(global_max_per_voxel,
                           std::max(min_per_voxel, int(v.pointCount / 8)));
        std::vector<char> already_selected(pts.size(), 0);
        for (std::size_t ei = 0;
             ei < extrema_candidates.size() && (int)v.selected.size() < cap;
             ++ei)
        {
            int idxPt = extrema_candidates[ei].second;
            SelectedPoint<T> sp;
            sp.pos = projs[idxPt].p_iso;
            sp.sdf = projs[idxPt].orig_sdf;
            sp.gx = pts[idxPt].ix;
            sp.gy = pts[idxPt].iy;
            sp.gz = pts[idxPt].iz;
            sp.type = is_peak[idxPt] ? SelectedPoint<T>::PEAK
                                     : SelectedPoint<T>::VALLEY;
            v.selected.push_back(sp);
            already_selected[idxPt] = 1;
        }

        // Stage 2: discontinuities based on projected-point jumps
        if ((int)v.selected.size() < cap)
        {
            // compute voxel-typical mean of mean_proj_dist for thresholding
            T mean_mean = 0.0;
            int cnt_valid = 0;
            for (std::size_t i = 0; i < pts.size(); ++i)
            {
                if (!projs[i].ok) continue;
                mean_mean += mean_proj_dist[i];
                ++cnt_valid;
            }
            mean_mean = (cnt_valid > 0) ? mean_mean / T(cnt_valid) : 0.0;
            T threshold = std::max<T>(min_discontinuity_threshold,
                                      discontinuity_multiplier * mean_mean);

            std::vector<std::pair<T, int>> disc_cands;
            for (std::size_t i = 0; i < pts.size(); ++i)
            {
                if (!projs[i].ok) continue;
                if (already_selected[i]) continue;
                if (max_proj_dist[i] >= threshold)
                    disc_cands.emplace_back(max_proj_dist[i], (int)i);
            }
            std::sort(disc_cands.begin(), disc_cands.end(),
                      [](const std::pair<double, int>& a,
                         const std::pair<double, int>& b)
                      { return a.first > b.first; });

            for (std::size_t di = 0;
                 di < disc_cands.size() && (int)v.selected.size() < cap; ++di)
            {
                int idxPt = disc_cands[di].second;
                SelectedPoint<T> sp;
                sp.pos = projs[idxPt].p_iso;
                sp.sdf = projs[idxPt].orig_sdf;
                sp.gx = pts[idxPt].ix;
                sp.gy = pts[idxPt].iy;
                sp.gz = pts[idxPt].iz;
                sp.type = SelectedPoint<T>::DISCONTINUITY;
                v.selected.push_back(sp);
                already_selected[idxPt] = 1;
            }
        }

        // Stage 3: fallback — if nothing selected, create representative isosurface point for this voxel
        if (v.selected.empty())
        {
            EigenVector3<T> center{0.5 * (v.vmin.x() + v.vmax.x()),
                                   0.5 * (v.vmin.y() + v.vmax.y()),
                                   0.5 * (v.vmin.z() + v.vmax.z())};
            T sc = grid::value_at_2(grid, center);
            if (std::isfinite(sc) && std::fabs(sc) <= 1e-8)
            {
                SelectedPoint<T> sp;
                sp.pos = center;
                sp.sdf = sc;
                sp.type = SelectedPoint<T>::CENTER_ISOSURFACE;
                sp.gx = sp.gy = sp.gz = -1;
                v.selected.push_back(sp);
            }
            else
            {
                // attempt gradient projection at center
                EigenVector3<T> grad
                    = grid::computeGradient_Working(center, grid);
                T gnorm = vecLen(grad);
                if (gnorm > 1e-12 && std::isfinite(sc))
                {
                    EigenVector3<T> n = normalize(grad);
                    EigenVector3<T> p_iso{center.x() - sc * n.x(),
                                          center.y() - sc * n.y(),
                                          center.z() - sc * n.z()};
                    SelectedPoint<T> sp;
                    sp.pos = p_iso;
                    sp.sdf = grid::value_at_2(grid, p_iso);
                    sp.type = SelectedPoint<T>::PROJECTED_ISOSURFACE;
                    sp.gx = sp.gy = sp.gz = -1;
                    v.selected.push_back(sp);
                }
                else
                {
                    // fallback: pick voxel grid point closest to zero
                    T bestAbs = 1e99;
                    int bestIdx = -1;
                    for (std::size_t i = 0; i < pts.size(); ++i)
                    {
                        T a = std::fabs(pts[i].sdf);
                        if (a < bestAbs)
                        {
                            bestAbs = a;
                            bestIdx = (int)i;
                        }
                    }
                    if (bestIdx >= 0)
                    {
                        // project that nearest grid point (we already computed projs for all)
                        if (projs[bestIdx].ok)
                        {
                            SelectedPoint<T> sp;
                            sp.pos = projs[bestIdx].p_iso;
                            sp.sdf = projs[bestIdx].orig_sdf;
                            sp.type = SelectedPoint<T>::NEAREST_GRID;
                            sp.gx = pts[bestIdx].ix;
                            sp.gy = pts[bestIdx].iy;
                            sp.gz = pts[bestIdx].iz;
                            v.selected.push_back(sp);
                        }
                        else
                        {
                            SelectedPoint<T> sp;
                            sp.pos = pts[bestIdx].p;
                            sp.sdf = pts[bestIdx].sdf;
                            sp.type = SelectedPoint<T>::NEAREST_GRID;
                            sp.gx = pts[bestIdx].ix;
                            sp.gy = pts[bestIdx].iy;
                            sp.gz = pts[bestIdx].iz;
                            v.selected.push_back(sp);
                        }
                    }
                }
            }
        }

        // final enforce cap and ordering: prefer extrema > discontinuity > projected/fallback
        if ((int)v.selected.size() > cap)
        {
            auto priority = [](const SelectedPoint<T>& sp) -> int
            {
                switch (sp.type)
                {
                case SelectedPoint<T>::PEAK:
                case SelectedPoint<T>::VALLEY:        return 3;
                case SelectedPoint<T>::DISCONTINUITY: return 2;
                default:                              return 1;
                }
            };
            std::sort(v.selected.begin(), v.selected.end(),
                      [&](const SelectedPoint<T>& a, const SelectedPoint<T>& b)
                      {
                          int pa = priority(a), pb = priority(b);
                          if (pa != pb) return pa > pb;
                          return std::fabs(a.sdf) < std::fabs(b.sdf);
                      });
            v.selected.resize(cap);
        }
    } // end voxels loop

    //Project all points to surface
    for (size_t i = 0; i < voxels.size(); ++i)
    {
        for (size_t j = 0; j < voxels[i].selected.size(); ++j)
        {
            //Project ALL points to surface!!!!! We are in local space so we
            // directly leverage all points
            voxels[i].selected[j].pos
                = voxels[i].selected[j].pos
                - grid::value_at_2<T, T>(grid, voxels[i].selected[j].pos)
                      * grid::computeGradient_Working<T, T>(
                          voxels[i].selected[j].pos, grid);
        }
    }

    //Dedupe all points
    std::cerr << "Deduping points...\n";
    T totalPoints = 0;
    for (size_t i = 0; i < voxels.size(); ++i)
    {
        T voxelSize = voxels[i].selected.size();
        totalPoints += voxelSize;
    }
    std::cerr << "Total voxels before dedup: " << totalPoints << "\n";
    for (size_t i = 0; i < voxels.size(); ++i)
    {
        std::vector<SelectedPoint<T>> val
            = filterClosePointsStable<T>(voxels[i].selected, 0.1);
        voxels[i].selected = val;
    }
    totalPoints = 0;
    for (size_t i = 0; i < voxels.size(); ++i)
    {
        T voxelSize = voxels[i].selected.size();
        totalPoints += voxelSize;
    }
    std::cerr << "Total voxels after dedup: " << totalPoints << "\n";
}

template <typename D, typename T>
void selectFeaturePointsPerVoxel2(std::vector<SDFVoxel<T>>& voxels,
                                  const EigenVector3<T>& global_mmin, T spacing,
                                  int nx, int ny, int nz,
                                  const grid::Grid<D, T>& grid)
{
    // Tunable parameters (safe defaults)
    const int SUBRES
        = 5; // samples per axis inside each voxel (5^3 = 125 samples)
    const double extremum_eps
        = 1e-8; // tolerance for strict extremum test on interpolated samples
    const double discontinuity_multiplier = 3.0;
    const double min_disc_threshold = spacing * 0.5;
    const int global_max_per_voxel = 32;
    const int min_per_voxel = 1;

    // Precompute neighbor offsets for the dense subgrid (26 neighbors)
    std::vector<std::array<int, 3>> neighOffsets;
    for (int dz = -1; dz <= 1; ++dz)
        for (int dy = -1; dy <= 1; ++dy)
            for (int dx = -1; dx <= 1; ++dx)
                if (!(dx == 0 && dy == 0 && dz == 0))
                    neighOffsets.push_back({dx, dy, dz});

    // Helper to compute a world position inside a voxel at subgrid index (sx,sy,sz)
    auto samplePosInVoxel
        = [&](const SDFVoxel<T>& v, int sx, int sy, int sz) -> EigenVector3<T>
    {
        // sx/ (SUBRES-1) ranges 0..1 inclusive
        double fx = (SUBRES == 1) ? 0.5 : double(sx) / double(SUBRES - 1);
        double fy = (SUBRES == 1) ? 0.5 : double(sy) / double(SUBRES - 1);
        double fz = (SUBRES == 1) ? 0.5 : double(sz) / double(SUBRES - 1);
        return EigenVector3<T>(v.vmin.x() + fx * (v.vmax.x() - v.vmin.x()),
                               v.vmin.y() + fy * (v.vmax.y() - v.vmin.y()),
                               v.vmin.z() + fz * (v.vmax.z() - v.vmin.z()));
    };

    // Loop over voxels
    for (SDFVoxel<T>& v : voxels)
    {
        v.selected.clear();
        if (v.pointCount == 0) continue;

        // Build dense sample grid inside voxel and evaluate sampler at each sample
        // We'll store values in a 1D vector with index formula sidx = (sx*SUBRES + sy)*SUBRES + sz
        const int S = SUBRES;
        const int S3 = S * S * S;
        std::vector<EigenVector3<T>> samplePos;
        samplePos.reserve(S3);
        std::vector<T> sampleVal;
        sampleVal.reserve(S3);
        for (int sx = 0; sx < S; ++sx)
        {
            for (int sy = 0; sy < S; ++sy)
            {
                for (int sz = 0; sz < S; ++sz)
                {
                    EigenVector3<T> p = samplePosInVoxel(v, sx, sy, sz);
                    T val = grid::value_at_2(grid, p);
                    // sampler is assumed to perform interpolation if needed
                    samplePos.push_back(p);
                    sampleVal.push_back(val);
                }
            }
        }

        // Helper to convert (sx,sy,sz) -> flat index
        auto sIndex = [&](int sx, int sy, int sz) -> int
        { return (sx * S + sy) * S + sz; };
        auto insideSample = [&](int sx, int sy, int sz) -> bool
        { return sx >= 0 && sx < S && sy >= 0 && sy < S && sz >= 0 && sz < S; };

        // Compute neighbor statistics on the dense samples
        std::vector<T> meanAbsDiff(S3, 0.0);
        std::vector<T> maxAbsDiff(S3, 0.0);
        std::vector<bool> is_peak(S3, false), is_valley(S3, false);

        for (int sx = 0; sx < S; ++sx)
        {
            for (int sy = 0; sy < S; ++sy)
            {
                for (int sz = 0; sz < S; ++sz)
                {
                    int idx = sIndex(sx, sy, sz);
                    double vcenter = sampleVal[idx];
                    double acc = 0.0;
                    int cnt = 0;
                    double maxd = 0.0;
                    bool all_lower = true;
                    bool all_higher = true;
                    for (auto& off : neighOffsets)
                    {
                        int nx_s = sx + off[0];
                        int ny_s = sy + off[1];
                        int nz_s = sz + off[2];
                        if (!insideSample(nx_s, ny_s, nz_s)) continue;
                        int nidx = sIndex(nx_s, ny_s, nz_s);
                        double nval = sampleVal[nidx];
                        double d = vcenter - nval;
                        acc += fabs(d);
                        ++cnt;
                        maxd = std::max(maxd, fabs(d));
                        if (!(d > extremum_eps)) all_lower = false;
                        if (!(d < -extremum_eps)) all_higher = false;
                    }
                    meanAbsDiff[idx] = (cnt > 0) ? acc / double(cnt) : 0.0;
                    maxAbsDiff[idx] = maxd;
                    is_peak[idx] = all_lower && (cnt > 0);
                    is_valley[idx] = all_higher && (cnt > 0);
                }
            }
        }

        // Stage 1: collect extrema (peaks and valleys) from dense samples (strongest first)
        std::vector<std::pair<double, int>> extremaCandidates; // (score, index)
        for (int i = 0; i < S3; ++i)
        {
            if (is_peak[i] || is_valley[i])
            {
                // Use min neighbor difference as score approximated by meanAbsDiff (good heuristic)
                extremaCandidates.emplace_back(meanAbsDiff[i], i);
            }
        }
        std::sort(extremaCandidates.begin(), extremaCandidates.end(),
                  [](const std::pair<T, int>& a, const std::pair<T, int>& b)
                  { return a.first > b.first; });

        int cap = std::min(global_max_per_voxel,
                           std::max(min_per_voxel, int(v.pointCount / 8)));
        // Accept extrema up to cap (avoid duplicates by position)
        auto addSelected =
            [&](const EigenVector3<T>& pos, T sdfv, SelectedPoint<T>::Type type)
        {
            // avoid duplicates (by close distance)
            const T dup_eps = spacing * 1e-6;
            for (auto& s : v.selected)
            {
                T dx = s.pos.x() - pos.x(), dy = s.pos.y() - pos.y(),
                  dz = s.pos.z() - pos.z();
                if (dx * dx + dy * dy + dz * dz <= dup_eps * dup_eps) return;
            }
            SelectedPoint<T> sp;
            sp.pos = pos;
            sp.sdf = sdfv;
            // map enum (we keep same enum names as original)
            if (type == SelectedPoint<T>::PEAK)
                sp.type = SelectedPoint<T>::PEAK;
            else if (type == SelectedPoint<T>::VALLEY)
                sp.type = SelectedPoint<T>::VALLEY;
            else if (type == SelectedPoint<T>::DISCONTINUITY)
                sp.type = SelectedPoint<T>::DISCONTINUITY;
            else if (type == SelectedPoint<T>::PROJECTED_ISOSURFACE)
                sp.type = SelectedPoint<T>::PROJECTED_ISOSURFACE;
            else if (type == SelectedPoint<T>::CENTER_ISOSURFACE)
                sp.type = SelectedPoint<T>::CENTER_ISOSURFACE;
            else
                sp.type = SelectedPoint<T>::NEAREST_GRID;
            v.selected.push_back(sp);
        };

        for (size_t ei = 0;
             ei < extremaCandidates.size() && (int)v.selected.size() < cap;
             ++ei)
        {
            int sidx = extremaCandidates[ei].second;
            EigenVector3<T> pos = samplePos[sidx];
            double val = sampleVal[sidx];
            if (is_peak[sidx])
                addSelected(pos, val, SelectedPoint<T>::PEAK);
            else
                addSelected(pos, val, SelectedPoint<T>::VALLEY);
        }

        // Stage 2: discontinuities on dense samples (adaptive threshold)
        if ((int)v.selected.size() < cap)
        {
            double mean_mean = 0.0;
            for (int i = 0; i < S3; ++i) mean_mean += meanAbsDiff[i];
            mean_mean /= double(S3);
            double threshold = std::max(min_disc_threshold,
                                        discontinuity_multiplier * mean_mean);

            std::vector<std::pair<T, int>> discCands;
            for (int i = 0; i < S3; ++i)
            {
                // skip if already picked (pos duplicate check via sampling pos)
                bool already = false;
                for (auto& s : v.selected)
                {
                    EigenVector3<T> dpos = samplePos[i];
                    double dx = s.pos.x() - dpos.x(), dy = s.pos.y() - dpos.y(),
                           dz = s.pos.z() - dpos.z();
                    if (dx * dx + dy * dy + dz * dz <= (spacing * 1e-10))
                    {
                        already = true;
                        break;
                    }
                }
                if (already) continue;
                if (maxAbsDiff[i] >= threshold)
                    discCands.emplace_back(maxAbsDiff[i], i);
            }
            sort(discCands.begin(), discCands.end(),
                 [](const std::pair<T, int>& a, const std::pair<T, int>& b)
                 { return a.first > b.first; });
            for (size_t di = 0;
                 di < discCands.size() && (int)v.selected.size() < cap; ++di)
            {
                int sidx = discCands[di].second;
                addSelected(samplePos[sidx], sampleVal[sidx],
                            SelectedPoint<T>::DISCONTINUITY);
            }
        }

        // Stage 3: fallback — if still nothing selected, project voxel center to iso or pick nearest sample to zero
        if (v.selected.empty())
        {
            EigenVector3<T> center(0.5 * (v.vmin.x() + v.vmax.x()),
                                   0.5 * (v.vmin.y() + v.vmax.y()),
                                   0.5 * (v.vmin.z() + v.vmax.z()));
            double sc = grid::value_at_2(grid, center);
            const double iso_eps = 1e-8;
            if (fabs(sc) <= iso_eps)
            {
                addSelected(center, sc, SelectedPoint<T>::CENTER_ISOSURFACE);
            }
            else
            {
                // estimate gradient via central differences (use small step related to voxel size)
                T h = std::max(1e-12,
                               std::min({spacing, (v.vmax.x() - v.vmin.x()),
                                         (v.vmax.y() - v.vmin.y()),
                                         (v.vmax.z() - v.vmin.z())})
                                   * 0.25);
                T vpx = grid::value_at_2(
                    grid,
                    EigenVector3<T>(center.x() + h, center.y(), center.z()));
                T vnx = grid::value_at_2(
                    grid,
                    EigenVector3<T>(center.x() - h, center.y(), center.z()));
                T vpy = grid::value_at_2(
                    grid,
                    EigenVector3<T>(center.x(), center.y() + h, center.z()));
                T vny = grid::value_at_2(
                    grid,
                    EigenVector3<T>(center.x(), center.y() - h, center.z()));
                T vpz = grid::value_at_2(
                    grid,
                    EigenVector3<T>(center.x(), center.y(), center.z() + h));
                T vnz = grid::value_at_2(
                    grid,
                    EigenVector3<T>(center.x(), center.y(), center.z() - h));
                T gx = (vpx - vnx) / (2.0 * h);
                T gy = (vpy - vny) / (2.0 * h);
                T gz = (vpz - vnz) / (2.0 * h);
                T gnorm = std::sqrt(gx * gx + gy * gy + gz * gz);
                if (gnorm > 1e-12)
                {
                    EigenVector3<T> n(gx / gnorm, gy / gnorm, gz / gnorm);
                    EigenVector3<T> p_iso(center.x() - sc * n.x(),
                                          center.y() - sc * n.y(),
                                          center.z() - sc * n.z());
                    T v_iso = grid::value_at_2(grid, p_iso);
                    addSelected(p_iso, v_iso,
                                SelectedPoint<T>::PROJECTED_ISOSURFACE);
                }
                else
                {
                    // fallback to pick the dense sample closest to zero
                    T bestAbs = std::numeric_limits<T>::infinity();
                    int besti = 0;
                    for (int i = 0; i < S3; ++i)
                    {
                        T a = std::fabs(sampleVal[i]);
                        if (a < bestAbs)
                        {
                            bestAbs = a;
                            besti = i;
                        }
                    }
                    addSelected(samplePos[besti], sampleVal[besti],
                                SelectedPoint<T>::NEAREST_GRID);
                }
            }
        }

        // enforce cap: if more than cap selected (unlikely), apply priority and trim
        if ((int)v.selected.size() > cap)
        {
            auto priority = [](const SelectedPoint<T>& sp) -> int
            {
                switch (sp.type)
                {
                case SelectedPoint<T>::PEAK:
                case SelectedPoint<T>::VALLEY:        return 3;
                case SelectedPoint<T>::DISCONTINUITY: return 2;
                default:                              return 1;
                }
            };
            std::sort(v.selected.begin(), v.selected.end(),
                      [&](const SelectedPoint<T>& a, const SelectedPoint<T>& b)
                      {
                          int pa = priority(a), pb = priority(b);
                          if (pa != pb) return pa > pb;
                          return fabs(a.sdf) < fabs(b.sdf);
                      });
            v.selected.resize(cap);
        }
    } // end voxels loop
}

} // namespace SDFSDFContact

#endif // GRID_SDF_SDF_CCD_GRADIENT_HPP
