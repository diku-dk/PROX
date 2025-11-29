#ifndef GRID_SDF_SDF_VOXELIZE_HPP
#define GRID_SDF_SDF_VOXELIZE_HPP

#include "grid_value_at.h"
#include <grid_grid.h>

namespace SDFSDFContact
{
template <typename T> struct SelectedPoint
{
    enum Type
    {
        PEAK,
        VALLEY,
        DISCONTINUITY,
        PROJECTED_ISOSURFACE,
        CENTER_ISOSURFACE,
        NEAREST_GRID
    } type;

    EigenVector3<T> pos;
    T sdf;
    // optional: store grid index if it originates from a grid sample
    int gx = -1, gy = -1, gz = -1;
};

template <typename T> struct SDFVoxel
{
    EigenVector3<T> vmin;
    EigenVector3<T> vmax;
    std::vector<EigenVector3<T>> vertices;
    uint32_t pointCount;
    uint32_t ixMin, ixMax;
    uint32_t iyMin, iyMax;
    uint32_t izMin, izMax;
    std::vector<SelectedPoint<T>> selected;
};

template <typename T>
std::vector<EigenVector3<T>> computeVertices(const EigenVector3<T>& vmin,
                                             const EigenVector3<T>& vmax)
{
    std::vector<EigenVector3<T>> corners;
    corners.push_back(
        EigenVector3<T>(vmin.x(), vmin.y(), vmin.z())); // (min, min, min)
    corners.push_back(
        EigenVector3<T>(vmax.x(), vmin.y(), vmin.z())); // (max, min, min)
    corners.push_back(
        EigenVector3<T>(vmin.x(), vmax.y(), vmin.z())); // (min, max, min)
    corners.push_back(
        EigenVector3<T>(vmax.x(), vmax.y(), vmin.z())); // (max, max, min)
    corners.push_back(
        EigenVector3<T>(vmin.x(), vmin.y(), vmax.z())); // (min, min, max)
    corners.push_back(
        EigenVector3<T>(vmax.x(), vmin.y(), vmax.z())); // (max, min, max)
    corners.push_back(
        EigenVector3<T>(vmin.x(), vmax.y(), vmax.z())); // (min, max, max)
    corners.push_back(
        EigenVector3<T>(vmax.x(), vmax.y(), vmax.z())); // (max, max, max)
    return corners;
}

template <typename T>
std::vector<SDFVoxel<T>>
createVoxels(const EigenVector3<T>& mmin, const EigenVector3<T>& mmax,
             uint32_t N, T point_spacing, uint32_t& out_total_grid_points_x,
             uint32_t& out_total_grid_points_y,
             uint32_t& out_total_grid_points_z)
{
    // Safety
    if (N <= 0) return {};

    // Compute voxel size in each axis
    const double eps = 1e-9;
    EigenVector3<T> fulls = mmax - mmin;
    EigenVector3<T> voxelSize{fulls.x() / T(N), fulls.y() / T(N),
                              fulls.z() / T(N)};

    // Compute how many sample points on each axis (discrete grid).
    // We assume points at positions: mmin + i * point_spacing for i = 0..(numPoints-1)
    auto computeNumPoints = [&](T len) -> uint32_t
    {
        if (len < -eps) return 0;
        // floor(len / spacing) gives the number of intervals; +1 gives number of points
        // add small epsilon to be robust to floating error
        return static_cast<uint32_t>(std::floor(len / point_spacing + 1e-9))
             + 1;
    };
    uint32_t nx = computeNumPoints(fulls.x());
    uint32_t ny = computeNumPoints(fulls.y());
    uint32_t nz = computeNumPoints(fulls.z());

    out_total_grid_points_x = nx;
    out_total_grid_points_y = ny;
    out_total_grid_points_z = nz;

    std::vector<SDFVoxel<T>> voxels;
    voxels.reserve(N * N * N);

    for (uint32_t iz = 0; iz < N; ++iz)
    {
        for (uint32_t iy = 0; iy < N; ++iy)
        {
            for (uint32_t ix = 0; ix < N; ++ix)
            {
                EigenVector3<T> vmin_vox{mmin.x() + ix * voxelSize.x(),
                                         mmin.y() + iy * voxelSize.y(),
                                         mmin.z() + iz * voxelSize.z()};
                EigenVector3<T> vmax_vox{mmin.x() + (ix + 1) * voxelSize.x(),
                                         mmin.y() + (iy + 1) * voxelSize.y(),
                                         mmin.z() + (iz + 1) * voxelSize.z()};

                // Convert voxel bounds to discrete grid index ranges [i_min, i_max]
                auto idx_min = [&](T coord_min_axis, T axis_mmin) -> uint32_t
                {
                    T t = (coord_min_axis - axis_mmin) / point_spacing;
                    return static_cast<uint32_t>(
                        ceil(t - 1e-9)); // inclusive lower index
                };
                auto idx_max = [&](T coord_max_axis, T axis_mmin) -> uint32_t
                {
                    T t = (coord_max_axis - axis_mmin) / point_spacing;
                    return static_cast<uint32_t>(
                        floor(t + 1e-9)); // inclusive upper index
                };

                uint32_t ix_min = idx_min(vmin_vox.x(), mmin.x());
                uint32_t iy_min = idx_min(vmin_vox.y(), mmin.y());
                uint32_t iz_min = idx_min(vmin_vox.z(), mmin.z());

                uint32_t ix_max = idx_max(vmax_vox.x(), mmin.x());
                uint32_t iy_max = idx_max(vmax_vox.y(), mmin.y());
                uint32_t iz_max = idx_max(vmax_vox.z(), mmin.z());

                // Clamp to valid global grid indices
                ix_min
                    = std::max<uint32_t>(0, std::min<uint32_t>(ix_min, nx - 1));
                iy_min
                    = std::max<uint32_t>(0, std::min<uint32_t>(iy_min, ny - 1));
                iz_min
                    = std::max<uint32_t>(0, std::min<uint32_t>(iz_min, nz - 1));

                ix_max
                    = std::max<uint32_t>(0, std::min<uint32_t>(ix_max, nx - 1));
                iy_max
                    = std::max<uint32_t>(0, std::min<uint32_t>(iy_max, ny - 1));
                iz_max
                    = std::max<uint32_t>(0, std::min<uint32_t>(iz_max, nz - 1));

                // If min > max after clamping, there are zero points on that axis
                uint32_t cx = (ix_min <= ix_max) ? (ix_max - ix_min + 1) : 0;
                uint32_t cy = (iy_min <= iy_max) ? (iy_max - iy_min + 1) : 0;
                uint32_t cz = (iz_min <= iz_max) ? (iz_max - iz_min + 1) : 0;

                uint32_t count = cx * cy * cz;

                SDFVoxel<T> v;
                v.vmin = vmin_vox;
                v.vmax = vmax_vox;
                v.vertices = computeVertices(v.vmin, v.vmax);
                v.pointCount = count;
                v.ixMin = ix_min;
                v.ixMax = ix_max;
                v.iyMin = iy_min;
                v.iyMax = iy_max;
                v.izMin = iz_min;
                v.izMax = iz_max;

                voxels.push_back(v);
            }
        }
    }
    return voxels;
}

template <typename D, typename T>
bool voxelContainsZeroSample(const SDFVoxel<T>& voxels,
                             const EigenVector3<T>& global_mmin,
                             T point_spacing, grid::Grid<D, T> gridSdf,
                             T zero_eps = 1e-9)
{
    if (voxels.pointCount == 0) return false; // no points to check

    //If there exists both a negative and a positive point, we know for a fact
    // the isosurface is contained somewhat within our voxel!
    bool negativeExists = false;
    bool positiveExists = false;
    // iterate explicitly through each point index in the 3D subgrid
    for (int iz = voxels.izMin; iz <= voxels.izMax; ++iz)
    {
        for (int iy = voxels.iyMin; iy <= voxels.iyMax; ++iy)
        {
            for (int ix = voxels.ixMin; ix <= voxels.ixMax; ++ix)
            {
                // compute world coordinate for this grid index
                EigenVector3<T> p{global_mmin.x() + ix * point_spacing,
                                  global_mmin.y() + iy * point_spacing,
                                  global_mmin.z() + iz * point_spacing};

                T val = grid::value_at_2(gridSdf, p);
                if (std::isnan(val)) continue;
                if (val > 0) positiveExists = true;
                if (val < 0) negativeExists = true;
                if (std::fabs(val) <= zero_eps)
                {
                    return true;
                }
            }
        }
    }
    return (negativeExists && positiveExists);
}

template <typename D, typename T>
std::vector<SDFVoxel<T>>
voxelFilterAll(const std::vector<SDFVoxel<T>>& voxels,
               const EigenVector3<T>& global_mmin, T point_spacing,
               grid::Grid<D, T> gridSdf, T zero_eps = 1e-9)
{
    std::vector<SDFVoxel<T>> returnVoxels;
    for (size_t i = 0; i < voxels.size(); ++i)
    {
        SDFVoxel<T> singleVoxel = voxels[i];

        bool shouldKeep = voxelContainsZeroSample(
            singleVoxel, global_mmin, point_spacing, gridSdf, zero_eps);
        if (shouldKeep) { returnVoxels.push_back(singleVoxel); }
    }
    return returnVoxels;
}
} // namespace SDFSDFContact

#endif // GRID_SDF_SDF_VOXELIZE_HPP
