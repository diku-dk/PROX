#ifndef BROAD_CCD_BUILD_HPP
#define BROAD_CCD_BUILD_HPP

#include <iostream>
#include <eigenhelperall.h>

namespace broad_ccd
{
template <typename T> class AABB
{
public:
    EigenVector3<T> min;
    EigenVector3<T> max;

    AABB() = default;

    AABB(const EigenVector3<T>& min, const EigenVector3<T>& max)
        : min(min)
        , max(max)
    {
    }

    // Expand AABB to include another AABB
    void expand(const AABB& other)
    {
        min.x() = std::min(min.x(), other.min.x());
        min.y() = std::min(min.y(), other.min.y());
        min.z() = std::min(min.z(), other.min.z());
        max.x() = std::max(max.x(), other.max.x());
        max.y() = std::max(max.y(), other.max.y());
        max.z() = std::max(max.z(), other.max.z());
    }

    // Expand AABB by a margin in all directions
    void expandByMargin(float margin)
    {
        min = min - EigenVector3<T>(margin, margin, margin);
        max = max + EigenVector3<T>(margin, margin, margin);
    }

    // Check if this AABB overlaps with another
    bool overlaps(const AABB& other) const
    {
        return (min.x() <= other.max.x() && max.x() >= other.min.x())
            && (min.y() <= other.max.y() && max.y() >= other.min.y())
            && (min.z() <= other.max.z() && max.z() >= other.min.z());
    }

    // Get center of AABB
    EigenVector3<T> getCenter() const
    {
        return {(min.x() + max.x()) * 0.5f, (min.y() + max.y()) * 0.5f,
                (min.z() + max.z()) * 0.5f};
    }

    // Get diagonal length
    float getDiagonal() const { return (max - min).length(); }
};

template <typename T> struct RigidBody
{
    EigenVector3<T> position;
    EigenVector3<T> linearVelocity;
    EigenVector3<T> angularVelocity;

    // Precomputed properties
    AABB<T> localAABB;  // AABB in local space (before transform)
    T maxDistanceFromCenter; // Maximum distance from center to any vertex

    // Get conservative AABB for time interval [t0, t1]
    AABB<T> getConservativeAABB(T t0, T t1) const
    {
        T dt = t1 - t0;

        // Start with local AABB translated to current position
        AABB<T> currentAABB(localAABB.min + position, localAABB.max + position);

        // Calculate expansion due to linear and angular motion
        T linearExpansion = linearVelocity.norm() * dt;
        T angularExpansion
            = maxDistanceFromCenter * angularVelocity.norm() * dt;
        T totalExpansion = linearExpansion + angularExpansion;

        // Expand AABB conservatively
        AABB<T> expandedAABB = currentAABB;
        expandedAABB.expandByMargin(totalExpansion);

        return expandedAABB;
    }
};
} // namespace broad_ccd

namespace broad_ccd
{
template <typename T> struct BVHNode
{
    AABB<T> bounds;
    int bodyIndex;  // -1 for internal nodes
    BVHNode<T>* left;
    BVHNode<T>* right;

    BVHNode()
        : bodyIndex(-1)
        , left(nullptr)
        , right(nullptr)
    {
    }
};

template <typename T> class BVH
{
public:
    BVH(std::vector<RigidBody<T>>& bodies, float t0, float t1)
        : bodies(bodies)
        , t0(t0)
        , t1(t1)
        , root(nullptr)
    {
    }

    ~BVH() { deleteRecursive(root); }

    void build()
    {
        std::vector<int> indices(bodies.size());
        for (int i = 0; i < bodies.size(); i++) { indices[i] = i; }
        root = buildRecursive(indices, 0, indices.size());
    }

    // Find all bodies that might collide with the given body
    std::vector<int> query(int bodyIndex)
    {
        std::vector<int> results;
        queryRecursive(root, bodyIndex, results);
        return results;
    }

    // Get all potential collision pairs
    std::vector<std::pair<int, int>> getAllPairs()
    {
        std::vector<std::pair<int, int>> pairs;

        for (int i = 0; i < bodies.size(); i++)
        {
            std::vector<int> candidates = query(i);
            for (int j : candidates)
            {
                if (i < j)
                {  // Avoid duplicates
                    pairs.emplace_back(i, j);
                }
            }
        }

        return pairs;
    }

private:
    BVHNode<T>* root;
    std::vector<RigidBody<T>>& bodies;
    float t0, t1;  // Current time interval

    BVHNode<T>* buildRecursive(std::vector<int>& indices, int start, int end)
    {
        if (start >= end) return nullptr;

        BVHNode<T>* node = new BVHNode<T>();

        if (end - start == 1)
        {
            // Leaf node
            node->bodyIndex = indices[start];
            node->bounds = bodies[indices[start]].getConservativeAABB(t0, t1);
            return node;
        }

        // Calculate bounds for all objects in this node
        AABB<T> centroidBounds;
        for (int i = start; i < end; i++)
        {
            AABB<T> bodyBounds = bodies[indices[i]].getConservativeAABB(t0, t1);
            centroidBounds.expand(bodyBounds);
        }

        // Choose split axis (longest axis)
        EigenVector3<T> diagonal = centroidBounds.max - centroidBounds.min;
        int axis = 0;
        if (diagonal.y() > diagonal.x()) axis = 1;
        if (diagonal.z() > diagonal[axis]) axis = 2;

        // Sort indices along chosen axis
        std::sort(indices.begin() + start, indices.begin() + end,
                  [&](int a, int b)
                  {
                      EigenVector3<T> centerA
                          = bodies[a].getConservativeAABB(t0, t1).getCenter();
                      EigenVector3<T> centerB
                          = bodies[b].getConservativeAABB(t0, t1).getCenter();
                      return centerA[axis] < centerB[axis];
                  });

        // Split at median
        int mid = start + (end - start) / 2;
        node->left = buildRecursive(indices, start, mid);
        node->right = buildRecursive(indices, mid, end);

        // Combine children bounds
        node->bounds = node->left->bounds;
        node->bounds.expand(node->right->bounds);

        return node;
    }

    void queryRecursive(BVHNode<T>* node, int queryIndex,
                        std::vector<int>& results)
    {
        if (!node || node->bodyIndex == queryIndex) return;

        if (!bodies[queryIndex].getConservativeAABB(t0, t1).overlaps(
                node->bounds))
        {
            return;
        }

        if (node->bodyIndex != -1)
        {
            // Leaf node
            results.push_back(node->bodyIndex);
        }
        else
        {
            queryRecursive(node->left, queryIndex, results);
            queryRecursive(node->right, queryIndex, results);
        }
    }

    void deleteRecursive(BVHNode<T>* node)
    {
        if (!node) return;
        deleteRecursive(node->left);
        deleteRecursive(node->right);
        delete node;
    }
};
} // namespace broad_ccd

#endif // BROAD_CCD_BUILD_HPP
