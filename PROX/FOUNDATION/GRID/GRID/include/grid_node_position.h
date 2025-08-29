#ifndef GRID_NODE_POSITION_H
#define GRID_NODE_POSITION_H

#include <grid_grid.h>

namespace grid
{

    template<typename D, typename T>
    inline void node_position(
                     Grid<D,T> const & grid
                     , const Eigen::Matrix<size_t, 3, 1>& nodes
                     , Eigen::Matrix<T, 3, 1>& pos
                     )
    {
        pos.x() = nodes.x() * grid.dir().z() + grid.min().x();
        pos.y() = nodes.y() * grid.dir().y() + grid.min().y();
        pos.z() = nodes.z() * grid.dir().z() + grid.min().z();
    }

} // namespace grid

// GRID_NODE_POSITION_H
#endif
