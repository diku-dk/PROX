#ifndef GRID_NODE_POSITION_H
#define GRID_NODE_POSITION_H

#include <grid_grid.h>

namespace grid
{

  template<typename D, typename T>
  inline
  void node_position(
                     Grid<D,T> const & grid
                     , size_t const & i
                     , size_t const & j
                     , size_t const & k
                     , T & x
                     , T & y
                     , T & z
                     )
  {
    x = i * grid.dx() + grid.min_x();
    y = j * grid.dy() + grid.min_y();
    z = k * grid.dz() + grid.min_z();
  }

} // namespace grid

// GRID_NODE_POSITION_H
#endif
