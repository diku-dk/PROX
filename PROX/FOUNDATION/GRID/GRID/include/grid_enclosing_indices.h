#ifndef GRID_ENCLOSING_INDICES_H
#define GRID_ENCLOSING_INDICES_H

#include <grid_grid.h>


#include <cmath>

namespace grid
{
  template <typename D, typename T>
  inline void enclosing_indices(
                                Grid<D,T> const & grid
                                , const Eigen::Matrix<T, 3, 1>& pos
                                , Eigen::Matrix<size_t, 3, 1>& nodes0
                                , Eigen::Matrix<size_t, 3, 1>& nodes1
                                )
  {
    using std::floor;


    /*T const diff_x = ( pos.x - grid.min().x ) / grid.dx();
    T const diff_y = ( pos.y - grid.min().x ) / grid.dy();
    T const diff_z = ( pos.z - grid.min().x ) / grid.dz();*/
    const Eigen::Matrix<T, 3, 1> diff = (pos - Eigen::Matrix<T, 3, 1>(grid.min().x))/grid.dir();

    nodes0.x() = static_cast<size_t>( floor( diff.x ) );
    nodes0.y() = static_cast<size_t>( floor( diff.y ) );
    nodes0.z() = static_cast<size_t>( floor( diff.z ) );

    nodes1.x() = ( nodes0.x() + 1 );
    nodes1.y() = ( nodes0.y() + 1 );
    nodes1.z() = ( nodes0.z() + 1 );

  }

} // namespace grid

// GRID_ENCLOSING_INDICES_H
#endif
