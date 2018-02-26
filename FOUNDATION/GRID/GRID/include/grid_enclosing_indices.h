#ifndef GRID_ENCLOSING_INDICES_H
#define GRID_ENCLOSING_INDICES_H

#include <grid_grid.h>

#include <cmath>

namespace grid
{
  template <typename D, typename T>
  inline void enclosing_indices(
                                Grid<D,T> const & grid
                                , T const & x
                                , T const & y
                                , T const & z
                                , size_t & i0
                                , size_t & j0
                                , size_t & k0
                                , size_t & i1
                                , size_t & j1
                                , size_t & k1
                                )
  {
    using std::floor;

    T const diff_x = ( x - grid.min_x() ) / grid.dx();
    T const diff_y = ( y - grid.min_x() ) / grid.dy();
    T const diff_z = ( z - grid.min_x() ) / grid.dz();

    i0 = static_cast<size_t>( floor( diff_x ) );
    j0 = static_cast<size_t>( floor( diff_y ) );
    k0 = static_cast<size_t>( floor( diff_z ) );

    i1 = ( i0 + 1 );
    j1 = ( j0 + 1 );
    k1 = ( k0 + 1 );

  }

} // namespace grid

// GRID_ENCLOSING_INDICES_H
#endif
