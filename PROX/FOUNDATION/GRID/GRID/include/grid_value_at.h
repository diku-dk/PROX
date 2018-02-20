#ifndef GRID_VALUE_AT_H
#define GRID_VALUE_AT_H

#include <grid_enclosing_indices.h>

#include <tiny_trillinear.h>
#include <tiny_is_number.h>
#include <tiny_is_finite.h>

#include <cmath>

namespace grid
{

  template<typename D,typename T>
  inline typename D value_at(Grid<D,T> const & grid, T const & x,  T const & y, T const & z )
  {
    assert( is_number(x) || !"value_at(): x was not a number");
    assert( is_finite(x) || !"value_at(): x was not finite"  );
    assert( is_number(y) || !"value_at(): y was not a number");
    assert( is_finite(y) || !"value_at(): y was not finite"  );
    assert( is_number(z) || !"value_at(): z was not a number");
    assert( is_finite(z) || !"value_at(): z was not finite"  );

    using std::min;
    using std::max;

    T const safe_x = max( grid.min_x(), min( x, grid.max_x()) );
    T const safe_y = max( grid.min_y(), min( y, grid.max_y()) );
    T const safe_z = max( grid.min_z(), min( z, grid.max_z()) );

    size_t i0, j0, k0, i1, j1, k1;

    enclosing_indices( grid, x, y, z, i0, j0, k0, i1, j1, k1 );

    D const d000 = grid( i0, j0, k0 );
    D const d001 = grid( i1, j0, k0 );
    D const d010 = grid( i0, j1, k0 );
    D const d011 = grid( i1, j1, k0 );
    D const d100 = grid( i0, j0, k1 );
    D const d101 = grid( i1, j0, k1 );
    D const d110 = grid( i0, j1, k1 );
    D const d111 = grid( i1, j1, k1 );

    T const s = ( x - ( i0*grid.dx() + grid.min_x() ) ) / grid.dx();
    T const t = ( y - ( j0*grid.dy() + grid.min_y() ) ) / grid.dy();
    T const u = ( z - ( k0*grid.dz() + grid.min_z() ) ) / grid.dz();

    return tiny::trillinear( d000, d001, d010, d011, d100, d101, d110, d111, s, t, u );
  }

} // namespace grid

// GRID_VALUE_AT_H
#endif
