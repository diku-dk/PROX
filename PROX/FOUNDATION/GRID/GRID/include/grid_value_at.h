#ifndef GRID_VALUE_AT_H
#define GRID_VALUE_AT_H

#include <grid_enclosing_indices.h>

#include <tiny_is_number.h>
#include <tiny_is_finite.h>
#include <boost/numeric/conversion/bounds.hpp>

#include <cmath>

template <typename T,typename T2>
inline T  trilinear(
    const T d000
    , const T d001
    , const T d010
    , const T d011
    , const T d100
    , const T d101
    , const T d110
    , const T d111
    , const Eigen::Matrix<T, 3, 1>& stu
    )
{

    T const x00 = boost::numeric_cast<T>( (d001 - d000)*stu.x() + d000 );
    T const x01 = boost::numeric_cast<T>( (d011 - d010)*stu.x() + d010 );
    T const x10 = boost::numeric_cast<T>( (d101 - d100)*stu.x() + d100 );
    T const x11 = boost::numeric_cast<T>( (d111 - d110)*stu.x() + d110 );
    T const y0  = boost::numeric_cast<T>( ( x01 -  x00)*stu.y() +  x00 );
    T const y1  = boost::numeric_cast<T>( ( x11 -  x10)*stu.y() +  x10 );
    return boost::numeric_cast<T>( (y1-y0)*stu.z() + y0 );
}

namespace grid
{

  template<typename D,typename T>
  inline typename D value_at(Grid<D,T> const & grid, const Eigen::Matrix<T, 3, 1>& pos)
  {
      assert( tiny::is_number(pos.x()) || !"value_at(): x was not a number");
      assert( tiny::is_finite(pos.x()) || !"value_at(): x was not finite"  );
      assert( tiny::is_number(pos.y()) || !"value_at(): y was not a number");
      assert( tiny::is_finite(pos.y()) || !"value_at(): y was not finite"  );
      assert( tiny::is_number(pos.z()) || !"value_at(): z was not a number");
      assert( tiny::is_finite(pos.z()) || !"value_at(): z was not finite"  );

    using std::min;
    using std::max;

    Eigen::Matrix<T, 3, 1> safe = max(grid.min(), min(grid.max()));
    /*T const safe_x = max( grid.min_x(), min( x, grid.max_x()) );
    T const safe_y = max( grid.min_y(), min( y, grid.max_y()) );
    T const safe_z = max( grid.min_z(), min( z, grid.max_z()) );*/

    //size_t i0, j0, k0, i1, j1, k1;
    //enclosing_indices( grid, x, y, z, i0, j0, k0, i1, j1, k1 );

    Eigen::Matrix<size_t, 3, 1> nodes0;
    Eigen::Matrix<size_t, 3, 1> nodes1;
    enclosing_indices( grid, pos, nodes0, nodes1);


/*    D const d000 = grid( i0, j0, k0 );
    D const d001 = grid( i1, j0, k0 );
    D const d010 = grid( i0, j1, k0 );
    D const d011 = grid( i1, j1, k0 );
    D const d100 = grid( i0, j0, k1 );
    D const d101 = grid( i1, j0, k1 );
    D const d110 = grid( i0, j1, k1 );
    D const d111 = grid( i1, j1, k1 );*/
    D const d000 = grid( nodes0);
    D const d001 = grid( Eigen::Matrix<size_t, 3, 1>(nodes1.x(), nodes0.y(), nodes0.z()) );
    D const d010 = grid( Eigen::Matrix<size_t, 3, 1>(nodes0.x(), nodes1.y(), nodes0.z()) );
    D const d011 = grid( Eigen::Matrix<size_t, 3, 1>(nodes1.x(), nodes1.y(), nodes0.z()) );
    D const d100 = grid( Eigen::Matrix<size_t, 3, 1>(nodes0.x(), nodes0.y(), nodes1.z()) );
    D const d101 = grid( Eigen::Matrix<size_t, 3, 1>(nodes1.x(), nodes0.y(), nodes1.z()) );
    D const d110 = grid( Eigen::Matrix<size_t, 3, 1>(nodes0.x(), nodes1.y(), nodes1.z()) );
    D const d111 = grid( Eigen::Matrix<size_t, 3, 1>(nodes1.x(), nodes1.y(), nodes1.z()) );

    /*T const s = ( x - ( i0*grid.dx() + grid.min_x() ) ) / grid.dx();
    T const t = ( y - ( j0*grid.dy() + grid.min_y() ) ) / grid.dy();
    T const u = ( z - ( k0*grid.dz() + grid.min_z() ) ) / grid.dz();*/
    const Eigen::Matrix<T, 3, 1> stu = (pos-(nodes0*grid.dir() + grid.min()))/grid.dir();

    return trillinear<D, T>( d000, d001, d010, d011, d100, d101, d110, d111, stu );
  }

} // namespace grid

// GRID_VALUE_AT_H
#endif
