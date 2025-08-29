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

    const Eigen::Matrix<T, 3, 1> diff = (pos - grid.min()).cwiseQuotient(grid.dir());

    nodes0.x() = static_cast<size_t>( floor( diff.x() ) );
    nodes0.y() = static_cast<size_t>( floor( diff.y() ) );
    nodes0.z() = static_cast<size_t>( floor( diff.z() ) );

    nodes1.x() = ( nodes0.x() + 1 );
    nodes1.y() = ( nodes0.y() + 1 );
    nodes1.z() = ( nodes0.z() + 1 );
}
}

// GRID_ENCLOSING_INDICES_H
#endif
