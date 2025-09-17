#ifndef GRID_VALUE_AT_H
#define GRID_VALUE_AT_H

#include <grid_enclosing_indices.h>

#include <boost/numeric/conversion/bounds.hpp>
#include <boost/cast.hpp>

#include <cmath>

namespace grid
{
    template<typename D,typename T>
    inline D trilinear(
        const D d000
        , const D d001
        , const D d010
        , const D d011
        , const D d100
        , const D d101
        , const D d110
        , const D d111
        , const Eigen::Matrix<T, 3, 1>& stu
        )
    {
        D const x00 = static_cast<D>( (d001 - d000)*stu.x() + d000 );
        D const x01 = static_cast<D>( (d011 - d010)*stu.x() + d010 );
        D const x10 = static_cast<D>( (d101 - d100)*stu.x() + d100 );
        D const x11 = static_cast<D>( (d111 - d110)*stu.x() + d110 );
        D const y0  = static_cast<D>( ( x01 -  x00)*stu.y() +  x00 );
        D const y1  = static_cast<D>( ( x11 -  x10)*stu.y() +  x10 );
        return static_cast<D>( (y1-y0)*stu.z() + y0 );
    }

    /*template<typename D,typename T>
    inline typename std::enable_if<!std::is_same<D,T>::value, D>::type value_at(Grid<D,T> const & grid, const Eigen::Matrix<T, 3, 1>& pos)
    {
        (void)0;
        return value_at<D,T>(grid,pos);
    }*/

    template<typename D,typename T>
    inline D value_at(Grid<D,T> const & grid, const Eigen::Matrix<T, 3, 1>& pos)
    {
        Eigen::Matrix<size_t, 3, 1> nodes0;
        Eigen::Matrix<size_t, 3, 1> nodes1;
        enclosing_indices( grid, pos, nodes0, nodes1);

        nodes0 = Eigen::Matrix<size_t, 3, 1>(
            std::max<size_t>(std::min<T>(nodes0.x(), grid.I() - 1), size_t(0)),
            std::max<size_t>(std::min<T>(nodes0.y(), grid.J() - 1), size_t(0)),
            std::max<size_t>(std::min<T>(nodes0.z(), grid.K() - 1), size_t(0)));
        nodes1 = Eigen::Matrix<size_t, 3, 1>(
            std::max<size_t>(std::min<size_t>(nodes1.x(), grid.I() - 1),
                             size_t(0)),
            std::max<size_t>(std::min<size_t>(nodes1.y(), grid.J() - 1),
                             size_t(0)),
            std::max<size_t>(std::min<size_t>(nodes1.z(), grid.K() - 1),
                             size_t(0)));
        D const d000 = grid(nodes0);
        D const d001 = grid( Eigen::Matrix<size_t, 3, 1>(nodes1.x(), nodes0.y(), nodes0.z()) );
        D const d010 = grid( Eigen::Matrix<size_t, 3, 1>(nodes0.x(), nodes1.y(), nodes0.z()) );
        D const d011 = grid( Eigen::Matrix<size_t, 3, 1>(nodes1.x(), nodes1.y(), nodes0.z()) );
        D const d100 = grid( Eigen::Matrix<size_t, 3, 1>(nodes0.x(), nodes0.y(), nodes1.z()) );
        D const d101 = grid( Eigen::Matrix<size_t, 3, 1>(nodes1.x(), nodes0.y(), nodes1.z()) );
        D const d110 = grid( Eigen::Matrix<size_t, 3, 1>(nodes0.x(), nodes1.y(), nodes1.z()) );
        D const d111 = grid( Eigen::Matrix<size_t, 3, 1>(nodes1.x(), nodes1.y(), nodes1.z()) );

        const Eigen::Matrix<T, 3, 1> stu = (pos - (nodes0.template cast<T>().cwiseProduct(grid.dir()) + grid.min())) .cwiseQuotient(grid.dir());

        return trilinear<D,T>( d000, d001, d010, d011, d100, d101, d110, d111, stu );
    }

}

// GRID_VALUE_AT_H
#endif
