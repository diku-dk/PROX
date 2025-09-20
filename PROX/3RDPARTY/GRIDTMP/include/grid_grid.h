#ifndef GRID_GRID_H
#define GRID_GRID_H

//#include <vector>
#include <iostream>

#include "eigen3/Eigen/Dense"

namespace grid
{

template <typename T> struct GridTriangle
{
    Eigen::Matrix<T, 3, 1> v0;
    Eigen::Matrix<T, 3, 1> v1;
    Eigen::Matrix<T, 3, 1> v2;
};
template < typename D, typename T>
class Grid
{
public:
    // Make core members accessible to existing iterator code that expects them.
    Eigen::Matrix<T, 3, 1> m_min;
    Eigen::Matrix<T, 3, 1> m_max;
    Eigen::Matrix<T, 3, 1> m_dir;
    Eigen::Matrix<size_t, 3, 1> m_nodes;
    std::vector<GridTriangle<T>> m_temporaryGridStructure;

    //vvv The data values stored at the grid nodes.
    //Internal data stored as Eigen::Matrix (rows = total nodes, cols = 1)
    Eigen::Matrix<D, Eigen::Dynamic, Eigen::Dynamic> m_data;

public:
    using value_type = D;
    Grid()
        : m_min(T(0),T(0),T(0))
        , m_max(T(0),T(0),T(0))
        , m_dir(T(0),T(0),T(0))
        , m_nodes(0u,0u,0u)
        , m_data()
    {}

    Grid(Grid<D,T> const & G) = default;
    Grid<D,T> & operator=(Grid<D,T> const & grid) = default;

protected:
    inline size_t index(const Eigen::Matrix<size_t, 3, 1>& nodes) const
    {
        assert( nodes.x() < this->m_nodes.x() && nodes.y() < this->m_nodes.y() && nodes.z() < this->m_nodes.z() || !"index(): i, j, or k was out of bounds");
        return (nodes.z()*this->m_nodes.y() + nodes.y())*this->m_nodes.x() + nodes.x();
    }

    inline size_t index(size_t i, size_t j, size_t k) const
    {
        assert( i < this->m_nodes.x() && j < this->m_nodes.y() && k < this->m_nodes.z() || !"index(): i, j, or k was out of bounds");
        return (k*this->m_nodes.y() + j)*this->m_nodes.x() + i;
    }

public:
    void create(const Eigen::Matrix<T, 3, 1>& min, const Eigen::Matrix<T, 3, 1>& max,
                const Eigen::Matrix<size_t, 3, 1>& nodes)
    {
        assert( (max.array() > min.array()).all() || !"create(): max was less than or equal min");
        assert( (nodes.array() > 0u).all() || !"create(): all nodes must be positive");

        this->m_min = min;
        this->m_max = max;

        this->m_dir.x() = (max.x()-min.x())/(nodes.x()-1);
        this->m_dir.y() = (max.y()-min.y())/(nodes.y()-1);
        this->m_dir.z() = (max.z()-min.z())/(nodes.z()-1);

        this->m_nodes = nodes;

        const size_t total = nodes.x()*nodes.y()*nodes.z();
        this->m_data.resize((Eigen::Index)total, 1);

        for (Eigen::Index r = 0; r < this->m_data.rows(); ++r) this->m_data(r,0) = D();
    }

    // Access by node index-vector
    D & operator()(const Eigen::Matrix<size_t, 3, 1>& nodes)
    {
        assert(this->m_data.size()>0 || !"operator(): no data");
        return this->m_data( (Eigen::Index)index(nodes), 0 );
    }

    D const & operator() (const Eigen::Matrix<size_t, 3, 1>& nodes) const
    {
        assert(this->m_data.size()>0 || !"operator(): no data");
        return this->m_data( (Eigen::Index)index(nodes), 0 );
    }


    Eigen::Matrix<T, 3, 1> dims()  const { return this->m_max - this->m_min; }
    size_t size()  const { return (size_t)this->m_data.rows()* (size_t)std::max<Eigen::Index>(1,this->m_data.cols()); }
    bool empty() const { return this->m_data.size()==0; }

    const Eigen::Matrix<T, 3, 1>& min() const { return this->m_min; }
    const Eigen::Matrix<T, 3, 1>& max() const { return this->m_max; }
    const Eigen::Matrix<T, 3, 1>& dir() const { return this->m_dir; }
    const Eigen::Matrix<size_t, 3, 1>& nodes() const { return this->m_nodes; }

    //For backwards compatib
    size_t I() const { return this->m_nodes.x(); }
    size_t J() const { return this->m_nodes.y(); }
    size_t K() const { return this->m_nodes.z(); }

    //Hopefully this will work with eigen matrix!
    D       * data_ptr()       { return this->m_data.data(); }
    D const * data_ptr() const { return this->m_data.data(); }

    //If we wish we can now do data() -- useful, i dont know?
    D* data() {return this->m_data.data(); }
    const D* data() const { return this->m_data.data(); }
};
} // namespace grid

// GRID_GRID_H
#endif
