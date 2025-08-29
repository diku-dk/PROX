#ifndef GRID_GRID_H
#define GRID_GRID_H

#include <tiny_value_traits.h>

//#include <vector>
#include <iostream>

#include "eigen3/Eigen/Dense"



namespace grid
{


  template < typename D, typename T>
  class Grid
  {
  public:

    //typedef  tiny::ValueTraits<T>  VT;

  protected:

    Eigen::Matrix<T, 3, 1> m_min;
    Eigen::Matrix<T, 3, 1> m_max;
    Eigen::Matrix<T, 3, 1> m_dir;
    Eigen::Matrix<size_t, 3, 1> m_nodes;
    //T      m_min_x;      ///< Minimum x coordinate (x coordinate of node (0,0,0)
    //T      m_min_y;      ///< Minimum y coordinate (y coordinate of node (0,0,0)
    //T      m_min_z;      ///< Minimum z coordinate (z coordinate of node (0,0,0)
    //T      m_max_x;      ///< Maximum x coordinate (x coordinate of node (I-1,J-1,K-1)
    //T      m_max_y;      ///< Maximum y coordinate (y coordinate of node (I-1,J-1,K-1)
    //T      m_max_z;      ///< Maximum z coordinate (z coordinate of node (I-1,J-1,K-1)
    //T      m_dx;         ///< Space between nodes in x direction.
    //T      m_dy;         ///< Space between nodes in y direction.
    //T      m_dz;         ///< Space between nodes in z direction.
    //size_t m_I;          ///< Number of nodes along x-axis.
    //size_t m_J;          ///< Number of nodes along y-axis.
    //size_t m_K;          ///< Number of nodes along z-axis.

    std::vector<D> m_data; ///< The data values stored at the grid nodes.

  protected:

    size_t index(const Eigen::Matrix<size_t, 3, 1>& nodes) const
    {
      assert( nodes < this->m_nodes|| !"index(): i, j, or k was out of bounds");

        return (nodes.z()*this->m_nodes.y() + nodes.y())*this->m_nodes.x() + nodes.x();
    }

  public:

    Grid()
    : m_min(T(0.0), T(0.0), T(0.0))
    , m_max(T(0.0), T(0.0), T(0.0))
    , m_dir(T(0.0), T(0.0), T(0.0))
    , m_nodes(0u, 0u, 0u)
    , m_data( )
    {}

    Grid(Grid<D,T,I> const & G)
    {
      *this = G;
    }

    ~Grid(){}

    Grid<D,T,I> & operator=(Grid<D,T,I> const & grid)
    {
      if( this != &grid)
      {
        this->m_min = grid.m_min;
        this->m_max = grid.m_max;
        this->m_dir = grid.m_dir;
        this->m_nodes = grid.m_nodes;
        this->m_data  = grid.m_data;
      }
      return (*this);
    }

  public:

    void create(const Eigen::Matrix<T, 3, 1>& min, const Eigen::Matrix<T, 3, 1>& max,
                const Eigen::Matrix<size_t, 3, 1>& nodes
                )
    {
      assert( max > min || !"create(): max was less than or equal min");

      assert( nodes>0u           || !"create(): all nodes must be positive");

      this->m_min = min;

      this->m_max = max;

      this->m_dir.x = (max.x()-min.x())/(nodes.x()-1);
      this->m_dir.y = (max.y()-min.y())/(nodes.y()-1);
      this->m_dir.z = (max.z()-min.z())/(nodes.z()-1);


      this->m_nodes = nodes;

      m_data.resize( nodes.x()*nodes.y()*nodes.z() );
    }

    D & operator()(const Eigen::Matrix<size_t, 3, 1>& nodes)
    {
      assert(this->m_data.size()>0 || !"operator(): no data");

        return this->m_data[ this->index(nodes.x(),nodes.y(),nodes.z()) ];
    }

    D const & operator() (const Eigen::Matrix<size_t, 3, 1>& nodes) const
    {
      assert(this->m_data.size()>0 || !"operator(): no data");

      return this->m_data[ this->index(nodes.x(),nodes.y(),nodes.z()) ];
    }

    Eigen::Matrix<T, 3, 1> dims()  const { return this->m_max - this->m_min; }
    /*T height() const { return this->m_max_y - this->m_min_y; }
    T depth()  const { return this->m_max_z - this->m_min_z; }*/

    size_t size()  const { return this->m_data.size(); }
    bool   empty() const { return this->m_data.empty(); }

    const Eigen::Matrix<T, 3, 1>& min() const { return this->m_min; }
    const Eigen::Matrix<T, 3, 1>& max() const { return this->m_max; }

    const Eigen::Matrix<T, 3, 1>& dir() const { return this->m_dir; }

    const Eigen::Matrix<size_t, 3, 1>& nodes() const { return this->m_nodes; }


    D       * data_ptr()       { return this->m_data[0]; }
    D const * data_ptr() const { return this->m_data[0]; }

  };

} // namespace grid

// GRID_GRID_H
#endif
