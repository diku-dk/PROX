#ifndef GRID_GRID_H
#define GRID_GRID_H

#include <tiny_value_traits.h>

#include <vector>

namespace grid
{


  template < typename D, typename T>
  class Grid
  {
  public:

    typedef  tiny::ValueTraits<T>  VT;

  protected:

    T      m_min_x;      ///< Minimum x coordinate (x coordinate of node (0,0,0)
    T      m_min_y;      ///< Minimum y coordinate (y coordinate of node (0,0,0)
    T      m_min_z;      ///< Minimum z coordinate (z coordinate of node (0,0,0)
    T      m_max_x;      ///< Maximum x coordinate (x coordinate of node (I-1,J-1,K-1)
    T      m_max_y;      ///< Maximum y coordinate (y coordinate of node (I-1,J-1,K-1)
    T      m_max_z;      ///< Maximum z coordinate (z coordinate of node (I-1,J-1,K-1)
    T      m_dx;         ///< Space between nodes in x direction.
    T      m_dy;         ///< Space between nodes in y direction.
    T      m_dz;         ///< Space between nodes in z direction.
    size_t m_I;          ///< Number of nodes along x-axis.
    size_t m_J;          ///< Number of nodes along y-axis.
    size_t m_K;          ///< Number of nodes along z-axis.

    std::vector<D> m_data; ///< The data values stored at the grid nodes.

  protected:

    size_t index(size_t const & i,size_t const & j,size_t const & k) const
    {
      assert( i < this->m_I || !"index(): i was out of bounds");
      assert( j < this->m_J || !"index(): j was out of bounds");
      assert( k < this->m_K || !"index(): k was out of bounds");

      return (k*this->m_J + j)*this->m_I + i;
    }

  public:

    Grid()
    : m_min_x(VT::zero())
    , m_min_y(VT::zero())
    , m_min_z(VT::zero())
    , m_max_x(VT::zero())
    , m_max_y(VT::zero())
    , m_max_z(VT::zero())
    , m_dx(VT::zero())
    , m_dy(VT::zero())
    , m_dz(VT::zero())
    , m_I(0u)
    , m_J(0u)
    , m_K(0u)
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
        this->m_min_x = grid.m_min_x;
        this->m_min_y = grid.m_min_y;
        this->m_min_z = grid.m_min_z;
        this->m_max_x = grid.m_max_x;
        this->m_max_y = grid.m_max_y;
        this->m_max_z = grid.m_max_z;
        this->m_dx    = grid.m_dx;
        this->m_dy    = grid.m_dy;
        this->m_dz    = grid.m_dz;
        this->m_I     = grid.m_I;
        this->m_J     = grid.m_J;
        this->m_K     = grid.m_K;
        this->m_data  = grid.m_data;
      }
      return (*this);
    }

  public:

    void create(
                  T const & min_x
                , T const & min_y
                , T const & min_z
                , T const & max_x
                , T const & max_y
                , T const & max_z
                , size_t const & I
                , size_t const & J
                , size_t const & K
                )
    {
      assert( max_x > min_x || !"create(): max_x was less than equal min_x");
      assert( max_y > min_y || !"create(): max_y was less than equal min_y");
      assert( max_z > min_z || !"create(): max_z was less than equal min_z");
      assert( I>0           || !"create(): I must be positive");
      assert( J>0           || !"create(): J must be positive");
      assert( K>0           || !"create(): K must be positive");

      this->m_min_x = min_x;
      this->m_min_y = min_y;
      this->m_min_z = min_z;

      this->m_max_x = max_x;
      this->m_max_y = max_y;
      this->m_max_z = max_z;

      this->m_dx = (max_x-min_x)/(I-1);
      this->m_dy = (max_y-min_y)/(J-1);
      this->m_dz = (max_z-min_z)/(K-1);

      this->m_I = I;
      this->m_J = J;
      this->m_K = K;

      m_data.resize( I*J*K );
    }

    D & operator()(size_t const & i,size_t const & j,size_t const & k)
    {
      assert(this->m_data.size()>0 || !"operator(): no data");

      return this->m_data[ this->index(i,j,k) ];
    }

    D const & operator() (size_t const & i,size_t const & j,size_t const & k) const
    {
      assert(this->m_data.size()>0 || !"operator(): no data");

      return this->m_data[ this->index(i,j,k) ];
    }

    T width()  const { return this->m_max_x - this->m_min_x; }
    T height() const { return this->m_max_y - this->m_min_y; }
    T depth()  const { return this->m_max_z - this->m_min_z; }

    size_t size()  const { return this->m_data.size(); }
    bool   empty() const { return this->m_data.empty(); }

    T const & min_x() const { return this->m_min_x; }
    T const & min_y() const { return this->m_min_y; }
    T const & min_z() const { return this->m_min_z; }

    T const & max_x() const { return this->m_max_x; }
    T const & max_y() const { return this->m_max_y; }
    T const & max_z() const { return this->m_max_z; }

    T const & dx() const { return this->m_dx; }
    T const & dy() const { return this->m_dy; }
    T const & dz() const { return this->m_dz; }

    size_t const & I() const { return this->m_I; }
    size_t const & J() const { return this->m_J; }
    size_t const & K() const { return this->m_K; }

    D       * data_ptr()       { return this->m_data[0]; }
    D const * data_ptr() const { return this->m_data[0]; }

  };

} // namespace grid

// GRID_GRID_H
#endif
