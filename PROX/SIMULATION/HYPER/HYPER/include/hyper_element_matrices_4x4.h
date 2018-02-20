#ifndef HYPER_ELEMENT_MATRICES_4X4_H
#define HYPER_ELEMENT_MATRICES_4X4_H

#include <mesh_array.h>
#include <tiny.h>

namespace hyper
{

  template<typename M>
  class ElementMatrices4x4
  {
  protected:

    M m_data[16u];

  public:

    ElementMatrices4x4()
    {
      for (unsigned int i=0;i<16u;++i)
      {
        this->m_data[i].clear();
      }
    }

    ElementMatrices4x4(ElementMatrices4x4 const & cpy)
    {
      *this = cpy;
    }

    ElementMatrices4x4 const & operator=(ElementMatrices4x4 const & rhs)
    {
      if( this != &rhs)
      {
        for (unsigned int i=0;i<16u;++i)
        {
          this->m_data[i] = rhs.m_data[i];
        }
      }
      return *this;
    }

  public:

    M const & operator()(unsigned int const & i,unsigned int const & j) const
    {
      assert( i < 4u || !"ElementMatrices4x4::operator(): index out of bound");
      assert( j < 4u || !"ElementMatrices4x4::operator(): index out of bound");

      return this->m_data[j*4u + i];
    }

    M       & operator()(unsigned int const & i,unsigned int const & j)
    {
      assert( i < 4u || !"ElementMatrices4x4::operator(): index out of bound");
      assert( j < 4u || !"ElementMatrices4x4::operator(): index out of bound");

      return this->m_data[j*4u + i];
    }

    M const & data(unsigned int const & i) const
    {
      assert( i < 16u || !"ElementMatrices4x4::data(): index out of bound");

      return this->m_data[i];
    }

    M & data(unsigned int const & i)
    {
      assert( i < 16u || !"ElementMatrices4x4::data(): index out of bound");

      return this->m_data[i];
    }

  public:

    static ElementMatrices4x4<M> zero()
    {
      ElementMatrices4x4<M> B;
      B(0,0) = M::zero();
      B(1,1) = M::zero();
      B(2,2) = M::zero();
      B(3,3) = M::zero();
      return B;
    }

  };

  template<typename M>
  inline ElementMatrices4x4<M> operator/(ElementMatrices4x4<M> const & A, typename M::real_type const & s)
  {
    ElementMatrices4x4<M> B;
    for (unsigned int i=0;i<16u;++i)
    {
      B.data(i) = A.data(i) / s;
    }
    return B;
  }

  template<typename M>
  inline ElementMatrices4x4<M> operator*(ElementMatrices4x4<M> const & A, typename M::real_type const & s)
  {
    ElementMatrices4x4<M> B;
    for (unsigned int i=0;i<16u;++i)
    {
      B.data(i) = A.data(i) * s;
    }
    return B;
  }

  template<typename M>
  inline ElementMatrices4x4<M> operator+(ElementMatrices4x4<M> const & A, ElementMatrices4x4<M> const & B)
  {
    ElementMatrices4x4<M> C;
    for (unsigned int i=0;i<16u;++i)
    {
      C.data(i) = B.data(i) + A.data(i);
    }
    return C;
  }

  template<typename M>
  inline ElementMatrices4x4<M> make_diag_identity()
  {
    ElementMatrices4x4<M> B;
    B(0,0) = M::identity();
    B(1,1) = M::identity();
    B(2,2) = M::identity();
    B(3,3) = M::identity();
    return B;
  }

  template<typename M>
  inline ElementMatrices4x4<M> make_all_identity()
  {
    ElementMatrices4x4<M> B;
    for (unsigned int i=0;i<4u;++i)
    {
      for (unsigned int j=0;j<4u;++j)
      {
        B(i,j) = M::identity();
      }
    }
    return B;
  }


  /**
   * Computes y = A x
   */
  template<typename MT>
  inline void mul(
                  mesh_array::T4Mesh  const & mesh
                  , mesh_array::Neighborhoods const & neighbors
                  , std::vector< typename MT::element_matrices_type > const & A
                  , typename MT::vector_block3x1_type const & x
                  , typename MT::vector_block3x1_type & y
                  )
  {
    assert( x.size() == y.size()                || !"mul(): x and y must have same size");
    assert( A.size() > 0u                       || !"mul(): A was empty");
    assert( mesh.tetrahedron_size() == A.size() || !"mul(): mesh did not match A?");
    assert( mesh.vertex_size() == x.size()      || !"mul(): mesh did not match x?");

//    typedef typename MT::vector_block3x1_type      vector_block3x1_type;
//    typedef typename MT::diagonal_block3x3_type    diagonal_block3x3_type;
//    typedef typename MT::compressed_block3x3_type  compressed_block3x3_type;
    typedef typename MT::element_matrices_type     element_matrices_type;

//    typedef typename MT::matrix3x3_type            M;
    typedef typename MT::vector3_type              V;
//    typedef typename MT::real_type                 T;
//    typedef typename MT::value_traits              VT;

    for(size_t i = 0u; i < y.size(); ++i)
    {

      V sum = V::zero();

      for (unsigned int entry = neighbors.m_offset[i]; entry < neighbors.m_offset[i+1];++entry)
      {
        unsigned int            const   idx         = neighbors.m_V2T[ entry ].second;
        mesh_array::Tetrahedron const & tetrahedron = mesh.tetrahedron(idx);
        unsigned int            const   local_idx   = tetrahedron.get_local_index(i);

        element_matrices_type const & Ae = A[idx];

        V const xi = MT::convert( x[tetrahedron.i()] );
        V const xj = MT::convert( x[tetrahedron.j()] );
        V const xk = MT::convert( x[tetrahedron.k()] );
        V const xm = MT::convert( x[tetrahedron.m()] );

        V const y0 = prod( Ae( local_idx, 0 ), xi );
        V const y1 = prod( Ae( local_idx, 1 ), xj );
        V const y2 = prod( Ae( local_idx, 2 ), xk );
        V const y3 = prod( Ae( local_idx, 3 ), xm );

        sum += y0 + y1 + y2 + y3;
      }

      y[i] = MT::convert( sum );

    }
  }

}// namespace hyper

// HYPER_ELEMENT_MATRICES_4X4_H
#endif
