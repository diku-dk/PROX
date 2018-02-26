#ifndef HYPER_MATH_POLICY_H
#define HYPER_MATH_POLICY_H

#include <hyper_element_matrices_4x4.h>

#include <sparse.h>
#include <tiny.h>
#include <mesh_array.h>

#include <vector>
#include <cassert>

namespace hyper
{

  template< typename T >
  class MathPolicy
  : public tiny::MathTypes<T>
  {
  public:

    typedef          tiny::MathTypes<T>                            base_type;
    typedef          T                                             real_type;
    typedef typename base_type::vector3_type                       vector3_type;
    typedef typename base_type::matrix3x3_type                     matrix3x3_type;

    typedef          sparse::Block<3,1,T>                          block3x1_type;
    typedef          sparse::Block<3,3,T>                          block3x3_type;

    typedef          sparse::Vector< block3x1_type >               vector_block3x1_type;
    typedef          sparse::DiagonalMatrix< block3x3_type >       diagonal_block3x3_type;
    typedef          sparse::CompressedRowMatrix< block3x3_type >  compressed_block3x3_type;
    typedef          ElementMatrices4x4< matrix3x3_type >          element_matrices_type;

  public:

    static block3x3_type convert( matrix3x3_type const & M)
    {
      block3x3_type B;

      B(0,0) = M(0,0);      B(0,1) = M(0,1);      B(0,2) = M(0,2);
      B(1,0) = M(1,0);      B(1,1) = M(1,1);      B(1,2) = M(1,2);
      B(2,0) = M(2,0);      B(2,1) = M(2,1);      B(2,2) = M(2,2);

      return B;
    }

    static matrix3x3_type convert( block3x3_type const & B)
    {
      matrix3x3_type M;

      M(0,0) = B(0,0);      M(0,1) = B(0,1);      M(0,2) = B(0,2);
      M(1,0) = B(1,0);      M(1,1) = B(1,1);      M(1,2) = B(1,2);
      M(2,0) = B(2,0);      M(2,1) = B(2,1);      M(2,2) = B(2,2);

      return M;
    }

    static block3x1_type convert( vector3_type const & V)
    {
      block3x1_type B;

      B(0) = V(0);
      B(1) = V(1);
      B(2) = V(2);

      return B;
    }

    static vector3_type convert( block3x1_type const & B)
    {
      vector3_type V;

      V(0) = B(0);
      V(1) = B(1);
      V(2) = B(2);

      return V;
    }

    static vector_block3x1_type convert( std::vector<vector3_type> const & data )
    {
      vector_block3x1_type V;

      V.resize(data.size());

      typename std::vector<vector3_type>::const_iterator v   = data.begin();
      typename std::vector<vector3_type>::const_iterator end = data.end();

      for(unsigned int i = 0u;v!=end;++v,++i)
        V(i) = convert(*v);
    }

    static std::vector<vector3_type> convert(vector_block3x1_type const & data )
    {
      std::vector<vector3_type> V;

      V.resize(data.size());

      typename std::vector<vector3_type>::const_iterator v   = data.begin();
      typename std::vector<vector3_type>::const_iterator end = data.end();
      
      for(unsigned int i = 0u;v!=end;++v,++i)
        V(i) = convert(*v);
    }

    static vector_block3x1_type convert(
                                        mesh_array::VertexAttribute<T,mesh_array::T4Mesh> const & x
                                        , mesh_array::VertexAttribute<T,mesh_array::T4Mesh> const & y
                                        , mesh_array::VertexAttribute<T,mesh_array::T4Mesh> const & z
                                        )
    {
      vector_block3x1_type result;

      unsigned int const N = x.size();

      result.resize(N);

      for(unsigned int i = 0u; i < N; ++i)
      {
        result(i)(0) = x(i);
        result(i)(1) = y(i);
        result(i)(2) = z(i);
      }

      return result;
    }

    static void convert(
                        vector_block3x1_type const & data
                        , mesh_array::VertexAttribute<T,mesh_array::T4Mesh>  & x
                        , mesh_array::VertexAttribute<T,mesh_array::T4Mesh>  & y
                        , mesh_array::VertexAttribute<T,mesh_array::T4Mesh>  & z
                        )
    {
      unsigned int const N = data.size();

      assert(x.size() == N || !"convert(): x vertex attribute size was not correct");
      assert(y.size() == N || !"convert(): y vertex attribute size was not correct");
      assert(z.size() == N || !"convert(): z vertex attribute size was not correct");

      for(unsigned int i = 0u; i < N; ++i)
      {
        x(i) = data(i)(0);
        y(i) = data(i)(1);
        z(i) = data(i)(2);
      }
    }

  };

} // namespace hyper

// HYPER_MATH_POLICY_H
#endif 
