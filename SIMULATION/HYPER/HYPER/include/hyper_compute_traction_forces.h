#ifndef HYPER_COMPUTE_TRACTION_FORCES_H
#define HYPER_COMPUTE_TRACTION_FORCES_H

#include <hyper_traction_info.h>
#include <hyper_math_policy.h>

#include <mesh_array.h>

#include <vector>
#include <cassert>

namespace hyper
{

  /**
   *
   * 2015_03_24 Kenny code review: This implementation is inherently sequential.
   * It will be "slow". One could exploit doing a scatter of traction forces
   * into an element-wise vector and then afterwards run a gather operation
   * accumulating element-wise forces into total nodal forces. This will result
   * in two embarrassingly parallel for-loops. Only down-side is that we need
   * to compute "neighborhood" information of the triangle surface. We can not
   * re-use tetrahedral neighborhood information, sadly:-(
   */
  template<typename MT>
  inline void compute_traction_forces(
                                      std::vector<TractionInfo<MT> > const & surface
                                      , mesh_array::T4Mesh  const & mesh
                                      , typename MT::vector_block3x1_type const & x
                                      , typename MT::vector_block3x1_type & f
                                      )
  {
    typedef typename MT::vector_block3x1_type    vector_block3x1_type;
    typedef typename MT::matrix3x3_type          M;
    typedef typename MT::vector3_type            V;
    typedef typename MT::real_type               T;
    typedef typename MT::value_traits            VT;

    unsigned int const N = mesh.vertex_size();

    f.resize( N );
    f.clear_data();

    unsigned int const K = surface.size();

    for(unsigned int e = 0u; e < K;++e)
    {
      unsigned int const i = surface[e].i();
      unsigned int const j = surface[e].j();
      unsigned int const k = surface[e].k();

      V const & t = surface[e].traction();

      V const & x_i = MT::convert( x[i] );
      V const & x_j = MT::convert( x[j] );
      V const & x_k = MT::convert( x[k] );

      V const A_vec = cross( x_j - x_i, x_k - x_i) / VT::two();
      T const A     = norm(A_vec);

      //--- Solving integral of products of shape functions for isoparametric ----
      //--- linear triangle element ----------------------------------------------
      //
      //   syms x y real
      //    Ni = 1 - x - y
      //    Nj = y
      //    Nk = x
      //    N = [Ni Nj Nk]
      //    A = N'* N
      //    AA = int(A,x,0,1-y)
      //    L = int(AA,y,0,1)
      //
      // We find the load (nodal traction distribution)  matrix
      //
      //  L = 1 / 24  *  [ 2*I_3x3    I_3x3   I_3x3;
      //                     I_3x3  2*I_3x3   I_3x3;
      //                     I_3x3    I_3x3 2*I_3x3; ]
      //--- The spatial load force -----------------------------------------------
      //
      //    lf(f) = A(f) * L * [Ti; Tj; Tk]
      //
      //  where T is the nodal surface traction.

      V const & ti = t;
      V const & tj = t;
      V const & tk = t;

      M const    I = M::identity();

      V const fi = A * ( VT::two()*I*ti +           I*tj +           I*tk )  / VT::numeric_cast(24.0);
      V const fj = A * (           I*ti + VT::two()*I*tj +           I*tk )  / VT::numeric_cast(24.0);
      V const fk = A * (           I*ti +           I*tj + VT::two()*I*tk )  / VT::numeric_cast(24.0);

      f[i] = MT::convert( fi + MT::convert(f[i]) );
      f[j] = MT::convert( fj + MT::convert(f[j]) );
      f[k] = MT::convert( fk + MT::convert(f[k]) );
    }
  }

}// namespace hyper

// HYPER_COMPUTE_TRACTION_FORCES_H
#endif
