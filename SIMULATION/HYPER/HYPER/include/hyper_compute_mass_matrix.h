#ifndef HYPER_COMPUTE_MASS_MATRIX_H
#define HYPER_COMPUTE_MASS_MATRIX_H

#include <hyper_constitutive_equation.h>
#include <hyper_math_policy.h>

#include <mesh_array.h>

#include <vector>
#include <cassert>

namespace hyper
{

  /**
   * Compute elementwise mass matrix information.
   *
   * @tparam MT    Math types type binder.
   */
  template<typename MT>
  inline void compute_mass_matrix(
                                  mesh_array::T4Mesh  const & mesh
                                  , typename MT::vector_block3x1_type const & x
                                  , ConstitutiveEquation<MT> const * model
                                  , std::vector< typename MT::element_matrices_type > & Me
                                  )
  {
    assert( model || !"compute_mass_matrix(): model was NULL");

    typedef typename MT::vector_block3x1_type    vector_block3x1_type;
    typedef typename MT::matrix3x3_type          M;
    typedef typename MT::vector3_type            V;
    typedef typename MT::real_type               T;
    typedef typename MT::value_traits            VT;
    typedef typename MT::element_matrices_type   element_matrices_type;

    element_matrices_type NN;

    if( model->material().lumped() )
    {
      // Lumped mass matrix
      //
      //  Me = rho V / 4  *  [ I_3x3    0       0       0;
      //                       0        I_3x3   0       0;
      //                       0        0       I_3x3   0;
      //                       0        0       0       I_3x3;]
      //
      NN = make_diag_identity<M>() / VT::four();
    }
    else
    {
      // Consistent mass matrix
      //
      //  Me = rho V / 20  *  [ 2*I_3x3    I_3x3   I_3x3   I_3x3;
      //                          I_3x3  2*I_3x3   I_3x3   I_3x3;
      //                          I_3x3    I_3x3 2*I_3x3   I_3x3;
      //                          I_3x3    I_3x3   I_3x3 2*I_3x3;]
      //
      NN = (make_all_identity<M>() + make_diag_identity<M>()) /  VT::numeric_cast(20.0) ;
    }

    T const rho = model->material().rho();

    //--- Allocate space for all mass matrix blocks ----------------------------
    unsigned int const K = mesh.tetrahedron_size();

    Me.resize( K );

    //--- Now compute element mass matrices ------------------------------------
    for(unsigned int idx=0u; idx < K; ++idx)
    {
      mesh_array::Tetrahedron const & tetrahedron = mesh.tetrahedron(idx);

      unsigned int const i  = tetrahedron.i();
      unsigned int const j  = tetrahedron.j();
      unsigned int const k  = tetrahedron.k();
      unsigned int const m  = tetrahedron.m();

      //--- Extract coordinates of tetrahedron vertices ------------------------
      V const xi = MT::convert( x[i] );
      V const xj = MT::convert( x[j] );
      V const xk = MT::convert( x[k] );
      V const xm = MT::convert( x[m] );

      V const u_ji = xj - xi;
      V const u_ki = xk - xi;
      V const u_mi = xm - xi;

      T const vol = inner_prod(u_mi ,cross(u_ji,u_ki) ) / VT::numeric_cast(6.0);

      assert( vol > VT::zero() || !"compute_mass_matrix(): degenerate tetrahedron found, negative volume");

      Me[idx] = NN * (vol * rho );
    }

  }

}// namespace hyper

// HYPER_COMPUTE_MASS_MATRIX_H
#endif
