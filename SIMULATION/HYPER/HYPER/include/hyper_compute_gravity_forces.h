#ifndef HYPER_COMPUTE_GRAVITY_FORCES_H
#define HYPER_COMPUTE_GRAVITY_FORCES_H

#include <hyper_math_policy.h>
#include <hyper_constitutive_equation.h>
#include <hyper_compute_mass_matrix.h>

#include <mesh_array.h>


#include <vector>
#include <cassert>

namespace hyper
{

  /**
   * Computes Gravity Forces
   *
   * Fext = M*g
   *
   */
  template<typename MT>
  inline void compute_gravity_forces(
                                     mesh_array::T4Mesh  const & mesh
                                     , mesh_array::Neighborhoods const & neighbors
                                     , typename MT::vector3_type const & g
                                     , typename MT::vector_block3x1_type const & x0
                                     , ConstitutiveEquation<MT> const * model
                                     , typename MT::vector_block3x1_type & Fext
                                     )
  {



    assert( x0.size()    == Fext.size() || !"compute_gravity_forces(): v and Fext must have same size");

    typedef typename MT::vector_block3x1_type      vector_block3x1_type;
    typedef typename MT::diagonal_block3x3_type    diagonal_block3x3_type;
    typedef typename MT::compressed_block3x3_type  compressed_block3x3_type;
    typedef typename MT::element_matrices_type     element_matrices_type;

    typedef typename MT::matrix3x3_type            M;
    typedef typename MT::vector3_type              V;
    typedef typename MT::real_type                 T;
    typedef typename MT::value_traits              VT;

    unsigned int const N = x0.size();  ///< Number of blocks

    std::vector< element_matrices_type > Me;
    compute_mass_matrix<MT>( mesh, x0, model, Me);

    vector_block3x1_type G(N);

    for(size_t i = 0u; i < N; ++i)
      G[i] = MT::convert( g );

    mul<MT>( mesh, neighbors, Me, G, Fext);
  }

}// namespace hyper

// HYPER_COMPUTE_GRAVITY_FORCES_H
#endif
