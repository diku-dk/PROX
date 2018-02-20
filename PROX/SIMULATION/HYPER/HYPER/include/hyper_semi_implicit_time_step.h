#ifndef HYPER_SEMI_IMPLICIT_STEP_H
#define HYPER_SEMI_IMPLICIT_STEP_H

#include <hyper_math_policy.h>
#include <hyper_constitutive_equation.h>
#include <hyper_compute_elastic_forces.h>
#include <hyper_compute_mass_matrix.h>
#include <hyper_compute_damping_matrix.h>
#include <hyper_assemble_diagonal.h>
#include <hyper_assemble_matrix.h>
#include <hyper_apply_dirichlet_conditions.h>

#include <mesh_array.h>

#include <sparse.h>
#include <sparse_preconditioners.h>
#include <sparse_conjugate_gradient.h>

#include <vector>
#include <cassert>

namespace hyper
{

  namespace details
  {
    
    /**
     * Computes:  b = M*v + dt*(Fext - f - C*v)
     */
    template<typename MT>
    inline void compute_rhs(
                            mesh_array::T4Mesh  const & mesh
                            , mesh_array::Neighborhoods const & neighbors
                            , typename MT::real_type const & dt
                            , std::vector< typename MT::element_matrices_type > const & Me
                            , std::vector< typename MT::element_matrices_type > const & Ce
                            , typename MT::vector_block3x1_type const & v
                            , typename MT::vector_block3x1_type const & Fext
                            , typename MT::vector_block3x1_type const & f
                            , typename MT::vector_block3x1_type & b
                            )
    {
      assert( v.size()    == Fext.size() || !"compute_rhs(): v and Fext must have same size");
      assert( Fext.size() == f.size()    || !"compute_rhs(): Fext and f must have same size");
      assert( v.size()    == b.size()    || !"compute_rhs(): v and b must have same size");

      typedef typename MT::vector_block3x1_type      vector_block3x1_type;
//      typedef typename MT::diagonal_block3x3_type    diagonal_block3x3_type;
//      typedef typename MT::compressed_block3x3_type  compressed_block3x3_type;
//      typedef typename MT::element_matrices_type     element_matrices_type;

//      typedef typename MT::matrix3x3_type            M;
      typedef typename MT::vector3_type              V;
//      typedef typename MT::real_type                 T;
//      typedef typename MT::value_traits              VT;

      unsigned int const N = v.size();  ///< Number of blocks

      vector_block3x1_type cv(N);
      vector_block3x1_type mv(N);

      mul<MT>( mesh, neighbors, Ce, v, cv);
      mul<MT>( mesh, neighbors, Me, v, mv);

      for(size_t i = 0u; i < N; ++i)
      {
        V const Fext_i  = MT::convert( Fext[i] );
        V const f_i     = MT::convert( f[i] );
        V const cv_i    = MT::convert( cv[i] );
        V const mv_i    = MT::convert( mv[i] );
        V const b_i     = mv_i + dt*(Fext_i - f_i - cv_i);
        b[i]            = MT::convert( b_i );
      }
    }

  }// end namespace details


  /**
   * Preconditioner types supported by the semi_implicit_time_step method.
   */
  typedef enum { identity_precond, jacobi_precond, block_jacobi_precond} preconditioner_type;

  /**
   *
   */
  template<typename MT>
  inline void semi_implicit_time_step(
                                      mesh_array::T4Mesh  const & mesh
                                      , mesh_array::Neighborhoods const & neighbors
                                      , std::vector<DirichletInfo<MT> > const & dirichlet_conditions
                                      , typename MT::real_type const & dt
                                      , typename MT::vector_block3x1_type & x
                                      , typename MT::vector_block3x1_type const & x0
                                      , typename MT::vector_block3x1_type & v
                                      , typename MT::vector_block3x1_type const & Fext
                                      , ConstitutiveEquation<MT> const * model
                                      , preconditioner_type const & precond
                                      )
  {
    assert( model || !"semi_implicit_time_step(): model was NULL");

//    typedef typename MT::matrix3x3_type            M;
//    typedef typename MT::vector3_type              V;
//    typedef typename MT::real_type                 T;
//    typedef typename MT::value_traits              VT;

    typedef typename MT::vector_block3x1_type      vector_block3x1_type;
    typedef typename MT::diagonal_block3x3_type    diagonal_block3x3_type;
    typedef typename MT::compressed_block3x3_type  compressed_block3x3_type;
    typedef typename MT::element_matrices_type     element_matrices_type;

    unsigned int const N = x.size();  ///< Number of blocks

    vector_block3x1_type f(N);
    vector_block3x1_type b(N);
    vector_block3x1_type tmp(N);

    std::vector< element_matrices_type > Me;
    std::vector< element_matrices_type > Ce;

    compute_elastic_forces<MT>( mesh, neighbors, x, x0, f, model );

    compute_mass_matrix<MT>( mesh, x0, model, Me);         // 2015-03-25 Kenny code review: Only needs to be recomputed if mesh changed
    compute_damping_matrix<MT>( mesh, x0, model, Ce);      // 2015-03-25 Kenny code review: Only needs to be recomputed if mesh changed

    //--- Velocity update ------------------------------------------------------
    //--- solve
    //---
    //---   A v^{t+1} = b
    //---
    //--- where
    //---
    //---  A = M
    //---  b = M*v^t + dt*(Fext - f - C*v^t)
    //---
    details::compute_rhs<MT>(mesh, neighbors, dt, Me, Ce, v, Fext, f, b );

    if( model->material().lumped() )
    {
      diagonal_block3x3_type D;

      assemble_diagonal<MT>(mesh ,neighbors, Me, D);

      apply_dirichlet_conditions<MT>(dirichlet_conditions, D, b);

      sparse::inverse(D);

      v.clear_data();

      sparse::prod(D,b,v);
    }
    else
    {
      compressed_block3x3_type A;

      assemble_matrix<MT>(mesh, neighbors, Me,  A);

      apply_dirichlet_conditions<MT>(dirichlet_conditions, A, b);

      compressed_block3x3_type P;

      switch(precond)
      {
        case identity_precond:     sparse::make_identity_preconditioner(P);       break;
        case jacobi_precond:       sparse::make_jacobi_preconditioner(A,P);       break;
        case block_jacobi_precond: sparse::make_block_jacobi_preconditioner(A,P); break;
      };

      sparse::preconditioned_conjugate_gradient(P,A,v,b);
    }

    //--- Position update ------------------------------------------------------
    //    x = x + dt*v;
    sparse::prod(dt, v, tmp);
    sparse::add(tmp, x);
  }
  
}// namespace hyper

// HYPER_SEMI_IMPLICIT_TIME_STEP_H
#endif
