#ifndef HYPER_COMPUTE_ELASTIC_FORCES_H
#define HYPER_COMPUTE_ELASTIC_FORCES_H

#include <hyper_constitutive_equation.h>

#include <hyper_strain_tensors.h>

#include <hyper_math_policy.h>

#include <mesh_array.h>
#include <tiny_is_finite.h>
#include <tiny_is_number.h>

#include <vector>
#include <cassert>

namespace hyper
{

  /**
   * Compute elastic forces for the given mesh.
   *
   * @tparam MT    Math types type binder.
   */
  template<typename MT>
  inline void compute_elastic_forces(
                                     mesh_array::T4Mesh  const & mesh
                                     , mesh_array::Neighborhoods const & neighbors
                                     , typename MT::vector_block3x1_type const & x
                                     , typename MT::vector_block3x1_type const & x0
                                     , typename MT::vector_block3x1_type & f
                                     , ConstitutiveEquation<MT> const * model
                                     )
  {
    assert( model || !"compute_elastic_forces(): model was NULL");

    typedef typename MT::vector_block3x1_type    vector_block3x1_type;
    typedef typename MT::matrix3x3_type          M;
    typedef typename MT::vector3_type            V;
    typedef typename MT::real_type               T;
    typedef typename MT::value_traits            VT;

    unsigned int const K = mesh.tetrahedron_size();

    std::vector<V> local_forces;

    local_forces.resize( K*4u );

    //--- Compute local elastic forces per tetrahedral element -----------------
    for ( unsigned int idx = 0u; idx < K; ++idx)  // 2014-01-22 Kenny: This for-loop is naive parallel
    {
      mesh_array::Tetrahedron const & tetrahedron = mesh.tetrahedron(idx);

      unsigned int const i  = tetrahedron.i();
      unsigned int const j  = tetrahedron.j();
      unsigned int const k  = tetrahedron.k();
      unsigned int const m  = tetrahedron.m();

      assert(i!=j || !"compute_elastic_forces(): i == j");
      assert(i!=k || !"compute_elastic_forces(): i == k");
      assert(i!=m || !"compute_elastic_forces(): i == m");
      assert(j!=k || !"compute_elastic_forces(): j == k");
      assert(j!=m || !"compute_elastic_forces(): j == m");
      assert(m!=k || !"compute_elastic_forces(): m == k");

      //--- Extract coordinates of tetrahedron vertices ------------------------
      V const xi = MT::convert( x[i] );
      V const xj = MT::convert( x[j] );
      V const xk = MT::convert( x[k] );
      V const xm = MT::convert( x[m] );

      V const x0i = MT::convert( x0[i] );
      V const x0j = MT::convert( x0[j] );
      V const x0k = MT::convert( x0[k] );
      V const x0m = MT::convert( x0[m] );

      //--- Compute deformation gradient ---------------------------------------
      V const u_ji = xj - xi;
      V const u_ki = xk - xi;
      V const u_mi = xm - xi;

      assert( is_number(u_ji(0)) || !"compute_elastic_forces(): u_ji(0) is not a number");
      assert( is_number(u_ji(1)) || !"compute_elastic_forces(): u_ji(1) is not a number");
      assert( is_number(u_ji(2)) || !"compute_elastic_forces(): u_ji(2) is not a number");

      assert( is_finite(u_ji(0)) || !"compute_elastic_forces(): u_ji(0) is not finite"  );
      assert( is_finite(u_ji(1)) || !"compute_elastic_forces(): u_ji(1) is not finite"  );
      assert( is_finite(u_ji(2)) || !"compute_elastic_forces(): u_ji(2) is not finite"  );

      assert( is_number(u_ki(0)) || !"compute_elastic_forces(): u_ki(0) is not a number");
      assert( is_number(u_ki(1)) || !"compute_elastic_forces(): u_ki(1) is not a number");
      assert( is_number(u_ki(2)) || !"compute_elastic_forces(): u_ki(2) is not a number");

      assert( is_finite(u_ki(0)) || !"compute_elastic_forces(): u_ki(0) is not finite"  );
      assert( is_finite(u_ki(1)) || !"compute_elastic_forces(): u_ki(1) is not finite"  );
      assert( is_finite(u_ki(2)) || !"compute_elastic_forces(): u_ki(2) is not finite"  );

      assert( is_number(u_mi(0)) || !"compute_elastic_forces(): u_mi(0) is not a number");
      assert( is_number(u_mi(1)) || !"compute_elastic_forces(): u_mi(1) is not a number");
      assert( is_number(u_mi(2)) || !"compute_elastic_forces(): u_mi(2) is not a number");

      assert( is_finite(u_mi(0)) || !"compute_elastic_forces(): u_mi(0) is not finite"  );
      assert( is_finite(u_mi(1)) || !"compute_elastic_forces(): u_mi(1) is not finite"  );
      assert( is_finite(u_mi(2)) || !"compute_elastic_forces(): u_mi(2) is not finite"  );

      assert( inner_prod( u_mi, cross(u_ji, u_ki)) > VT::zero() || !"compute_elastic_forces(): degenerate tetrahedron in spatial space");

      V const u0_ji = x0j - x0i;
      V const u0_ki = x0k - x0i;
      V const u0_mi = x0m - x0i;

      assert( is_number(u0_ji(0)) || !"compute_elastic_forces(): u0_ji(0) is not a number");
      assert( is_number(u0_ji(1)) || !"compute_elastic_forces(): u0_ji(1) is not a number");
      assert( is_number(u0_ji(2)) || !"compute_elastic_forces(): u0_ji(2) is not a number");

      assert( is_finite(u0_ji(0)) || !"compute_elastic_forces(): u0_ji(0) is not finite"  );
      assert( is_finite(u0_ji(1)) || !"compute_elastic_forces(): u0_ji(1) is not finite"  );
      assert( is_finite(u0_ji(2)) || !"compute_elastic_forces(): u0_ji(2) is not finite"  );

      assert( is_number(u0_ki(0)) || !"compute_elastic_forces(): u0_ki(0) is not a number");
      assert( is_number(u0_ki(1)) || !"compute_elastic_forces(): u0_ki(1) is not a number");
      assert( is_number(u0_ki(2)) || !"compute_elastic_forces(): u0_ki(2) is not a number");

      assert( is_finite(u0_ki(0)) || !"compute_elastic_forces(): u0_ki(0) is not finite"  );
      assert( is_finite(u0_ki(1)) || !"compute_elastic_forces(): u0_ki(1) is not finite"  );
      assert( is_finite(u0_ki(2)) || !"compute_elastic_forces(): u0_ki(2) is not finite"  );

      assert( is_number(u0_mi(0)) || !"compute_elastic_forces(): u0_mi(0) is not a number");
      assert( is_number(u0_mi(1)) || !"compute_elastic_forces(): u0_mi(1) is not a number");
      assert( is_number(u0_mi(2)) || !"compute_elastic_forces(): u0_mi(2) is not a number");

      assert( is_finite(u0_mi(0)) || !"compute_elastic_forces(): u0_mi(0) is not finite"  );
      assert( is_finite(u0_mi(1)) || !"compute_elastic_forces(): u0_mi(1) is not finite"  );
      assert( is_finite(u0_mi(2)) || !"compute_elastic_forces(): u0_mi(2) is not finite"  );

      assert( inner_prod( u0_mi, cross(u0_ji, u0_ki)) > VT::zero() || !"compute_elastic_forces(): degenerate tetrahedron in spatial space");

      M const D = M::make(
                          u_ji(0), u_ki(0), u_mi(0)
                        , u_ji(1), u_ki(1), u_mi(1)
                        , u_ji(2), u_ki(2), u_mi(2)
                        );

      M const D0 = M::make(
                          u0_ji(0), u0_ki(0), u0_mi(0)
                          , u0_ji(1), u0_ki(1), u0_mi(1)
                          , u0_ji(2), u0_ki(2), u0_mi(2)
                          );

      M const invD0 = inverse(D0);

      assert( is_number(invD0(0,0)) || !"compute_elastic_forces(): F(0,0) is not a number");
      assert( is_number(invD0(0,1)) || !"compute_elastic_forces(): F(0,1) is not a number");
      assert( is_number(invD0(0,2)) || !"compute_elastic_forces(): F(0,2) is not a number");
      assert( is_number(invD0(1,0)) || !"compute_elastic_forces(): F(1,0) is not a number");
      assert( is_number(invD0(1,1)) || !"compute_elastic_forces(): F(1,1) is not a number");
      assert( is_number(invD0(1,2)) || !"compute_elastic_forces(): F(1,2) is not a number");
      assert( is_number(invD0(2,0)) || !"compute_elastic_forces(): F(2,0) is not a number");
      assert( is_number(invD0(2,1)) || !"compute_elastic_forces(): F(2,1) is not a number");
      assert( is_number(invD0(2,2)) || !"compute_elastic_forces(): F(2,2) is not a number");

      assert( is_finite(invD0(0,0)) || !"compute_elastic_forces(): invD0(0,0) is not finite");
      assert( is_finite(invD0(0,1)) || !"compute_elastic_forces(): invD0(0,1) is not finite");
      assert( is_finite(invD0(0,2)) || !"compute_elastic_forces(): invD0(0,2) is not finite");
      assert( is_finite(invD0(1,0)) || !"compute_elastic_forces(): invD0(1,0) is not finite");
      assert( is_finite(invD0(1,1)) || !"compute_elastic_forces(): invD0(1,1) is not finite");
      assert( is_finite(invD0(1,2)) || !"compute_elastic_forces(): invD0(1,2) is not finite");
      assert( is_finite(invD0(2,0)) || !"compute_elastic_forces(): invD0(2,0) is not finite");
      assert( is_finite(invD0(2,1)) || !"compute_elastic_forces(): invD0(2,1) is not finite");
      assert( is_finite(invD0(2,2)) || !"compute_elastic_forces(): invD0(2,2) is not finite");

      M const F  = prod(D, invD0 );

      assert( is_number(F(0,0)) || !"compute_elastic_forces(): F(0,0) is not a number");
      assert( is_number(F(0,1)) || !"compute_elastic_forces(): F(0,1) is not a number");
      assert( is_number(F(0,2)) || !"compute_elastic_forces(): F(0,2) is not a number");
      assert( is_number(F(1,0)) || !"compute_elastic_forces(): F(1,0) is not a number");
      assert( is_number(F(1,1)) || !"compute_elastic_forces(): F(1,1) is not a number");
      assert( is_number(F(1,2)) || !"compute_elastic_forces(): F(1,2) is not a number");
      assert( is_number(F(2,0)) || !"compute_elastic_forces(): F(2,0) is not a number");
      assert( is_number(F(2,1)) || !"compute_elastic_forces(): F(2,1) is not a number");
      assert( is_number(F(2,2)) || !"compute_elastic_forces(): F(2,2) is not a number");

      assert( is_finite(F(0,0)) || !"compute_elastic_forces(): F(0,0) is not finite");
      assert( is_finite(F(0,1)) || !"compute_elastic_forces(): F(0,1) is not finite");
      assert( is_finite(F(0,2)) || !"compute_elastic_forces(): F(0,2) is not finite");
      assert( is_finite(F(1,0)) || !"compute_elastic_forces(): F(1,0) is not finite");
      assert( is_finite(F(1,1)) || !"compute_elastic_forces(): F(1,1) is not finite");
      assert( is_finite(F(1,2)) || !"compute_elastic_forces(): F(1,2) is not finite");
      assert( is_finite(F(2,0)) || !"compute_elastic_forces(): F(2,0) is not finite");
      assert( is_finite(F(2,1)) || !"compute_elastic_forces(): F(2,1) is not finite");
      assert( is_finite(F(2,2)) || !"compute_elastic_forces(): F(2,2) is not finite");

      //--- compute stress tensor ----------------------------------------------
      M const sigma = model->sigma(F);

      assert( is_number(sigma(0,0)) || !"compute_elastic_forces(): sigma(0,0) is not a number");
      assert( is_number(sigma(0,1)) || !"compute_elastic_forces(): sigma(0,1) is not a number");
      assert( is_number(sigma(0,2)) || !"compute_elastic_forces(): sigma(0,2) is not a number");
      assert( is_number(sigma(1,0)) || !"compute_elastic_forces(): sigma(1,0) is not a number");
      assert( is_number(sigma(1,1)) || !"compute_elastic_forces(): sigma(1,1) is not a number");
      assert( is_number(sigma(1,2)) || !"compute_elastic_forces(): sigma(1,2) is not a number");
      assert( is_number(sigma(2,0)) || !"compute_elastic_forces(): sigma(2,0) is not a number");
      assert( is_number(sigma(2,1)) || !"compute_elastic_forces(): sigma(2,1) is not a number");
      assert( is_number(sigma(2,2)) || !"compute_elastic_forces(): sigma(2,2) is not a number");

      assert( is_finite(sigma(0,0)) || !"compute_elastic_forces(): sigma(0,0) is not finite");
      assert( is_finite(sigma(0,1)) || !"compute_elastic_forces(): sigma(0,1) is not finite");
      assert( is_finite(sigma(0,2)) || !"compute_elastic_forces(): sigma(0,2) is not finite");
      assert( is_finite(sigma(1,0)) || !"compute_elastic_forces(): sigma(1,0) is not finite");
      assert( is_finite(sigma(1,1)) || !"compute_elastic_forces(): sigma(1,1) is not finite");
      assert( is_finite(sigma(1,2)) || !"compute_elastic_forces(): sigma(1,2) is not finite");
      assert( is_finite(sigma(2,0)) || !"compute_elastic_forces(): sigma(2,0) is not finite");
      assert( is_finite(sigma(2,1)) || !"compute_elastic_forces(): sigma(2,1) is not finite");
      assert( is_finite(sigma(2,2)) || !"compute_elastic_forces(): sigma(2,2) is not finite");

      //--- Compute 6 times the volume of the tetrahedron ----------------------
      // T const vol6    = dot( (u_mi , cross( u_ji, u_ki ) );
      //--- compute the spatial gradient of the linear shape functions (volume -
      //--- weighted coordiantes)
      //---
      //---  By definition the volume weighted coordinates of node m is given by
      //---
      //---     w_ijk(x)  = vol(i,j,k,x) / vol(i,j,k,m)
      //---
      //--- The gradient is then
      //---
      //---   Nabla_w_ijk(x) =   nabla_x vol(i,j,k,x) / vol(i,j,k,m)
      //---
      //---  where
      //---
      //---       nabla_x vol(i,j,k,x) =  (x_j-x_i) \times (x_k-i) / 6
      //---
      // V const nabla_w_jkm   = cross( u_mj, u_kj  ) / vol6;
      // V const nabla_w_ikm   = cross( u_ki, u_mi  ) / vol6;
      // V const nabla_w_ijm   = cross( u_mi, u_ji  ) / vol6;
      // V const nabla_w_ijk   = cross( u_ji, u_ki  ) / vol6;
      //
      //--- compute the nodal element forces -----------------------------------
      //---
      //--- f_i = int_V sigma nabla_w_jkm(x) dV = sigma * nabla_w_jkm * V
      //---
      //--- which can be simplied to -------------------------------------------
      //---
      //--- f_i = sigma cross( u_mj, u_kj  ) / 6
      //---
      //--- Similar for nodes j, k and m.
      //---
      V const u_mj = xm - xj;
      V const u_kj = xk - xj;

      assert( is_number(u_mj(0)) || !"compute_elastic_forces(): u_mj(0) is not a number");
      assert( is_number(u_mj(1)) || !"compute_elastic_forces(): u_mj(1) is not a number");
      assert( is_number(u_mj(2)) || !"compute_elastic_forces(): u_mj(2) is not a number");

      assert( is_finite(u_mj(0)) || !"compute_elastic_forces(): u_mj(0) is not finite"  );
      assert( is_finite(u_mj(1)) || !"compute_elastic_forces(): u_mj(1) is not finite"  );
      assert( is_finite(u_mj(2)) || !"compute_elastic_forces(): u_mj(2) is not finite"  );

      assert( is_number(u_kj(0)) || !"compute_elastic_forces(): u_kj(0) is not a number");
      assert( is_number(u_kj(1)) || !"compute_elastic_forces(): u_kj(1) is not a number");
      assert( is_number(u_kj(2)) || !"compute_elastic_forces(): u_kj(2) is not a number");

      assert( is_finite(u_kj(0)) || !"compute_elastic_forces(): u_kj(0) is not finite"  );
      assert( is_finite(u_kj(1)) || !"compute_elastic_forces(): u_kj(1) is not finite"  );
      assert( is_finite(u_kj(2)) || !"compute_elastic_forces(): u_kj(2) is not finite"  );

      V const fi = prod( sigma,  cross(u_mj, u_kj) ) / VT::numeric_cast(6.0);
      V const fj = prod( sigma,  cross(u_ki, u_mi) ) / VT::numeric_cast(6.0);
      V const fk = prod( sigma,  cross(u_mi, u_ji) ) / VT::numeric_cast(6.0);
      V const fm = prod( sigma,  cross(u_ji, u_ki) ) / VT::numeric_cast(6.0);

      assert( is_number(fi(0)) || !"compute_elastic_forces(): fi(0) is not a number");
      assert( is_number(fi(1)) || !"compute_elastic_forces(): fi(1) is not a number");
      assert( is_number(fi(2)) || !"compute_elastic_forces(): fi(2) is not a number");

      assert( is_finite(fi(0)) || !"compute_elastic_forces(): fi(0) is not finite"  );
      assert( is_finite(fi(1)) || !"compute_elastic_forces(): fi(1) is not finite"  );
      assert( is_finite(fi(2)) || !"compute_elastic_forces(): fi(2) is not finite"  );

      assert( is_number(fj(0)) || !"compute_elastic_forces(): fj(0) is not a number");
      assert( is_number(fj(1)) || !"compute_elastic_forces(): fj(1) is not a number");
      assert( is_number(fj(2)) || !"compute_elastic_forces(): fj(2) is not a number");

      assert( is_finite(fj(0)) || !"compute_elastic_forces(): fj(0) is not finite"  );
      assert( is_finite(fj(1)) || !"compute_elastic_forces(): fj(1) is not finite"  );
      assert( is_finite(fj(2)) || !"compute_elastic_forces(): fj(2) is not finite"  );

      assert( is_number(fk(0)) || !"compute_elastic_forces(): fk(0) is not a number");
      assert( is_number(fk(1)) || !"compute_elastic_forces(): fk(1) is not a number");
      assert( is_number(fk(2)) || !"compute_elastic_forces(): fk(2) is not a number");

      assert( is_finite(fk(0)) || !"compute_elastic_forces(): fk(0) is not finite"  );
      assert( is_finite(fk(1)) || !"compute_elastic_forces(): fk(1) is not finite"  );
      assert( is_finite(fk(2)) || !"compute_elastic_forces(): fk(2) is not finite"  );

      assert( is_number(fm(0)) || !"compute_elastic_forces(): fm(0) is not a number");
      assert( is_number(fm(1)) || !"compute_elastic_forces(): fm(1) is not a number");
      assert( is_number(fm(2)) || !"compute_elastic_forces(): fm(2) is not a number");

      assert( is_finite(fm(0)) || !"compute_elastic_forces(): fm(0) is not finite"  );
      assert( is_finite(fm(1)) || !"compute_elastic_forces(): fm(1) is not finite"  );
      assert( is_finite(fm(2)) || !"compute_elastic_forces(): fm(2) is not finite"  );

      local_forces[ idx*4u + 0u ] = fi;
      local_forces[ idx*4u + 1u ] = fj;
      local_forces[ idx*4u + 2u ] = fk;
      local_forces[ idx*4u + 3u ] = fm;
    }

    //--- Gather tetrahedral elements forces per vertex ------------------------
    f.resize( mesh.vertex_size() );

    for(size_t i = 0u; i < mesh.vertex_size(); ++i)  // 2014-01-22 Kenny: This outer for-loop is naive parallel
    {
      V sum = V::zero();

      for (unsigned int entry = neighbors.m_offset[i]; entry < neighbors.m_offset[i+1];++entry)
      {
        unsigned int            const   idx         = neighbors.m_V2T[ entry ].second;
        mesh_array::Tetrahedron const & tetrahedron = mesh.tetrahedron(idx);
        unsigned int            const   local_idx   = tetrahedron.get_local_index(i);

        sum += local_forces[ idx*4u + local_idx ];
      }

      f[i] = MT::convert(sum );
    }
    
  }
  
}// namespace hyper

// HYPER_COMPUTE_ELASTIC_FORCES_H
#endif
