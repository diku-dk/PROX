#ifndef HYPER_COMPUTE_CFL_TIME_STEP_SIZE_H
#define HYPER_COMPUTE_CFL_TIME_STEP_SIZE_H

#include <hyper_engine.h>
#include <hyper_math_policy.h>

#include <mesh_array.h>
#include <tiny_is_finite.h>
#include <tiny_is_number.h>

#include <cassert>

namespace hyper
{

  /**
   * Compute CFL time-step size.
   *
   * @tparam MT    Math types type binder.
   */
  template<typename MT>
  inline typename MT::real_type compute_CFL_time_step_size(
                                                      mesh_array::T4Mesh  const & mesh
                                                      , typename MT::vector_block3x1_type const & x
                                                      , typename MT::vector_block3x1_type const & v
                                                      , typename MT::real_type const & dt_wanted
                                                      )
  {
    typedef typename MT::vector3_type            V;
    typedef typename MT::real_type               T;
    typedef typename MT::value_traits            VT;

    unsigned int const K = mesh.tetrahedron_size();

    T dt = dt_wanted;

    //--- Compute local elastic forces per tetrahedral element -----------------
    for ( unsigned int idx = 0u; idx < K; ++idx)  // 2014-01-22 Kenny: This for-loop is naive parallel
    {
      mesh_array::Tetrahedron const & tetrahedron = mesh.tetrahedron(idx);

      unsigned int const i  = tetrahedron.i();
      unsigned int const j  = tetrahedron.j();
      unsigned int const k  = tetrahedron.k();
      unsigned int const m  = tetrahedron.m();

      V const xi = MT::convert( x[i] );
      V const xj = MT::convert( x[j] );
      V const xk = MT::convert( x[k] );
      V const xm = MT::convert( x[m] );

      V const vi = MT::convert( v[i] );
      V const vj = MT::convert( v[j] );
      V const vk = MT::convert( v[k] );
      V const vm = MT::convert( v[m] );

      V u_ji  = xj - xi;
      V u_ki  = xk - xi;
      V u_mi  = xm - xi;

      T const vol6  = inner_prod( u_mi,  cross(u_ji, u_ki)   );

      assert( vol6 > VT::zero() || !"compute_CFL_time_step(): initial degenerate tetrahedron");

      T dt_searh = dt;

      u_ji = (xj + dt_searh*vj) - (xi + dt_searh*vi);
      u_ki = (xk + dt_searh*vk) - (xi + dt_searh*vi);
      u_mi = (xm + dt_searh*vm) - (xi + dt_searh*vi);

      T vol6_search = inner_prod( u_mi, cross(u_ji, u_ki));  // Objective function, find small enought dt such that volume is scritly positive

      T const alpha       = VT::numeric_cast(0.63);   // Step reduction parameter
      T const beta        = VT::numeric_cast(0.9);    // Sufficient decrease paramter

      //--- We could use vol6_serach < VT::zero(), but to make deformation
      //--- rate limited numericall we use a relaxed version, requiring only
      //-- a fractional allowed change of the original positive volume
      while( vol6_search < vol6*beta  ) // sufficient decrease conditions
      {
        dt_searh = alpha*dt_searh;  // Step reduction step

        u_ji = (xj + dt_searh*vj) - (xi + dt_searh*vi);
        u_ki = (xk + dt_searh*vk) - (xi + dt_searh*vi);
        u_mi = (xm + dt_searh*vm) - (xi + dt_searh*vi);

        vol6_search = inner_prod( u_mi, cross(u_ji, u_ki));  // update objective
      }

      dt = dt_searh;

    }

    return dt;
    
  }



  /**
   * Use a CFL condition to find a "safe" time step that
   * will not cause immediate degenerate tetrahedra
   * in the mesh.
   */
  template<typename MT>
  inline typename MT::real_type compute_CFL_time_step_size(
                                                           Engine<MT> & engine
                                                           , typename MT::real_type const & time_wanted
                                                           )
  {
    using std::min;

    typedef typename MT::real_type       T;

    typedef          Engine<MT>                            engine_type;
    typedef typename engine_type::body_iterator            body_iterator;
    typedef typename engine_type::const_body_iterator      const_body_iterator;

    T time_safe = time_wanted;

    for (body_iterator body = engine.body_begin();body != engine.body_end(); ++body)
    {
      T const time_cfl = compute_CFL_time_step_size<MT>(
                                               body->m_mesh
                                               , body->m_x
                                               , body->m_v
                                               , time_wanted
                                               );

      time_safe = min ( time_safe, time_cfl);
    }

    return time_safe;
  }


}// namespace hyper

// HYPER_COMPUTE_CFL_TIME_STEP_SIZE_H
#endif
