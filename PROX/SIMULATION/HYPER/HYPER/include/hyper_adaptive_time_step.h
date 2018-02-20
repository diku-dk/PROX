#ifndef HYPER_ADAPTIVE_TIME_STEP_H
#define HYPER_ADAPTIVE_TIME_STEP_H

#include <hyper_body.h>
#include <hyper_engine.h>
#include <hyper_semi_implicit_time_step.h>

#include <util_log.h>

#include <algorithm>

namespace hyper
{

  template< typename MT >
  class State
  {
  public:

    typedef typename MT::vector_block3x1_type  vector_block3x1_type;

  public:

    vector_block3x1_type  m_x;
    vector_block3x1_type  m_v;

  private:

    void copy(State const & state)
    {
      if(this == &state)
        return;

      this->m_x = state.m_x;
      this->m_v = state.m_v;
    }

  public:

    State()
    : m_x()
    , m_v()
    {
      this->clear();
    }

    virtual ~State()
    {
      this->clear();
    }

    State(State const & state)
    : m_x()
    , m_v()
    {
      this->copy(state);
    }

    State & operator= (State const & state)
    {
      this->copy(state);
      return *this;
    }

  public:

    void clear()
    {
      this->m_x.clear();
      this->m_v.clear();
    }

  };

  template< typename MT >
  inline void copy(State<MT> const & state, Body<MT> & body )
  {
    body.m_x = state.m_x;
    body.m_v = state.m_v;
  }

  template< typename MT >
  inline void copy(Body<MT> const & body, State<MT> & state )
  {
    state.m_x = body.m_x;
    state.m_v = body.m_v;
  }

  template< typename MT >
  inline void do_time_step(Body<MT> const & body, State<MT>  & state, typename MT::real_type const & dt )
  {
    semi_implicit_time_step(
                            body.m_mesh
                            , body.m_neighbors
                            , body.m_dirichlet_conditions
                            , dt
                            , state.m_x
                            , body.m_x0
                            , state.m_v
                            , body.m_F
                            , body.m_model
                            , jacobi_precond
                            );
  }

  template< typename MT >
  inline typename MT::real_type accuracy(State<MT> const & A, State<MT> const & B)
  {
    typedef typename MT::real_type             T;
    typedef typename MT::value_traits          VT;
    typedef typename MT::vector_block3x1_type  vector_block3x1_type;

    using std::abs;
    using std::max;

    vector_block3x1_type y;
    sparse::sub(A.m_x,B.m_x, y);

    unsigned int const N = y.size();

    T error = VT::zero();

    for(unsigned int i = 0u; i < N; ++i)
    {
      error = max( abs( y(i)(0) ), error )  ;
      error = max( abs( y(i)(1) ), error )  ;
      error = max( abs( y(i)(2) ), error )  ;
    }

    return error;
  }


  template< typename MT >
  inline void adaptive_time_step(
                                   typename MT::real_type const & dt
                                 , Body<MT> & body
                                 , Params<MT> const & params
                                 )
  {
    typedef typename MT::real_type             T;
    typedef typename MT::value_traits          VT;

    using std::min;
    using std::max;


    T dt_left = dt;

    // If adaptive time stepper is not initialized pick some sensible time-step value
    if (body.m_adaptive_dt <= VT::zero() )
    {
      body.m_adaptive_dt = dt;
      body.m_adaptive_dt = max( params.adaptive_min_dt(), min( params.adaptive_max_dt(), body.m_adaptive_dt) ); // Clamp between min and max allowed time step
      body.m_adaptive_dt = min( body.m_adaptive_dt, dt_left);                                 // Make sure we do not step longer than needed
    }

    while (dt_left> VT::zero() )
    {
      // First we try if we are allowed to increase the time-step size
      if( body.m_adaptive_unchanged > params.adaptive_doubling_count() )
      {
        body.m_adaptive_dt = body.m_adaptive_dt*VT::two();

        body.m_adaptive_dt = max( params.adaptive_min_dt(), min( params.adaptive_max_dt(), body.m_adaptive_dt) ); // Clamp between min and max allowed time step
        body.m_adaptive_dt = min( body.m_adaptive_dt, dt_left);                                 // Make sure we do not step longer than needed

        body.m_adaptive_unchanged = 0u;
      }

      // Then we make a test time-integration
      State<MT> A;
      State<MT> B;

      copy(body,A);
      copy(body,B);

      do_time_step(body, A, (body.m_adaptive_dt*VT::half())  );
      do_time_step(body, A, (body.m_adaptive_dt*VT::half())  );
      do_time_step(body, B, body.m_adaptive_dt               );

      // If accuracy is acceptable we accept the time-step or if it is impossible to do schrinking to get better accuracy
      if( accuracy(A,B) <= params.adaptive_halving_tolerance()  || dt_left <= params.adaptive_min_dt() )
      {
        copy(A,body);

        ++body.m_adaptive_unchanged;

      }else{

        // While accuracy is unacceptable we keep on schrinking the step-size to get better accuracy
        do{

          body.m_adaptive_dt = body.m_adaptive_dt*VT::half();

          body.m_adaptive_dt = max( params.adaptive_min_dt(), min( params.adaptive_max_dt(), body.m_adaptive_dt) ); // Clamp between min and max allowed time step
          body.m_adaptive_dt = min( body.m_adaptive_dt, dt_left);                                 // Make sure we do not step longer than needed

          copy(body,A);
          copy(body,B);

          do_time_step(body, A, body.m_adaptive_dt*VT::half()  );
          do_time_step(body, A, body.m_adaptive_dt*VT::half()  );
          do_time_step(body, B, body.m_adaptive_dt             );

          // Impossible to schrink more, we have to give up and accept the accuracy
          if( body.m_adaptive_dt == params.adaptive_min_dt())
            break;
          
        }while( accuracy(A,B) > params.adaptive_halving_tolerance() );

        copy(A,body);

        body.m_adaptive_unchanged = 1u;
      }

      util::Log logging;
      logging << "hyper::adaptive_time_step(): " << "Adaptive: accepted dt = " << body.m_adaptive_dt << util::Log::newline();
      
      dt_left -= body.m_adaptive_dt;
    }
  }

}// end of namespace hyper

// HYPER_ADAPTIVE_TIME_STEP_H
#endif
