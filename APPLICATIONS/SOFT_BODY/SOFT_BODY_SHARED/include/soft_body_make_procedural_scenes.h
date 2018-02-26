#ifndef SOFT_BODY_MAKE_PROCEDURAL_SCENES_H
#define SOFT_BODY_MAKE_PROCEDURAL_SCENES_H

#include <hyper.h>

#include <tiny.h>

#include <util_string_helper.h>
#include <util_config_file.h>
#include <util_log.h>

namespace soft_body
{

  template<typename MT>
  inline void make_base_scene(
                                hyper::Engine<MT> & engine
                                , util::ConfigFile const & config_file
                                , size_t & vid0
                                , size_t & vid1
                                , size_t & vid2
                                )
  {
    typedef typename MT::real_type          T;
    typedef typename MT::vector3_type       V;
    typedef typename MT::matrix3x3_type     M;
    typedef typename MT::quaternion_type    Q;
    typedef typename MT::coordsys_type      C;
    typedef typename MT::value_traits       VT;

    util::Log logging;

    //------------------------------------------------------------------------
    //--- First we create visual materials
    //------------------------------------------------------------------------
    vid0 = engine.create_visual_material();
    vid1 = engine.create_visual_material();
    vid2 = engine.create_visual_material();

    engine.set_visual_material_name(vid0, "Stone");
    engine.set_visual_material_name(vid1, "Ground");
    engine.set_visual_material_name(vid2, "Cannonball");

    //------------------------------------------------------------------------
    //--- Second we create "world" (global) stuff, like simulation
    //--- parameters such as gravity and materials etc..
    //------------------------------------------------------------------------
    T const gravity_x = util::to_value<T>( config_file.get_value("gravity_x", "0.0")   );
    T const gravity_y = util::to_value<T>( config_file.get_value("gravity_y", "-9.82") );
    T const gravity_z = util::to_value<T>( config_file.get_value("gravity_z", "0.0")   );

    V const gravity = V::make(gravity_x, gravity_y, gravity_z );

    engine.gravity() = gravity;

    logging << "soft_body::make_base_scene(): Using gravity = " << gravity << util::Log::newline();

    std::string const model_type = config_file.get_value("material_model", "saint_vernant_kirchoff");

    logging << "soft_body::make_base_scene(): Using "<< model_type << " model" << util::Log::newline();

    engine.set_default_model_type(model_type);

    std::string const material_name = config_file.get_value("material_name", "default");

    if( material_name.compare("custom")==0)
    {
      logging << "soft_body::make_base_scene(): Using custom material" << util::Log::newline();

      bool const use_lumped  = util::to_value<bool>( config_file.get_value("use_lumped", "true")     );
      T    const E           = util::to_value<T>( config_file.get_value("young_modulus", "1000.0")   );
      T    const nu          = util::to_value<T>( config_file.get_value("poisson_ratio", "0.3")      );
      T    const rho         = util::to_value<T>( config_file.get_value("mass_density", "10.0")      );
      T    const alpha       = util::to_value<T>( config_file.get_value("mass_damping", "0.0")       );
      T    const beta        = util::to_value<T>( config_file.get_value("stiffness_damping", "0.0")  );
      T    const c           = util::to_value<T>( config_file.get_value("viscous_damping", "0.0001") );

      engine.default_model()->material().lumped() = use_lumped;
      engine.default_model()->material().E()      = E;
      engine.default_model()->material().nu()     = nu;
      engine.default_model()->material().rho()    = rho;
      engine.default_model()->material().alpha()  = alpha;
      engine.default_model()->material().beta()   = beta;
      engine.default_model()->material().c()      = c;

    }
    else
    {
      logging << "soft_body::make_base_scene(): Using" << material_name << "material" << util::Log::newline();

      engine.default_model()->material() = hyper::make_material<T>(material_name);
    }
  }

  template<typename MT>
  inline void make_default_scene( hyper::Engine<MT> & engine, util::ConfigFile const & config_file )
  {
    typedef typename MT::real_type          T;
    typedef typename MT::vector3_type       V;
    typedef typename MT::matrix3x3_type     M;
    typedef typename MT::quaternion_type    Q;
    typedef typename MT::coordsys_type      C;
    typedef typename MT::value_traits       VT;

    size_t vid0;
    size_t vid1;
    size_t vid2;

    make_base_scene(engine, config_file, vid0, vid1, vid2);

    //------------------------------------------------------------------------
    //--- Third we create the bodies that occupy the world
    //------------------------------------------------------------------------

    size_t const body_idx0 = engine.create_body();
    size_t const body_idx1 = engine.create_body();
    size_t const body_idx2 = engine.create_body();
    size_t const body_idx3 = engine.create_body();
    size_t const body_idx4 = engine.create_body();

    //------------------------------------------------------------------------
    //--- Fourth we create the geometries of the bodies in the world
    //------------------------------------------------------------------------

    engine.get_body(body_idx0).set_visual_idx( vid0 );
    engine.get_body(body_idx1).set_visual_idx( vid1 );
    engine.get_body(body_idx2).set_visual_idx( vid2 );
    engine.get_body(body_idx3).set_visual_idx( vid2 );
    engine.get_body(body_idx4).set_visual_idx( vid2 );

    engine.set_cone_geometry(body_idx0, 4.0f, 5.0f);
    engine.set_capsule_geometry(body_idx1, 2.0f, 8.0f);
    engine.set_box_geometry(body_idx2, 2.0f, 1.0f, 8.0f);
    engine.set_box_geometry(body_idx3, 2.0f, 1.0f, 8.0f);
    engine.set_box_geometry(body_idx4, 2.0f, 1.0f, 8.0f);

    scale( engine.get_body(body_idx0), V::make(1.0, 1.1, 1.0) );
    translate( engine.get_body(body_idx0), V::make(-6.0, 0.0, 1.0) );

    bend( engine.get_body(body_idx1), V::j(), V::i(), VT::pi()*0.25 );

    twist( engine.get_body(body_idx2), V::k(), VT::pi()*0.25 );
    translate( engine.get_body(body_idx2), V::make(4.0, 0.0, -1.0) );

    bend( engine.get_body(body_idx3), V::k(), V::i(), VT::pi()*0.25 );
    translate( engine.get_body(body_idx3), V::make(0.0, 10.0, 0.0) );

    //------------------------------------------------------------------------
    //--- Fifth we specify any Dirichlet or Traction forces on the bodies
    //--- in the world
    //------------------------------------------------------------------------

    hyper::make_dirichlet_conditions(
                              engine.get_body(body_idx3)
                                     , hyper::make_box<MT>(V::make(0.0, 0.0, -4.0), V::make(3.0, 2.0, 0.01))
                              , hyper::use_material_coordinates()
                              );

    hyper::make_traction_conditions(
                             engine.get_body(body_idx3)
                             , V::make(0.0,-1.0,0.0)
                             , hyper::make_box<MT>(V::make(0.0, 0.0, 4.0),V::make(3.0, 2.0, 0.01))
                             , hyper::use_material_coordinates()
                             );


    //------------------------------------------------------------------------
    //--- Sixth we specify any scripted motions and attach them to the bodies
    //--- they should act on
    //------------------------------------------------------------------------

    size_t const twist_idx0 = engine.create_twister_motion();

    engine.get_body(body_idx4).m_scripted_motion = &(engine.get_twister_motion(twist_idx0));

  }

  template<typename MT>
  inline void make_cantilever_tower_scene( hyper::Engine<MT> & engine, util::ConfigFile const & config_file )
  {
    typedef typename MT::real_type          T;
    typedef typename MT::vector3_type       V;
    typedef typename MT::matrix3x3_type     M;
    typedef typename MT::quaternion_type    Q;
    typedef typename MT::coordsys_type      C;
    typedef typename MT::value_traits       VT;

    size_t vid0;
    size_t vid1;
    size_t vid2;

    make_base_scene(engine, config_file, vid0, vid1, vid2);

    size_t const segments       = util::to_value<size_t>( config_file.get_value("segments", "3")           );
    size_t const slices         = util::to_value<size_t>( config_file.get_value("slices",   "5")           );
    T      const radius         = util::to_value<T>(      config_file.get_value("radius",       "2.5")     );
    T      const beam_depth     = util::to_value<T>(      config_file.get_value("beam_depth",   "5.0")     );
    T      const beam_height    = util::to_value<T>(      config_file.get_value("beam_height",  "1.0")     );
    bool   const flip_dirichlet = util::to_value<bool>(   config_file.get_value("flip_dirichlet", "false") );


    assert( slices      > 2u         || !"make_cantilever_tower_scene(): slices must be larger than two");
    assert( segments    > 0u         || !"make_cantilever_tower_scene(): segments must be positive");
    assert( radius      > VT::zero() || !"make_cantilever_tower_scene(): radius must be positive");
    assert( beam_depth  > VT::zero() || !"make_cantilever_tower_scene(): beam depth  must be positive");
    assert( beam_height > VT::zero() || !"make_cantilever_tower_scene(): beam height must be positive");

    T const delta_theta    = VT::two()* VT::pi()/ slices;
    T const center_radius  = radius + beam_depth*VT::half();
    T const beam_width     = VT::numeric_cast(1.9)*radius*sin( delta_theta*VT::half() );


    for(size_t i = 0u; i < segments; ++i)
    {
      for(size_t j = 0u; j < slices; ++j)
      {
        T const theta = j*delta_theta + (i%2)*( delta_theta*VT::half() );
        T const x     = center_radius * cos(  theta  );
        T const y     = (i + VT::half()) * beam_height;
        T const z     = - center_radius * sin(  theta  );

        V const disp = V::make( x, y, z );
        //Q const rot  = Q::Ru( theta - VT::pi_half(),  V::j() );

        size_t const body_idx = engine.create_body();

        engine.get_body(body_idx).set_visual_idx( vid0 );

        engine.set_box_geometry(
                                body_idx
                                , beam_width
                                , beam_height
                                , beam_depth
                                );

        rotate(
               engine.get_body(body_idx)
               , V::j()
               , (theta - VT::pi_half())
               );

        translate(
                  engine.get_body(body_idx)
                  , disp
                  );

        T const sign           = flip_dirichlet ?  -VT::one() : VT::one() ;
        V const box_origo      = V::make( VT::zero() , VT::zero(), sign*beam_depth*VT::half());
        V const box_dimensions = V::make(beam_width+VT::numeric_cast(0.01), beam_height+VT::numeric_cast(0.01), VT::numeric_cast(0.01) );

        hyper::make_dirichlet_conditions(
                                  engine.get_body(body_idx)
                                  , hyper::make_box<MT>(
                                                 box_origo
                                                 , box_dimensions
                                                 )
                                  , hyper::use_material_coordinates()
                                  );

      }
    }

  }

  template<typename MT>
  inline void make_plate_stack_scene( hyper::Engine<MT> & engine, util::ConfigFile const & config_file )
  {
    typedef typename MT::real_type          T;
    typedef typename MT::vector3_type       V;
    typedef typename MT::matrix3x3_type     M;
    typedef typename MT::quaternion_type    Q;
    typedef typename MT::coordsys_type      C;
    typedef typename MT::value_traits       VT;

    size_t vid0;
    size_t vid1;
    size_t vid2;

    make_base_scene(engine, config_file, vid0, vid1, vid2);

    size_t const layers          = util::to_value<size_t>( config_file.get_value("layers", "10")           );
    T      const plate_width     = util::to_value<T>(      config_file.get_value("plate_width",   "5.0")   );
    T      const plate_height    = util::to_value<T>(      config_file.get_value("plate_height",  "1.0")   );
    T      const plate_depth     = util::to_value<T>(      config_file.get_value("plate_depth",   "5.0")   );
    T      const angle           = util::to_value<T>(      config_file.get_value("angle",   "0.0")         );
    bool   const use_scripted    = util::to_value<bool>(   config_file.get_value("use_scripted", "false")  );

    assert( layers       > 0u         || !"make_plate_stack_scene(): layers must be positive");
    assert( plate_width  > VT::zero() || !"make_plate_stack_scene(): plate width must be positive");
    assert( plate_depth  > VT::zero() || !"make_plate_stack_scene(): plate depth  must be positive");
    assert( plate_height > VT::zero() || !"make_plate_stack_scene(): plate height must be positive");

    size_t const motion_idx0 = engine.create_twister_motion();
    //size_t const motion_idx1 = engine.create_twister_motion();
    //size_t const motion_idx2 = engine.create_twister_motion();

    T theta = VT::zero();

    for(size_t i = 0u; i < layers; ++i)
    {
      T const x     = VT::zero();
      T const y     = (i + VT::half()) * plate_height;
      T const z     = VT::zero();

      V const disp = V::make( x, y, z );

      size_t const body_idx = engine.create_body();

      engine.get_body(body_idx).set_visual_idx( vid0 );

      engine.set_box_geometry(
                              body_idx
                              , plate_width
                              , plate_height
                              , plate_depth
                              );

      rotate(
             engine.get_body(body_idx)
             , V::j()
             , theta
             );

      theta += angle;

      translate(
                engine.get_body(body_idx)
                , disp
                );



      if(use_scripted)
      {
        engine.get_body(body_idx).m_scripted_motion = &(engine.get_twister_motion(motion_idx0));
      }
      else if(i==0)
      {
        V const box_origo      = V::make( VT::zero() , -plate_height*VT::half(), VT::zero() );
        V const box_dimensions = V::make(plate_width+VT::numeric_cast(0.01), VT::numeric_cast(0.01), plate_depth+VT::numeric_cast(0.01) );

        hyper::make_dirichlet_conditions(
                                  engine.get_body(body_idx)
                                  , hyper::make_box<MT>(
                                                 box_origo
                                                 , box_dimensions
                                                 )
                                  , hyper::use_material_coordinates()
                                  );
      }
      
    }
    
  }

} // end of namespace soft_body

// SOFT_BODY_MAKE_PROCEDURAL_SCENES_H
#endif
