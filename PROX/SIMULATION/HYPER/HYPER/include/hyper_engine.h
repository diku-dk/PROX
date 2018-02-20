#ifndef HYPER_ENGINE_H
#define HYPER_ENGINE_H

#include <hyper_body.h>
#include <hyper_twister_motion.h>
#include <hyper_contact_point.h>
#include <hyper_params.h>

#include <hyper_saint_vernant_kirchhoff.h>
#include <hyper_neo_hookean.h>
#include <hyper_corotational_elasticity.h>

#include <mesh_array.h>

#include <tiny.h>
#include <util_profiling.h>
#include <util_string_helper.h>
#include <util_log.h>

#include <vector>
#include <cassert>

namespace hyper
{

  // 2015-03-25 Kenny code review: No support for  Dirichlet conditions or Traction conditions, materials is limited to one default material, and external forces are only gravity.

  template<typename MT>
  class Engine
  {
  public:

    typedef typename MT::real_type             T;
    typedef typename MT::vector3_type          V;
    typedef typename MT::matrix3x3_type        M;
    typedef typename MT::quaternion_type       Q;
    typedef typename MT::value_traits          VT;

    typedef          Body< MT >                body_type;
    typedef          ContactPoint< MT >        contact_type;
    typedef          Params< MT >              params_type;

    typedef          CorotationalElasticity< MT >  corotational_type;
    typedef          SaintVernantKirchhoff< MT >   saint_vernant_kirchoff_type;
    typedef          NeoHookean< MT >              neohookean_type;
    typedef          ConstitutiveEquation< MT >    model_type;

    typedef typename std::vector<body_type>::iterator          body_iterator;
    typedef typename std::vector<body_type>::const_iterator    const_body_iterator;
    typedef typename std::vector<contact_type>::iterator       contact_iterator;
    typedef typename std::vector<contact_type>::const_iterator const_contact_iterator;

  protected:

    std::vector<body_type>     m_bodies;
    std::vector<contact_type>  m_contacts;
    params_type                m_params;
    T                          m_time;
    V                          m_gravity;

    saint_vernant_kirchoff_type   m_saint_vernant_kirchoff_model;
    neohookean_type               m_neohookean_model;
    corotational_type             m_corotational_model;
    model_type                  * m_default_model;

    std::vector<TwisterMotion<MT> > m_twister_motions;

  public:

    Engine()
    : m_bodies()
    , m_contacts()
    , m_params()
    , m_time(VT::zero())
    , m_gravity(V::zero())
    , m_saint_vernant_kirchoff_model()
    , m_neohookean_model()
    , m_corotational_model()
    , m_default_model()
    , m_twister_motions()
    {
      this->clear();
      this->m_default_model = & (this->m_saint_vernant_kirchoff_model);
    }

    virtual ~Engine()
    {
      this->clear();
    }

  public:

    void clear()
    {
      this->m_bodies.clear();
      this->m_contacts.clear();
      this->m_params.clear();
      this->m_visual_names.clear();
      this->m_time = VT::zero();

      this->m_gravity = V::zero();

      this->m_saint_vernant_kirchoff_model = saint_vernant_kirchoff_type();
      this->m_neohookean_model             = neohookean_type();
      this->m_corotational_model             = corotational_type();
      this->m_default_model                = & (this->m_saint_vernant_kirchoff_model);

      this->m_twister_motions.clear();

      util::Profiling::reset();
    }

    bool empty() const
    {
      return this->m_bodies.empty();
    }

    size_t create_body()
    {
      body_type b;
      b.set_idx( m_bodies.size() );
      b.set_name( "hyper" + util::to_string(b.get_idx()) );

      b.m_model = this->m_default_model;  // 2015-03-25 Kenny code review: by default all bodies get the default material.

      this->m_bodies.push_back( b );
      return b.get_idx();
    }

    body_type & get_body(size_t const & body_idx)
    {
      assert( body_idx <= this->m_bodies.size()              || !"get_body(): index out of range");
      assert( this->m_bodies[body_idx].get_idx() == body_idx || !"get_body(): internal error"    );

      return this->m_bodies[body_idx];
    }

    body_type const & get_body(size_t const & body_idx) const
    {
      assert( body_idx <= this->m_bodies.size()              || !"get_body(): index out of range");
      assert( this->m_bodies[body_idx].get_idx() == body_idx || !"get_body(): internal error"    );

      return this->m_bodies[body_idx];
    }

  public:

    void set_geometry(
                      size_t const & body_idx
                      , mesh_array::T3Mesh                                 const & surface
                      , mesh_array::VertexAttribute<T, mesh_array::T3Mesh> const & surface_X
                      , mesh_array::VertexAttribute<T, mesh_array::T3Mesh> const & surface_Y
                      , mesh_array::VertexAttribute<T, mesh_array::T3Mesh> const & surface_Z
                      )
    {
      mesh_array::T4Mesh                                 volume;
      mesh_array::VertexAttribute<T, mesh_array::T4Mesh> volume_X;
      mesh_array::VertexAttribute<T, mesh_array::T4Mesh> volume_Y;
      mesh_array::VertexAttribute<T, mesh_array::T4Mesh> volume_Z;

      mesh_array::tetgen(surface
                         , surface_X
                         , surface_Y
                         , surface_Z
                         , volume
                         , volume_X
                         , volume_Y
                         , volume_Z
                         , this->m_params.tetgen_settings()
                         );

      m_bodies[body_idx].init(this->m_params
                              , volume
                              , volume_X
                              , volume_Y
                              , volume_Z
                              );
    }

    void set_box_geometry(
                          size_t const & body_idx
                          , float const & width
                          , float const & height
                          , float const & depth
                          )
    {
      mesh_array::T3Mesh                                 surface;
      mesh_array::VertexAttribute<T, mesh_array::T3Mesh> surface_X;
      mesh_array::VertexAttribute<T, mesh_array::T3Mesh> surface_Y;
      mesh_array::VertexAttribute<T, mesh_array::T3Mesh> surface_Z;

      mesh_array::make_box<typename MT::base_type>(
                               width
                               , height
                               , depth
                               , surface
                               , surface_X
                               , surface_Y
                               , surface_Z
                               );

      this->set_geometry(
                         body_idx
                         , surface
                         , surface_X
                         , surface_Y
                         , surface_Z
                         );
    }

    void set_capsule_geometry(
                              size_t const & body_idx
                              , float const & radius
                              , float const & height
                              )
    {
      mesh_array::T3Mesh                                 surface;
      mesh_array::VertexAttribute<T, mesh_array::T3Mesh> surface_X;
      mesh_array::VertexAttribute<T, mesh_array::T3Mesh> surface_Y;
      mesh_array::VertexAttribute<T, mesh_array::T3Mesh> surface_Z;

      size_t const segments = 12u;
      size_t const slices   = 12u;

      mesh_array::make_capsule<typename MT::base_type>(
                                   radius
                                   , height
                                   , slices
                                   , segments
                                   , surface
                                   , surface_X
                                   , surface_Y
                                   , surface_Z
                                   );
      this->set_geometry(
                         body_idx
                         , surface
                         , surface_X
                         , surface_Y
                         , surface_Z
                         );

    }

    void set_cone_geometry(
                           size_t const & body_idx
                           , float const & radius
                           , float const & height
                           )
    {
      mesh_array::T3Mesh                                 surface;
      mesh_array::VertexAttribute<T, mesh_array::T3Mesh> surface_X;
      mesh_array::VertexAttribute<T, mesh_array::T3Mesh> surface_Y;
      mesh_array::VertexAttribute<T, mesh_array::T3Mesh> surface_Z;

      size_t const slices = 12u;

      mesh_array::make_cone<typename MT::base_type>(
                                radius
                                , height
                                , slices
                                , surface
                                , surface_X
                                , surface_Y
                                , surface_Z
                                );
      this->set_geometry(
                         body_idx
                         , surface
                         , surface_X
                         , surface_Y
                         , surface_Z
                         );

    }

    void set_cylinder_geometry(
                               size_t const & body_idx
                               , float const & radius
                               , float const & height
                               )
    {
      mesh_array::T3Mesh                                 surface;
      mesh_array::VertexAttribute<T, mesh_array::T3Mesh> surface_X;
      mesh_array::VertexAttribute<T, mesh_array::T3Mesh> surface_Y;
      mesh_array::VertexAttribute<T, mesh_array::T3Mesh> surface_Z;

      size_t const slices = 12u;

      mesh_array::make_cylinder<typename MT::base_type>(
                                    radius
                                    , height
                                    , slices
                                    , surface
                                    , surface_X
                                    , surface_Y
                                    , surface_Z
                                    );

      this->set_geometry(
                         body_idx
                         , surface
                         , surface_X
                         , surface_Y
                         , surface_Z
                         );
    }

    void set_ellipsoid_geometry(
                                size_t const & body_idx
                                , float const & sx
                                , float const & sy
                                , float const & sz
                                )
    {
      mesh_array::T3Mesh                                 surface;
      mesh_array::VertexAttribute<T, mesh_array::T3Mesh> surface_X;
      mesh_array::VertexAttribute<T, mesh_array::T3Mesh> surface_Y;
      mesh_array::VertexAttribute<T, mesh_array::T3Mesh> surface_Z;

      size_t const segments = 12u;
      size_t const slices   = 12u;

      mesh_array::make_ellipsoid<typename MT::base_type>(
                                     sx
                                     , sy
                                     , sz
                                     , slices
                                     , segments
                                     , surface
                                     , surface_X
                                     , surface_Y
                                     , surface_Z
                                     );

      this->set_geometry(
                         body_idx
                         , surface
                         , surface_X
                         , surface_Y
                         , surface_Z
                         );
    }

    void set_sphere_geometry(
                             size_t const & body_idx
                             , float const & radius
                             )
    {
      mesh_array::T3Mesh                                 surface;
      mesh_array::VertexAttribute<T, mesh_array::T3Mesh> surface_X;
      mesh_array::VertexAttribute<T, mesh_array::T3Mesh> surface_Y;
      mesh_array::VertexAttribute<T, mesh_array::T3Mesh> surface_Z;

      size_t const segments = 12u;
      size_t const slices   = 12u;

      mesh_array::make_sphere<typename MT::base_type>(
                                  radius
                                  , slices
                                  , segments
                                  , surface
                                  , surface_X
                                  , surface_Y
                                  , surface_Z
                                  );

      this->set_geometry(
                         body_idx
                         , surface
                         , surface_X
                         , surface_Y
                         , surface_Z
                         );
    }

  public:

    std::vector<std::string>  m_visual_names;

    std::string get_visual_material_name( size_t const & visual_idx )
    {
      return this->m_visual_names[visual_idx];
    }

    void set_visual_material_name( size_t const & visual_idx, std::string const & name)
    {
      this->m_visual_names[visual_idx] = name;
    }

    size_t create_visual_material()
    {
      m_visual_names.push_back("");
      return m_visual_names.size()-1u;
    }

    size_t get_number_of_visual_materials()
    {
      return m_visual_names.size();
    }

    void get_visual_material_indices(size_t * visual_idx_array)
    {
      for(size_t i = 0u; i < this->m_visual_names.size(); ++i)
        visual_idx_array[i] = i;
    }

  public:

    size_t create_twister_motion()
    {
      TwisterMotion<MT> motion;
      this->m_twister_motions.push_back(motion);
      return this->m_twister_motions.size() - 1u;
    }

    TwisterMotion<MT> & get_twister_motion(size_t const & idx)
    {
      return this->m_twister_motions[idx];
    }

  public:

    body_iterator       body_begin()       { return this->m_bodies.begin(); }
    body_iterator       body_end()         { return this->m_bodies.end();   }
    const_body_iterator body_begin() const { return this->m_bodies.begin(); }
    const_body_iterator body_end()   const { return this->m_bodies.end();   }

    contact_iterator       contact_begin()       { return this->m_contacts.begin(); }
    contact_iterator       contact_end()         { return this->m_contacts.end();   }
    const_contact_iterator contact_begin() const { return this->m_contacts.begin(); }
    const_contact_iterator contact_end()   const { return this->m_contacts.end();   }

    params_type const & params() const { return this->m_params; }
    params_type       & params()       { return this->m_params; }

    std::vector<contact_type> const & contacts() const { return this->m_contacts; }
    std::vector<contact_type>       & contacts()       { return this->m_contacts; }

    T const & time() const { return this->m_time; }
    T       & time()       { return this->m_time; }

    V const & gravity() const { return this->m_gravity; }
    V       & gravity()       { return this->m_gravity; }

    model_type const * default_model() const { return this->m_default_model; }
    model_type       * default_model()       { return this->m_default_model; }

    /**
     * @note Changing the model type will reset material settings to default material settings.
     */
    void set_default_model_type(std::string const & type)
    {
      if(type.compare("saint_vernant_kirchoff")==0)
      {
        this->m_default_model = & (this->m_saint_vernant_kirchoff_model);
      }
      else if(type.compare("neohookean")==0)
      {
        this->m_default_model = & (this->m_neohookean_model);
      }
      else if(type.compare("corotational")==0)
      {
        this->m_default_model = & (this->m_corotational_model);
      }
      else
      {
        util::Log logging;

        logging << "hyper::Engine::set_default_model_type() : unregnozed material model "
            << type
            << " supported types are"
            << util::Log::newline();

        logging << "\t"
            << "saint_vernant_kirchoff"
            << util::Log::newline();

        logging << "\t"
            << "neohookean"
            << util::Log::newline();
      }

      this->m_default_model->material() = make_default_material<T>();
    }

  };
  
} // namespace hyper

// HYPER_ENGINE_H
#endif 
