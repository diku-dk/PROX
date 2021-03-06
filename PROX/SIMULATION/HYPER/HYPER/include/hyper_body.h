#ifndef HYPER_BODY_H
#define HYPER_BODY_H

#include <hyper_params.h>
#include <hyper_constitutive_equation.h>
#include <hyper_dirichlet_info.h>
#include <hyper_traction_info.h>
#include <hyper_scripted_motion.h>

#include <tiny.h>
#include <mesh_array.h>
#include <kdop.h>

#include <cassert>
#include <vector>

namespace hyper
{

  template< typename MT >
  class Body
  {
  public:
    
    typedef typename MT::real_type             T;
    typedef typename MT::vector3_type          V;
    typedef typename MT::matrix3x3_type        M;
    typedef typename MT::quaternion_type       Q;
    typedef typename MT::value_traits          VT;
    typedef typename MT::vector_block3x1_type  vector_block3x1_type;


  protected:
    
    size_t       m_idx;           ///< A body index.
    size_t       m_visual_idx;    ///< Auxiliary index that can be used by an application to associate visualization data with this body.
    std::string  m_name;          ///< A name.

  public:

    T                   m_radius;        ///< Radius of a bounding sphere of the body.
    V                   m_center;        ///< Center of bounding sphere of the body.
    V                   m_min_point;     ///< Minimum corner point of bounding AABB of the body
    V                   m_max_point;     ///< Maximum corner point of bounding AABB of the body

    mesh_array::T4Mesh  m_mesh;      ///< Tetrahedral mesh

    mesh_array::VertexAttribute<T,mesh_array::T4Mesh> m_X0;     ///< Undeformed (material) x-coordinate
    mesh_array::VertexAttribute<T,mesh_array::T4Mesh> m_Y0;     ///< Undeformed (material) y-coordinate
    mesh_array::VertexAttribute<T,mesh_array::T4Mesh> m_Z0;     ///< Undeformed (material) z-coordinate
    mesh_array::VertexAttribute<T,mesh_array::T4Mesh> m_X;      ///< Deformed (spatial) x-coordinate
    mesh_array::VertexAttribute<T,mesh_array::T4Mesh> m_Y;      ///< Deformed (spatial) y-coordinate
    mesh_array::VertexAttribute<T,mesh_array::T4Mesh> m_Z;      ///< Deformed (spatial) z-coordinate

    mesh_array::TetrahedronAttribute<mesh_array::TetrahedronSurfaceInfo,mesh_array::T4Mesh> m_opposite_face_on_surf;

    kdop::Tree<T,8> m_tree;  ///< kDOP BVH tree used for collision detection of deformed (spatial) mesh

  public:

    mesh_array::Neighborhoods         m_neighbors;
    std::vector<DirichletInfo<MT> >   m_dirichlet_conditions;
    std::vector<TractionInfo<MT> >    m_traction_conditions;

    vector_block3x1_type                   m_x0;
    vector_block3x1_type                   m_x;
    vector_block3x1_type                   m_v;
    vector_block3x1_type                   m_Fext;   ///< External body forces
    vector_block3x1_type                   m_F;      ///< Total applied forces (extenal body forces and surface traction (load) forces)

    ConstitutiveEquation<MT>             * m_model;

  public:

    T            m_adaptive_dt;          ///< When using adaptive time-stepping
                                         ///< this holds the current time-step value
                                         ///< used to time step this body.
    unsigned int m_adaptive_unchanged;   ///< When using adaptive time-stepping
                                         ///< this member holds of consecutive
                                         ///< times in a row that the adaptive
                                         ///< time stepper used the same time-step
                                         ///< size.

  public:

    ScriptedMotion<MT>                  * m_scripted_motion;

  private:

    void copy (Body const & body)
    {
      if(this == &body)
        return;

      this->m_idx          = body.m_idx;
      this->m_visual_idx   = body.m_visual_idx;
      this->m_name         = body.m_name;

      this->m_radius       = body.m_radius;
      this->m_center       = body.m_center;
      this->m_min_point    = body.m_min_point;
      this->m_max_point    = body.m_max_point;

      this->m_mesh         = body.m_mesh;
      this->m_X0           = body.m_X0;
      this->m_Y0           = body.m_Y0;
      this->m_Z0           = body.m_Z0;
      this->m_X            = body.m_X;
      this->m_Y            = body.m_Y;
      this->m_Z            = body.m_Z;

      this->m_opposite_face_on_surf = body.m_opposite_face_on_surf;

      this->m_tree         = body.m_tree;

      this->m_neighbors            = body.m_neighbors;
      this->m_dirichlet_conditions = body.m_dirichlet_conditions;
      this->m_traction_conditions  = body.m_traction_conditions;
      this->m_x0                   = body.m_x0;
      this->m_x                    = body.m_x;
      this->m_v                    = body.m_v;
      this->m_Fext                 = body.m_Fext;
      this->m_F                    = body.m_F;
      this->m_model                = body.m_model;
      this->m_adaptive_dt          = body.m_adaptive_dt;
      this->m_adaptive_unchanged   = body.m_adaptive_unchanged;
      this->m_scripted_motion      = body.m_scripted_motion;
    }
    
  public:
    
    Body()
    : m_idx()
    , m_visual_idx()
    , m_name()
    , m_radius()
    , m_center()
    , m_min_point()
    , m_max_point()
    , m_mesh()
    , m_X0()
    , m_Y0()
    , m_Z0()
    , m_X()
    , m_Y()
    , m_Z()
    , m_opposite_face_on_surf()
    , m_tree()
    , m_neighbors()
    , m_dirichlet_conditions()
    , m_traction_conditions()
    , m_x0()
    , m_x()
    , m_v()
    , m_Fext()
    , m_F()
    , m_model(0)
    , m_adaptive_dt()
    , m_adaptive_unchanged()
    , m_scripted_motion(0)
    {
      this->clear();
    }
    
    virtual ~Body()
    {
      this->clear();
    }
    
    Body (Body const & body)
    : m_idx()
    , m_visual_idx()
    , m_name()
    , m_radius()
    , m_center()
    , m_min_point()
    , m_max_point()
    , m_mesh()
    , m_X0()
    , m_Y0()
    , m_Z0()
    , m_X()
    , m_Y()
    , m_Z()
    , m_opposite_face_on_surf()
    , m_tree()
    , m_neighbors()
    , m_dirichlet_conditions()
    , m_traction_conditions()
    , m_x0()
    , m_x()
    , m_v()
    , m_Fext()
    , m_F()
    , m_model(0)
    , m_adaptive_dt()
    , m_adaptive_unchanged()
    , m_scripted_motion(0)
    {
      this->copy(body);
    }
    
    Body & operator= (Body const & body)
    {
      this->copy(body);
      return *this;
    }
    
  public:

    void init(
              Params<MT> const & params
              , mesh_array::T4Mesh const & mesh
              , mesh_array::VertexAttribute<T,mesh_array::T4Mesh> const & X
              , mesh_array::VertexAttribute<T,mesh_array::T4Mesh> const & Y
              , mesh_array::VertexAttribute<T,mesh_array::T4Mesh> const & Z
              )
    {
      using tiny::norm;

      kdop::mesh_reorder( mesh
                         , X
                         , Y
                         , Z
                         , this->m_mesh
                         , this->m_X0
                         , this->m_Y0
                         , this->m_Z0
                         );

      mesh_array::compute_surface_map( this->m_mesh
                                , this->m_X0
                                , this->m_Y0
                                , this->m_Z0
                                , this->m_opposite_face_on_surf
                                );


      size_t mem_bytes = 8000u;  // 8000 bytes is default chunk size

      this->m_tree = kdop::make_tree<V,8,T>(
                                             mem_bytes
                                             , this->m_mesh
                                             , this->m_X0
                                             , this->m_Y0
                                             , this->m_Z0
                                             , kdop::sequential()
                                             );

      this->m_X.bind(this->m_mesh);
      this->m_Y.bind(this->m_mesh);
      this->m_Z.bind(this->m_mesh);

      this->m_X = this->m_X0;
      this->m_Y = this->m_Y0;
      this->m_Z = this->m_Z0;

      mesh_array::compute_neighbors(this->m_mesh, this->m_neighbors);

      this->m_x0 = MT::convert( this->m_X0, this->m_Y0, this->m_Z0);
      this->m_x  = MT::convert( this->m_X, this->m_Y, this->m_Z);

      this->m_v.resize(this->m_x0.size() );
      this->m_v.clear_data();

      this->m_Fext.resize(this->m_x0.size() );
      this->m_Fext.clear_data();

      this->m_F.resize(this->m_x0.size() );
      this->m_F.clear_data();
    }

    bool empty() const
    {
      return ((this->m_mesh.vertex_size() == 0u) || (this->m_mesh.tetrahedron_size() == 0u));
    }

    void clear()
    {
      this->m_idx = 0u;
      this->m_visual_idx = 0u;
      this->m_name = "";

      this->m_radius = VT::one();
      this->m_center.clear();
      this->m_min_point.clear();
      this->m_max_point.clear();

      this->m_mesh.clear();

      this->m_X.release();
      this->m_Y.release();
      this->m_Z.release();

      this->m_X0.release();
      this->m_Y0.release();
      this->m_Z0.release();
      this->m_opposite_face_on_surf.release();

      this->m_tree.clear();

      this->m_neighbors.clear();
      this->m_dirichlet_conditions.clear();
      this->m_traction_conditions.clear();
      this->m_x0.clear();
      this->m_x.clear();
      this->m_v.clear();
      this->m_Fext.clear();
      this->m_F.clear();
      this->m_model = 0;
      this->m_scripted_motion = 0;

      this->m_adaptive_dt = VT::zero();
      this->m_adaptive_unchanged = 0u;
    }
    
    void set_idx( size_t const & idx )
    {
      this->m_idx = idx;
    }

    size_t const & get_idx() const
    {
      return this->m_idx;
    }

    void set_visual_idx( size_t const & visual_idx )
    {
      this->m_visual_idx = visual_idx;
    }

    size_t const & get_visual_idx() const
    {
      return this->m_visual_idx;
    }

    void set_name(std::string const & name)
    {
      this->m_name = name;
    }

    std::string const & get_name() const
    {
      return this->m_name;
    }

  };

} // namespace hyper

// HYPER_BODY_H
#endif 
