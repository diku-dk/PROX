#ifndef PROCEDURAL_MAKE_INTERNAL_EDGE_H
#define PROCEDURAL_MAKE_INTERNAL_EDGE_H

#include <procedural_factory.h>

#include <mesh_array.h>

#include <tiny_math_types.h>

#include <vector>

namespace procedural
{


  /**
   * This factory function creates a box-shape that is tesselated
   * in such a way that if one place a box-shape on top of it then
   * it will generate many internal-edge and internal-vertex cases.
   *
   * @param scene_size    This parameter is an overall scene size
   *                      scale. The default value is usually 20.
   *
   * @param mesh          Upon return this argument holds the generated mesh.
   * @param X             Upon return this argument holds the generated
   *                      x-coordinates of the mesh vertices.
   * @param Y             Upon return this argument holds the generated
   *                      y-coordinates of the mesh vertices.
   * @param Z             Upon return this argument holds the generated
   *                      z-coordinates of the mesh vertices.
   */
  template<typename T>
  inline void make_surface_geometry(
                                    T const & scene_size
                                    , mesh_array::T3Mesh & mesh
                                    , mesh_array::VertexAttribute<T,mesh_array::T3Mesh> & X
                                    , mesh_array::VertexAttribute<T,mesh_array::T3Mesh> & Y
                                    , mesh_array::VertexAttribute<T,mesh_array::T3Mesh> & Z
                                    )
  {
    size_t const V      = 14u;
    size_t const F      = 24u;

    X.release();
    Y.release();
    Z.release();

    mesh.clear();
    mesh.set_capacity(V,F);

    X.bind(mesh);
    Y.bind(mesh);
    Z.bind(mesh);

    mesh_array::Vertex const v1      = mesh.push_vertex();
    mesh_array::Vertex const v2      = mesh.push_vertex();
    mesh_array::Vertex const v3      = mesh.push_vertex();
    mesh_array::Vertex const v4      = mesh.push_vertex();
    mesh_array::Vertex const v5      = mesh.push_vertex();
    mesh_array::Vertex const v6      = mesh.push_vertex();
    mesh_array::Vertex const v7      = mesh.push_vertex();
    mesh_array::Vertex const v8      = mesh.push_vertex();
    mesh_array::Vertex const v9      = mesh.push_vertex();
    mesh_array::Vertex const v10    = mesh.push_vertex();
    mesh_array::Vertex const v11    = mesh.push_vertex();
    mesh_array::Vertex const v12    = mesh.push_vertex();
    mesh_array::Vertex const v13    = mesh.push_vertex();
    mesh_array::Vertex const v14    = mesh.push_vertex();

    T const W = 2.0*scene_size*0.05;
    T const H = 2.0*scene_size*0.05;
    T const D = scene_size*0.05;

    X(v1)  = -W;      Y(v1)  = -H;       Z(v1)  =  D;
    X(v2)  =  0.0;    Y(v2)  = -H;       Z(v2)  =  D;
    X(v3)  =  W;      Y(v3)  = -H;       Z(v3)  =  D;
    X(v4)  =  W;      Y(v4)  =  0.0;     Z(v4)  =  D;
    X(v5)  =  W;      Y(v5)  =  H;       Z(v5)  =  D;
    X(v6)  =  0.0;    Y(v6)  =  H;       Z(v6)  =  D;
    X(v7)  = -W;      Y(v7)  =  H;       Z(v7)  =  D;
    X(v8)  = -W;      Y(v8)  =  0.0;     Z(v8)  =  D;
    X(v9)  = -W;      Y(v9)  = -H;       Z(v9)  = -D;
    X(v10) =  W;      Y(v10) = -H;       Z(v10) = -D;
    X(v11) =  W;      Y(v11) =  H;       Z(v11) = -D;
    X(v12) = -W;      Y(v12) =  H;       Z(v12) = -D;
    X(v13) = -W/2.0;  Y(v13) =  -H/2.0;  Z(v13) =  D;
    X(v14) =  W/2.0;  Y(v14) =  -H/2.0;  Z(v14) =  D;

    // top face
    mesh.push_triangle( v1, v2, v13 );
    mesh.push_triangle( v1, v13, v8 );
    mesh.push_triangle( v2, v3, v14 );
    mesh.push_triangle( v3, v4, v14 );
    mesh.push_triangle( v4, v5, v6 );
    mesh.push_triangle( v6, v7, v8 );
    mesh.push_triangle( v2, v14, v6 );
    mesh.push_triangle( v14, v4, v6 );
    mesh.push_triangle( v2, v6, v13 );
    mesh.push_triangle( v13, v6, v8 );

    // bottom face
    mesh.push_triangle( v10, v9, v12 );
    mesh.push_triangle( v10, v12, v11 );

    // front face
    mesh.push_triangle( v1, v9, v2 );
    mesh.push_triangle( v9, v10, v2 );
    mesh.push_triangle( v2, v10, v3 );

    // back face
    mesh.push_triangle( v11, v12, v6 );
    mesh.push_triangle( v11, v6, v5 );
    mesh.push_triangle( v6, v12, v7 );

    // right face
    mesh.push_triangle( v3, v10, v4 );
    mesh.push_triangle( v10, v11, v4 );
    mesh.push_triangle( v11, v5, v4 );

    // left face
    mesh.push_triangle( v12, v8, v7 );
    mesh.push_triangle( v12, v9, v8 );
    mesh.push_triangle( v9, v1, v8 );
  }



  template<typename MT>
  inline void make_internal_edge(
                                 content::API *  engine
                                 , typename MT::vector3_type const & position
                                 , typename MT::quaternion_type const & orientation
                                 , typename MT::real_type const & scene_size
                                 , MaterialInfo<typename MT::real_type> const & mat_info
                                 , mesh_array::TetGenSettings tetset = mesh_array::tetgen_default_settings()
                                 )
  {
    using std::floor;
    using std::ceil;

    typedef typename MT::real_type       T;
    typedef typename MT::vector3_type    V;
    typedef typename MT::quaternion_type Q;
    typedef typename MT::value_traits    VT;

    T      const stone_density   = get_material_density<MT>(mat_info, "Stone");
    size_t const mid             = get_material_id<MT>(mat_info, "Stone");

    mesh_array::T3Mesh mesh;
    mesh_array::VertexAttribute<T,mesh_array::T3Mesh> X;
    mesh_array::VertexAttribute<T,mesh_array::T3Mesh> Y;
    mesh_array::VertexAttribute<T,mesh_array::T3Mesh> Z;

    make_surface_geometry(scene_size, mesh, X, Y, Z);

    GeometryHandle<MT> const base = create_geometry_handle_mesh<MT>(engine, mesh, X, Y, Z, tetset);

    T const H = 0.05*scene_size;
    std::vector< V > top_vertices;
    top_vertices.resize( 8u );
    top_vertices[0] = V::make( -H,  H,  H);
    top_vertices[1] = V::make(  H,  H,  H);
    top_vertices[2] = V::make( -H,  H, -H);
    top_vertices[3] = V::make(  H,  H, -H);
    top_vertices[4] = V::make( -H, -H,  H);
    top_vertices[5] = V::make(  H, -H,  H);
    top_vertices[6] = V::make( -H, -H, -H);
    top_vertices[7] = V::make(  H, -H, -H);

    GeometryHandle<MT> const top = create_geometry_handle_convex<MT>(engine, top_vertices);

    {
      V const T_b2m = base.Tb2m();
      Q const Q_b2m = base.Qb2m();

      V const T_m2l = V::make( 0.0, 0.0, 0.0 );
      //Q const Q_m2l = Q::identity();
      Q const Q_m2l = Q::Rx( - VT::pi_half() );

      V const T_l2w = position;
      Q const Q_l2w = orientation;

      V T_b2w;
      Q Q_b2w;

      compute_body_to_world_transform<MT>(
                                          T_b2m
                                          , Q_b2m
                                          , T_m2l
                                          , Q_m2l
                                          , T_l2w
                                          , Q_l2w
                                          , T_b2w
                                          , Q_b2w
                                          );

      size_t const rid = create_rigid_body<MT>(
                                               engine
                                               , T_b2w
                                               , Q_b2w
                                               , base
                                               , mid
                                               , stone_density
                                               );
      engine->set_rigid_body_fixed(rid, true);
    }

    {
      V const T_b2m = top.Tb2m();
      Q const Q_b2m = top.Qb2m();

      V const T_m2l = V::make( 0.0, 2.0*H, 0.0 );
      Q const Q_m2l = Q::identity();

      V const T_l2w = position;
      Q const Q_l2w = orientation;

      V T_b2w;
      Q Q_b2w;

      compute_body_to_world_transform<MT>(
                                          T_b2m
                                          , Q_b2m
                                          , T_m2l
                                          , Q_m2l
                                          , T_l2w
                                          , Q_l2w
                                          , T_b2w
                                          , Q_b2w
                                          );

      size_t const rid = create_rigid_body<MT>(
                                               engine
                                               , T_b2w
                                               , Q_b2w
                                               , top
                                               , mid
                                               , stone_density
                                               );

      engine->set_rigid_body_velocity(rid, 0.0, 0.0, 0.0);
    }
    
  }
  
} //namespace procedural

// PROCEDURAL_MAKE_INTERNAL_EDGE_H
#endif
