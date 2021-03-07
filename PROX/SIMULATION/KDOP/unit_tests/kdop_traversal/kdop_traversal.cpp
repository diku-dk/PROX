#include <kdop_tandem_traversal.h>
#include <kdop_refit_tree.h>
#include <kdop_mesh_reorder.h>
#include <kdop_make_tree.h>
#include <kdop_io.h>

#include <mesh_array.h>
#include <tiny.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

typedef tiny::MathTypes<float> MT;
typedef MT::vector3_type       V;
typedef MT::real_type          T;

class GeometryInfo
{
public:

  kdop::Tree<T,6>                                   m_tree;
  mesh_array::T4Mesh                                m_mesh;
  mesh_array::VertexAttribute<T,mesh_array::T4Mesh> m_X;
  mesh_array::VertexAttribute<T,mesh_array::T4Mesh> m_Y;
  mesh_array::VertexAttribute<T,mesh_array::T4Mesh> m_Z;

  mesh_array::TetrahedronAttribute<mesh_array::TetrahedronSurfaceInfo,mesh_array::T4Mesh> m_surface_map;

};

class TestCallback : public geometry::ContactsCallback<V>
{
public:
  
  void operator()( V const & p, V const & n, V::real_type const & d )
  {
    // Do something?
  }
};

/**
 * Creates a test mesh which is basically a tetrahedron and its x-y plane mirrored counter part.
 *
 * @param p   Initial translation of the test mesh
 */
void make_test_mesh(
                      mesh_array::T3Mesh & mesh
                      , mesh_array::VertexAttribute<T,mesh_array::T3Mesh> & X
                      , mesh_array::VertexAttribute<T,mesh_array::T3Mesh> & Y
                      , mesh_array::VertexAttribute<T,mesh_array::T3Mesh> & Z
                      , V const & p
                      )
{
  mesh.clear();

  mesh.set_capacity( 5u, 6u );

  X.bind(mesh);
  Y.bind(mesh);
  Z.bind(mesh);

  mesh_array::Vertex const v1 = mesh.push_vertex();
  mesh_array::Vertex const v2 = mesh.push_vertex();
  mesh_array::Vertex const v3 = mesh.push_vertex();
  mesh_array::Vertex const vp = mesh.push_vertex();
  mesh_array::Vertex const vm = mesh.push_vertex();

  X(v1) = -1.0f + p(0);    Y(v1) = -1.0f + p(1);    Z(v1) =  0.0f + p(2);
  X(v2) =  1.0f + p(0);    Y(v2) = -1.0f + p(1);    Z(v2) =  0.0f + p(2);
  X(v3) =  0.0f + p(0);    Y(v3) =  1.0f + p(1);    Z(v3) =  0.0f + p(2);
  X(vp) =  0.0f + p(0);    Y(vp) =  0.0f + p(1);    Z(vp) =  2.0f + p(2);
  X(vm) =  0.0f + p(0);    Y(vm) =  0.0f + p(1);    Z(vm) = -2.0f + p(2);

  mesh.push_triangle( v1, v2, vp );
  mesh.push_triangle( v2, v3, vp );
  mesh.push_triangle( v3, v1, vp );

  mesh.push_triangle( v2, v1, vm );
  mesh.push_triangle( v3, v2, vm );
  mesh.push_triangle( v1, v3, vm );
}


void make_geometry( GeometryInfo & info, V const & p )
{
  mesh_array::T3Mesh surface;
  mesh_array::VertexAttribute<T,mesh_array::T3Mesh> sX;
  mesh_array::VertexAttribute<T,mesh_array::T3Mesh> sY;
  mesh_array::VertexAttribute<T,mesh_array::T3Mesh> sZ;

  make_test_mesh(surface,sX,sY,sZ, p);

  mesh_array::T4Mesh mesh_in;
  mesh_array::VertexAttribute<T,mesh_array::T4Mesh> X_in;
  mesh_array::VertexAttribute<T,mesh_array::T4Mesh> Y_in;
  mesh_array::VertexAttribute<T,mesh_array::T4Mesh> Z_in;
  mesh_array::tetgen(surface, sX, sY, sZ, mesh_in, X_in, Y_in, Z_in);
  
  kdop::mesh_reorder( mesh_in, X_in, Y_in, Z_in, info.m_mesh, info.m_X, info.m_Y, info.m_Z );

  mesh_array::compute_surface_map( info.m_mesh, info.m_X, info.m_Y, info.m_Z, info.m_surface_map );

  info.m_tree = kdop::make_tree<V,6,T>( 32000, info.m_mesh, info.m_X, info.m_Y, info.m_Z );

  BOOST_CHECK( info.m_surface_map.size()>0u );

}

void make_geometry( GeometryInfo & info )
{
  make_geometry(info, V::zero() );
}


BOOST_AUTO_TEST_SUITE(kdop);

BOOST_AUTO_TEST_CASE(kdop_tandem_traversal)
{
  using namespace kdop::io;
  
  GeometryInfo A;
  GeometryInfo B;

  //
  // After building 6DOP tree there should be a total of 3 nodes. the
  // root is AABB( (-1,-1,-2), (1,1,2) ) and the root should have two
  // children nodes AABB( (-1,-1,-2), (1,1,0) ) and AABB( (-1,-1,0), (1,1,2) )
  // and both children are leafs
  //
  make_geometry(A, V::make(0.0f,0.0f,-1.5f));
  make_geometry(B, V::make(0.0f,0.0f, 1.5f));

  TestCallback my_test_callback;

  std::vector<kdop::TestPair<V, 6, T> > test_pairs;

  kdop::TestPair<V, 6, T> test =  kdop::TestPair<V, 6, T> (
                                                           A.m_tree
                                                           , B.m_tree
                                                           , A.m_mesh
                                                           , B.m_mesh
                                                           , A.m_X
                                                           , B.m_X
                                                           , A.m_Y
                                                           , B.m_Y
                                                           , A.m_Z
                                                           , B.m_Z
                                                           , A.m_surface_map
                                                           , B.m_surface_map
                                                           , my_test_callback
                                                           );

  test_pairs.push_back(test);

  kdop::tandem_traversal<V,6,T>( test_pairs );
}

BOOST_AUTO_TEST_SUITE_END();
