#include <tiny.h>
#include <narrow.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

typedef tiny::MathTypes<float>   M;
typedef M::real_type             T;
typedef M::vector3_type          V;
typedef M::quaternion_type       Q; 
typedef M::coordsys_type         C;

class MyObject : public narrow::Object< M > {};

typedef narrow::Geometry< M >          geometry_type;

typedef geometry_type::box_type        box_type;
typedef geometry_type::convex_type     convex_type;
typedef geometry_type::sphere_type     sphere_type;
typedef geometry_type::tetramesh_type  tetramesh_type;



class MyCallback
  : public geometry::ContactsCallback<V>
{
public:

  bool m_hit;

  MyCallback()
  : m_hit(false)
  {}


  void operator()(
                  V const & point
                  , V const & normal
                  , V::real_type const & distance
                  )
  {
    m_hit = true;
  }


};




/**
 * Test Two Shapes Function.
 * This is a convenience function that makes it more easy to write test cases.
 */
template<typename shape_typeA, typename shape_typeB>
inline void do_test(
                    V posA
                    , V posB
                    , Q rotA
                    , Q rotB
                    , shape_typeA const & shapeA
                    , shape_typeB const & shapeB
                    , bool const & should_overlap
                    )
{
  narrow::System<M> system;

  size_t const geoA_idx = system.create_geometry();
  size_t const geoB_idx = system.create_geometry();

  MyObject objA;
  MyObject objB;

  objA.set_geometry_idx( geoA_idx );
  objB.set_geometry_idx( geoB_idx );

  geometry_type  & geoA = system.get_geometry(objA);
  geometry_type  & geoB = system.get_geometry(objB);
  
  geoA.add_shape( shapeA ); 
  geoB.add_shape( shapeB );


  // The part in the block below is specific for tetrahedra  meshes.
  //
  // Their BVH needs to be created in the "object" prior to using
  // them (as pre-processing before simulation) by using make_kdop_bvh. During
  // run-time while simulating their BVH must be updated to reflect positional
  // changes. This is done using the update_kdop_bvh.
  {
    std::vector< narrow::KDopBvhUpdateWorkItem< M > > work_pool;

    if( geoA.m_tetramesh.has_data())
    {
      narrow::make_kdop_bvh( system.params(), objA, geoA);
      narrow::KDopBvhUpdateWorkItem<M> const work = narrow::KDopBvhUpdateWorkItem<M>( objA, geoA, posA, rotA );

      work_pool.push_back(work);
    }
    if( geoB.m_tetramesh.has_data())
    {
      narrow::make_kdop_bvh( system.params(), objB, geoB);

      narrow::KDopBvhUpdateWorkItem<M> const work = narrow::KDopBvhUpdateWorkItem<M>( objB, geoB, posB, rotB );
      work_pool.push_back(work);
    }

    if( work_pool.size() )
    {
      narrow::update_kdop_bvh(work_pool );
    }
  }

  MyCallback callback;

  narrow::TestPair<M> test_pair = narrow::TestPair<M>( objA, objB, posA, rotA, posB, rotB, callback );

  std::vector<narrow::TestPair<M> > test_pairs;
  test_pairs.push_back(test_pair);

  narrow::dispatch_collision_handlers( system, test_pairs );

  BOOST_CHECK_EQUAL(callback.m_hit, should_overlap);

}


BOOST_AUTO_TEST_SUITE(narrow);

BOOST_AUTO_TEST_CASE(dispatch_sphere_box_test)
{	
  sphere_type shapeA;
  shapeA.radius()    = 1.0f;
  shapeA.transform() = C::identity();    
  
  box_type shapeB;
  shapeB.half_extent()(0) = 1.0f;
  shapeB.half_extent()(1) = 1.0f;
  shapeB.half_extent()(2) = 1.0f;
  
  shapeB.transform() = C::identity();    
  
  {
    V posA = V::make(-1.0f, 0.0f, 0.0f);
    V posB = V::make( 1.0f, 0.0f, 0.0f);
    Q rotA = Q::identity();
    Q rotB = Q::identity();
    do_test( posA, posB, rotA, rotB, shapeA, shapeB, true);
  }
  
  {
    V posA = V::make(-2.0f, 0.0f, 0.0f);
    V posB = V::make( 2.0f, 0.0f, 0.0f);
    Q rotA = Q::identity();
    Q rotB = Q::identity();
    do_test( posA, posB, rotA, rotB, shapeA, shapeB, false);
  }  
}

BOOST_AUTO_TEST_CASE(dispatch_sphere_sphere_test)
{	
  sphere_type shapeA;
  shapeA.radius()    = 1.0f;
  shapeA.transform() = C::identity();    
  
  sphere_type shapeB;
  shapeB.radius()    = 1.0f;
  shapeB.transform() = C::identity();    
  
  {
    V posA = V::make(-1.0f, 0.0f, 0.0f);
    V posB = V::make( 1.0f, 0.0f, 0.0f);
    Q rotA = Q::identity();
    Q rotB = Q::identity();
    do_test( posA, posB, rotA, rotB, shapeA, shapeB, true);
  }
  
  {
    V posA = V::make(-2.0f, 0.0f, 0.0f);
    V posB = V::make( 2.0f, 0.0f, 0.0f);
    Q rotA = Q::identity();
    Q rotB = Q::identity();
    do_test( posA, posB, rotA, rotB, shapeA, shapeB, false);
  }  
}

BOOST_AUTO_TEST_CASE(dispatch_box_box_test)
{	
  box_type shapeA;
  shapeA.half_extent()(0) = 1.0f;
  shapeA.half_extent()(1) = 1.0f;
  shapeA.half_extent()(2) = 1.0f;
  shapeA.transform()      = C::identity();    
  
  box_type shapeB;
  shapeB.half_extent()(0) = 1.0f;
  shapeB.half_extent()(1) = 1.0f;
  shapeB.half_extent()(2) = 1.0f;
  shapeB.transform()      = C::identity();    
  
  {
    V posA = V::make(-1.0f, 0.0f, 0.0f);
    V posB = V::make( 1.0f, 0.0f, 0.0f);
    Q rotA = Q::identity();
    Q rotB = Q::identity();
    do_test( posA, posB, rotA, rotB, shapeA, shapeB, true);
  }
  
  {
    V posA = V::make(-2.0f, 0.0f, 0.0f);
    V posB = V::make( 2.0f, 0.0f, 0.0f);
    Q rotA = Q::identity();
    Q rotB = Q::identity();
    do_test( posA, posB, rotA, rotB, shapeA, shapeB, false);
  }  
}

BOOST_AUTO_TEST_CASE(dispatch_box_sphere_test)
{	
  box_type shapeA;
  shapeA.half_extent()(0) = 1.0f;
  shapeA.half_extent()(1) = 1.0f;
  shapeA.half_extent()(2) = 1.0f;
  shapeA.transform()      = C::identity();    
  
  sphere_type shapeB;
  shapeB.radius()    = 1.0f;
  shapeB.transform() = C::identity();    
  
  {
    V posA = V::make(-1.0f, 0.0f, 0.0f);
    V posB = V::make( 1.0f, 0.0f, 0.0f);
    Q rotA = Q::identity();
    Q rotB = Q::identity();
    do_test( posA, posB, rotA, rotB, shapeA, shapeB, true);
  }
  
  {
    V posA = V::make(-2.0f, 0.0f, 0.0f);
    V posB = V::make( 2.0f, 0.0f, 0.0f);
    Q rotA = Q::identity();
    Q rotB = Q::identity();
    do_test( posA, posB, rotA, rotB, shapeA, shapeB, false);
  }  
}

BOOST_AUTO_TEST_CASE(dispatch_tetmesh_tetmesh_test)
{
  mesh_array::T3Mesh surfaceA;
  mesh_array::VertexAttribute<T,mesh_array::T3Mesh> sXA;
  mesh_array::VertexAttribute<T,mesh_array::T3Mesh> sYA;
  mesh_array::VertexAttribute<T,mesh_array::T3Mesh> sZA;

  mesh_array::make_box<M>( 2.0f, 2.0f, 2.0f, surfaceA, sXA, sYA, sZA);

  mesh_array::T4Mesh mesh;
  mesh_array::VertexAttribute<T,mesh_array::T4Mesh> X;
  mesh_array::VertexAttribute<T,mesh_array::T4Mesh> Y;
  mesh_array::VertexAttribute<T,mesh_array::T4Mesh> Z;

  mesh_array::tetgen(surfaceA, sXA, sYA, sZA, mesh, X, Y, Z);

  tetramesh_type shapeA;
  tetramesh_type shapeB;

  shapeA.set_tetramesh_shape(mesh, X, Y, Z);
  shapeB.set_tetramesh_shape(mesh, X, Y, Z);


  {
    V posA = V::make(-1.0f, 0.0f, 0.0f);
    V posB = V::make( 1.0f, 0.0f, 0.0f);
    Q rotA = Q::identity();
    Q rotB = Q::identity();
    do_test( posA, posB, rotA, rotB, shapeA, shapeB, true);
  }

  {
    V posA = V::make(-2.0f, 0.0f, 0.0f);
    V posB = V::make( 2.0f, 0.0f, 0.0f);
    Q rotA = Q::identity();
    Q rotB = Q::identity();
    do_test( posA, posB, rotA, rotB, shapeA, shapeB, false);
  }
}

BOOST_AUTO_TEST_CASE(dispatch_sphere_tetmesh_test)
{
  sphere_type shapeA;
  shapeA.radius()    = 1.0f;
  shapeA.transform() = C::identity();

  tetramesh_type shapeB;

  mesh_array::T3Mesh surface;
  mesh_array::VertexAttribute<T,mesh_array::T3Mesh> sX;
  mesh_array::VertexAttribute<T,mesh_array::T3Mesh> sY;
  mesh_array::VertexAttribute<T,mesh_array::T3Mesh> sZ;

  mesh_array::make_box<M>( 2.0f, 2.0f, 2.0f, surface, sX, sY, sZ);

  mesh_array::T4Mesh mesh;
  mesh_array::VertexAttribute<T,mesh_array::T4Mesh> X;
  mesh_array::VertexAttribute<T,mesh_array::T4Mesh> Y;
  mesh_array::VertexAttribute<T,mesh_array::T4Mesh> Z;

  mesh_array::tetgen(surface, sX, sY, sZ, mesh, X, Y, Z);

  shapeB.set_tetramesh_shape(mesh, X, Y, Z);

  {
    V posA = V::make(  1.0f, 0.0f, 0.0f);
    V posB = V::make( -0.9f, 0.0f, 0.0f);
    Q rotA = Q::identity();
    Q rotB = Q::identity();
    do_test( posA, posB, rotA, rotB, shapeA, shapeB, true);
  }

  {
    V posA = V::make( 1.1f, 0.0f, 0.0f);
    V posB = V::make(-1.1f, 0.0f, 0.0f);
    Q rotA = Q::identity();
    Q rotB = Q::identity();
    do_test( posA, posB, rotA, rotB, shapeA, shapeB, false);
  }
}


BOOST_AUTO_TEST_CASE(dispatch_convex_convex_test)
{	
  convex_type shapeA;
  shapeA.data().add_point( V::make( 1.0f,  1.0f,  1.0f ) );
  shapeA.data().add_point( V::make(-1.0f,  1.0f,  1.0f ) );
  shapeA.data().add_point( V::make( 1.0f, -1.0f,  1.0f ) );
  shapeA.data().add_point( V::make(-1.0f, -1.0f,  1.0f ) );
  shapeA.data().add_point( V::make( 1.0f,  1.0f, -1.0f ) );
  shapeA.data().add_point( V::make(-1.0f,  1.0f, -1.0f ) );
  shapeA.data().add_point( V::make( 1.0f, -1.0f, -1.0f ) );
  shapeA.data().add_point( V::make(-1.0f, -1.0f, -1.0f ) );
  shapeA.transform() = C::identity();    
  
  convex_type shapeB;
  shapeB.data().add_point( V::make( 1.0f,  1.0f,  1.0f ) );
  shapeB.data().add_point( V::make(-1.0f,  1.0f,  1.0f ) );
  shapeB.data().add_point( V::make( 1.0f, -1.0f,  1.0f ) );
  shapeB.data().add_point( V::make(-1.0f, -1.0f,  1.0f ) );
  shapeB.data().add_point( V::make( 1.0f,  1.0f, -1.0f ) );
  shapeB.data().add_point( V::make(-1.0f,  1.0f, -1.0f ) );
  shapeB.data().add_point( V::make( 1.0f, -1.0f, -1.0f ) );
  shapeB.data().add_point( V::make(-1.0f, -1.0f, -1.0f ) );
  shapeB.transform() = C::identity();    
  
  {
    V posA = V::make(-1.0f, 0.0f, 0.0f);
    V posB = V::make( 1.0f, 0.0f, 0.0f);
    Q rotA = Q::identity();
    Q rotB = Q::identity();
    // 2015-11-30 Kenny: Convex types not yet supported
    //    do_test( posA, posB, rotA, rotB, shapeA, shapeB, true);
  }
  
  {
    V posA = V::make(-2.0f, 0.0f, 0.0f);
    V posB = V::make( 2.0f, 0.0f, 0.0f);
    Q rotA = Q::identity();
    Q rotB = Q::identity();
    // 2015-11-30 Kenny: Convex types not yet supported
    // do_test( posA, posB, rotA, rotB, shapeA, shapeB, false);
  }  
}


BOOST_AUTO_TEST_SUITE_END();
