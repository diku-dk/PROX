#include <procedural_factory.h>

#include <tiny_math_types.h>

#include <cmath>

namespace procedural
{
	
  template<typename MT>
	void make_arch(
                 content::API *  engine
                 , typename MT::vector3_type const & position
                 , typename MT::quaternion_type const & orientation
                 , typename MT::real_type const & r_outer
                 , typename MT::real_type const & r_inner
                 , typename MT::real_type const & pillar_height
                 , typename MT::real_type const & stone_depth
                 , size_t const & arch_slices
                 , size_t const & pillar_segments
                 , MaterialInfo<typename MT::real_type> mat_info
                 )
	{
    using std::cos;
    using std::sin;
    
    typedef typename MT::real_type       T;
    typedef typename MT::vector3_type    V;
    typedef typename MT::quaternion_type  Q;
    typedef typename MT::value_traits    VT;
    
		T      const delta_theta         = VT::pi()/arch_slices;
		T      const pillar_stone_width  = pillar_height / pillar_segments;
		T      const pillar_stone_height = r_outer - r_inner;
		T      const center_radius       = ( r_outer + r_inner )*VT::half();
		T      const stone_density       = get_material_density<MT>(mat_info, "Stone");
		size_t const mid                 = get_material_id<MT>(mat_info, "Stone");
		
		std::vector<V>  arch_vertices(8u);
    
		compute_arch_stone_vertices<MT>(
                                      delta_theta
                                    , stone_depth
                                    , r_outer
                                    , r_inner
                                    , &arch_vertices[0]
                                    );
		
		GeometryHandle<MT> arch_stone   = create_geometry_handle_cuboid<MT>( engine, &arch_vertices[0] );
		GeometryHandle<MT> pillar_stone = create_geometry_handle_box<MT>( engine, pillar_stone_width, pillar_stone_height, stone_depth );
		
		for( size_t i = 0; i<pillar_segments; ++i )
		{
			T const y = ( VT::half() + i )*pillar_stone_width;

      V const T_b2m       = pillar_stone.Tb2m();
      Q const Q_b2m       = pillar_stone.Qb2m();

			V const T_left_m2l  = V::make( -center_radius, y, VT::zero() );
			Q const Q_left_m2l  = Q::Ru( VT::pi_half(),  V::k() );

			V const T_l2w       = position;
			Q const Q_l2w       = orientation;

      V T_left_b2w;
      Q Q_left_b2w;

      compute_body_to_world_transform<MT>(
                                          T_b2m
                                          , Q_b2m
                                          , T_left_m2l
                                          , Q_left_m2l
                                          , T_l2w
                                          , Q_l2w
                                          , T_left_b2w
                                          , Q_left_b2w
                                          );

			create_rigid_body<MT>(
                            engine
                            , T_left_b2w
                            , Q_left_b2w
                            , pillar_stone
                            , mid
                            , stone_density
                            );

      V const T_right_m2l = V::make(  center_radius, y, VT::zero() );
      Q const Q_right_m2l = Q::Ru( VT::pi_half(),  V::k() );

      V T_right_b2w;
      Q Q_right_b2w;

      compute_body_to_world_transform<MT>(
                                          T_b2m
                                          , Q_b2m
                                          , T_right_m2l
                                          , Q_right_m2l
                                          , T_l2w
                                          , Q_l2w
                                          , T_right_b2w
                                          , Q_right_b2w
                                          );

			create_rigid_body<MT>(
                            engine
                            , T_right_b2w
                            , Q_right_b2w
                            , pillar_stone
                            , mid
                            , stone_density
                            );
			
		}


		for(size_t i = 0;i<arch_slices;++i)
		{
			T const theta = delta_theta*( i+VT::half() );
			T const x     = center_radius * cos( theta );
			T const y     = center_radius * sin( theta ) + pillar_height;
			
      V const T_b2m = arch_stone.Tb2m();
      Q const Q_b2m = arch_stone.Qb2m();

      V const T_m2l = V::make( x, y, VT::zero() );
			Q const Q_m2l = Q::Ru(  theta - VT::pi_half() , V::k() );

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

			create_rigid_body<MT>(  engine
                        , T_b2w
                        , Q_b2w
                        , arch_stone
                        , mid
                        , stone_density
                        );
		}
	}
  
  typedef tiny::MathTypes<float>  MTf;
  
  template
	void make_arch<MTf>(
                        content::API *  engine
                        , MTf::vector3_type const & position
                        , MTf::quaternion_type const & orientation
                        , MTf::real_type const & r_outer
                        , MTf::real_type const & r_inner
                        , MTf::real_type const & pillar_height
                        , MTf::real_type const & stone_depth
                        , size_t const & arch_slices
                        , size_t const & pillar_segments
                        , MaterialInfo<MTf::real_type> mat_info
                        );
  
  
} //namespace procedural
