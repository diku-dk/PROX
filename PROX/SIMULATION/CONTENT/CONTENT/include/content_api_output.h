#ifndef CONTENT_API_OUTPUT_H
#define CONTENT_API_OUTPUT_H

#include <string>

namespace content
{
  
  class Output
    {
    public:
      
      Output(){}
      virtual ~Output(){}
      
    public:
      
      // Getters
      /**
       * Get the size of the collision envelope.
       *
       * @return   The value of the collision envelope.
       */
      virtual float get_collision_envelope() = 0;
      
      /**
       * Get the simulation time-step size.
       * A simulator may use an internal time step size which is much smaller
       * than for instance the time-step between two consecutive frames in
       * an animation. Say one wants 30 fps then the time-step in the animation
       * is ~0.01333 seconds. However, one may wish to use a simulation time
       * step of say 0.001 seconds.
       *
       * @return      The value of the time step.
       */
      virtual float get_time_step() = 0;
      
      /**
       * Get total simulated time.
       *
       * @return simulated time.
       */
      virtual float get_time() = 0;
      
      /**
       * Get gravity up direction parameters.
       *
       * @param x             Upon return this argument holds the x-coordinate of the up direction.
       * @param y             Upon return this argument holds the y-coordinate of the up direction.
       * @param z             Upon return this argument holds the z-coordinate of the up direction.
       */
      virtual void get_gravity_up(
                                float & x
                               , float & y
                               , float & z
                               ) = 0;

      /**
       * Get gravity acceleration.
       */
      virtual void get_gravity_acceleration( float & acceleration ) = 0;

      /**
       * Get linear and angular damping parameters.
       */
      virtual void get_damping_parameters( float & linear, float & angular ) = 0;


      /**
       * Get material name.
       *
       * @param material_index    An unique index value among the collection of materials which represents the material.
       * @return                A human readable name for the material or the empty string if no such name is known.
       */
      virtual std::string get_material_name(size_t const & material_index) = 0;
      
      /**
       * Get the number of materials.
       *
       * @return   The number of materials.
       */
      virtual size_t get_number_of_materials() = 0;
      
      /**
       * Get Material Indices.
       *
       * @param index_array     Upon return this array will hold the indices of all
       *                        materials. The array must be pre-allocated by the caller
       *                        and must contain at least 'get_number_of_materials()' entries.
       */
      virtual void get_material_indices( size_t * index_array ) = 0;
      
      /**
       * Get Number of Material Properties.
       *
       * @return    The number of defined material properties.
       */
      virtual size_t get_number_of_properties() = 0;
      
      /**
       * Get Material Property Indices.
       * A material property is defined uniquely by a pair of un-ordered material indices. 
       *
       * @param first_index_array     Upon return this array will hold the indices of all
       *                            the first materials. The array must be pre-allocated by the caller
       *                            and must contain at least 'get_number_of_properties()' entries.
       * @param second_index_array    Upon return this array will hold the indices of all
       *                            the first materials. The array must be pre-allocated by the caller
       *                            and must contain at least 'get_number_of_properties()' entries.
       */
      virtual void get_material_property_indices(  size_t * first_index_array
                                                 , size_t * second_index_array 
                                                 ) = 0;
      
      /**
       * Get Master Material.
       *
       * @param first_index    The first material index identifying the material property.
       * @param second_index   The second material index identifying the material property.
       * @return             The master material index.
       */
      virtual size_t get_master( size_t const & first_index
                                , size_t const & second_index 
                                ) = 0;
      
      /**
       * Get Friction Coefficients.
       *
       * @param first_index    The first material index identifying the material property.
       * @param second_index   The second material index identifying the material property.
       * @param x            Upon return this argument holds the value of the x-direction friction coefficient.
       * @param y            Upon return this argument holds the value of the y-direction friction coefficient.
       * @param z            Upon return this argument holds the value of the z-direction friction coefficient.
       */
      virtual void get_friction( size_t const & first_index
                                , size_t const & second_index
                                , float & x
                                , float & y
                                , float & z 
                                ) = 0;
      
      /**
       * Get Master Direction.
       *
       * @param first_index    The first material index identifying the material property.
       * @param second_index   The second material index identifying the material property.
       * @param x            Upon return this argument holds the value of the x-direction.
       * @param y            Upon return this argument holds the value of the y-direction.
       * @param z            Upon return this argument holds the value of the z-direction.
       */
      virtual void get_master_direction( size_t const & first_index
                                        , size_t const & second_index
                                        , float & x
                                        , float & y
                                        , float & z 
                                        ) = 0;
      
      /**
       * Get Restitution Coefficient.
       *
       * @param first_index    The first material index identifying the material property.
       * @param second_index   The second material index identifying the material property.
       * @return             The resitution coefficient.
       */
      virtual float get_restitution( size_t const & first_index
                                    , size_t const & second_index 
                                    ) = 0;
      
      /**
       * Get the number of geometries.
       *
       * @return    The number of geometries.
       */
      virtual size_t get_number_of_geometries() = 0;
      
      /**
       * Get geometry indices.
       *
       * @param index_array     Upon return this array holds all the indices representing
       *                        the geometries. The array must be preallocated by the caller
       *                        and contain at least 'get_number_of_geometries()' entries.
       */
      virtual void get_geometry_indices( size_t * index_array ) = 0;
      
      /**
       * Get Geometry Name.
       *
       * @param geometry_index      The unique geometry index.
       * @return                  A human readable string corresponding to the geometry index or
       *                          the empty string if the name is unknown.
       */
      virtual std::string get_geometry_name(size_t const & geometry_index) = 0;
      
      /**
       * Get number of boxes in the geometry.
       *
       * @param geometry_index      The index of the geometry.
       * @return                  The number of boxes in the geometry.
       */
      virtual size_t get_number_of_boxes( size_t const & geometry_index ) = 0;
      
      /**
       * Get Box Shape Parameters.
       *
       * @param geoemtry_index       The index of the geometry.
       * @param box_number         The number of the box.
       * @param width              Upon return this argument holds the value of the width of the box (the x-extent).
       * @param height             Upon return this argument holds the value of the height of the box (the y-extent).
       * @param depth              Upon return this argument holds the value of the depth of the box (the z-extent).
       */
      virtual void get_box_shape( size_t const & geometry_index
                                 , size_t const & box_number
                                 ,  float & width
                                 , float & height
                                 , float & depth
                                 ) = 0;
      
      /**
       * Get Box Shape Position.
       * 
       * @param geoemtry_index       The index of the geometry.
       * @param box_number         The number of the box.
       * @param x                  Upon return this argument holds the x-coordinate.
       * @param y                  Upon return this argument holds the y-coordinate.
       * @param z                  Upon return this argument holds the z-coordinate.
       */
      virtual void get_box_position( size_t const & geometry_index
                                    , size_t const & box_number
                                    , float & x
                                    , float & y
                                    , float & z
                                    ) = 0;
      
      /**
       * Get Box Shape Orientation.
       * 
       * @param geoemtry_index      The index of the geometry.
       * @param box_number          The number of the box.
       * @param Qs                  Upon return this argument holds the real part of the rotation quaternion.
       * @param Qx                  Upon return this argument holds the first imaginary part of the rotation quaternion.
       * @param Qy                  Upon return this argument holds the second imaginary part of the rotation quaternion.
       * @param Qz                  Upon return this argument holds the third imaginary part of the rotation quaternion.
       */
      virtual void get_box_orientation( size_t const & geometry_index
                                       , size_t const & box_number
                                       , float & Qs
                                       , float & Qx
                                       , float & Qy
                                       , float & Qz
                                       ) = 0;
      
      /**
       * Get number of cones in the geometry.
       *
       * @param geometry_index      The index of the geometry.
       * @return                  The number of cones in the geometry.
       */
      virtual size_t get_number_of_cones( size_t const & geometry_index ) = 0;
      
      /**
       * Get Cone Shape Parameters.
       *
       * @param geoemtry_index       The index of the geometry.
       * @param cone_number        The number of the cone.
       * @param radius             Upon return this argument holds the value of the base radius of the cone.
       * @param height             Upon return this argument holds the value of the height of the cone.
       */
      virtual void get_cone_shape( size_t const & geometry_index
                                  , size_t const & cone_number
                                  , float & radius
                                  , float & height
                                  ) = 0;
      
      /**
       * Get Cone Shape Position.
       * 
       * @param geoemtry_index       The index of the geometry.
       * @param cone_number        The number of the cone.
       * @param x                  Upon return this argument holds the x-coordinate.
       * @param y                  Upon return this argument holds the y-coordinate.
       * @param z                  Upon return this argument holds the z-coordinate.
       */
      virtual void get_cone_position(  size_t const & geometry_index
                                     , size_t const & cone_number
                                     , float & x
                                     , float & y
                                     , float & z
                                     ) = 0;
      
      /**
       * Get Cone Shape Orientation.
       * 
       * @param geoemtry_index    The index of the geometry.
       * @param cone_number     The number of the cone.
       * @param Qs                  Upon return this argument holds the real part of the rotation quaternion.
       * @param Qx                  Upon return this argument holds the first imaginary part of the rotation quaternion.
       * @param Qy                  Upon return this argument holds the second imaginary part of the rotation quaternion.
       * @param Qz                  Upon return this argument holds the third imaginary part of the rotation quaternion.
       */
      virtual void get_cone_orientation(  size_t const & geometry_index
                                        , size_t const & cone_number
                                        , float & Qs
                                        , float & Qx
                                        , float & Qy
                                        , float & Qz
                                        ) = 0;     
      
      /**
       * Get number of capsules in the geometry.
       *
       * @param geometry_index      The index of the geometry.
       * @return                  The number of capsules in the geometry.
       */
      virtual size_t get_number_of_capsules( size_t const & geometry_index ) = 0;
      
      /**
       * Get Capsule Shape Parameters.
       *
       * @param geoemtry_index       The index of the geometry.
       * @param capsule_number     The number of the capsule.
       * @param radius             Upon return this argument holds the value of the radius of the capsule.
       * @param height             Upon return this argument holds the value of the height of the capsule.
       */
      virtual void get_capsule_shape(  size_t const & geometry_index
                                     , size_t const & capsule_number
                                     , float & radius
                                     , float & height
                                     ) = 0;
      
      /**
       * Get Capsule Shape Position.
       * 
       * @param geoemtry_index       The index of the geometry.
       * @param capsule_number     The number of the capsule.
       * @param x                  Upon return this argument holds the x-coordinate.
       * @param y                  Upon return this argument holds the y-coordinate.
       * @param z                  Upon return this argument holds the z-coordinate.
       */
      virtual void get_capsule_position(  size_t const & geometry_index
                                        , size_t const & capsule_number
                                        , float & x
                                        , float & y
                                        , float & z
                                        ) = 0;
      
      /**
       * Get Capsule Shape Orientation.
       * 
       * @param geoemtry_index    The index of the geometry.
       * @param capsule_number  The number of the capsule.
       * @param Qs                  Upon return this argument holds the real part of the rotation quaternion.
       * @param Qx                  Upon return this argument holds the first imaginary part of the rotation quaternion.
       * @param Qy                  Upon return this argument holds the second imaginary part of the rotation quaternion.
       * @param Qz                  Upon return this argument holds the third imaginary part of the rotation quaternion.
       */
      virtual void get_capsule_orientation(  size_t const & geometry_index
                                           , size_t const & capsule_number
                                           , float & Qs
                                           , float & Qx
                                           , float & Qy
                                           , float & Qz
                                           ) = 0;       
      
      
      /**
       * Get number of cylinders in the geometry.
       *
       * @param geometry_index      The index of the geometry.
       * @return                  The number of cylinders in the geometry.
       */
      virtual size_t get_number_of_cylinders( size_t const & geometry_index ) = 0;
      
      /**
       * Get Cylinder Shape Parameters.
       *
       * @param geoemtry_index       The index of the geometry.
       * @param cylinder_number     The number of the cylinder.
       * @param radius             Upon return this argument holds the value of the radius of the cylinder.
       * @param height             Upon return this argument holds the value of the height of the cylinder.
       */
      virtual void get_cylinder_shape( size_t const & geometry_index
                                      , size_t const & cylinder_number
                                      , float & radius
                                      , float & height
                                      ) = 0;
      
      /**
       * Get Cylinder Shape Position.
       * 
       * @param geoemtry_index       The index of the geometry.
       * @param cylinder_number    The number of the cylinder.
       * @param x                  Upon return this argument holds the x-coordinate.
       * @param y                  Upon return this argument holds the y-coordinate.
       * @param z                  Upon return this argument holds the z-coordinate.
       */
      virtual void get_cylinder_position(size_t const & geometry_index
                                         , size_t const & cylinder_number
                                         , float & x
                                         , float & y
                                         , float & z
                                         ) = 0;
      
      /**
       * Get Cylinder Shape Orientation.
       * 
       * @param geoemtry_index    The index of the geometry.
       * @param cylinder_number The number of the cylinder.
       * @param Qs                  Upon return this argument holds the real part of the rotation quaternion.
       * @param Qx                  Upon return this argument holds the first imaginary part of the rotation quaternion.
       * @param Qy                  Upon return this argument holds the second imaginary part of the rotation quaternion.
       * @param Qz                  Upon return this argument holds the third imaginary part of the rotation quaternion.
       */
      virtual void get_cylinder_orientation(  size_t const & geometry_index
                                            , size_t const & cylinder_number
                                            , float & Qs
                                            , float & Qx
                                            , float & Qy
                                            , float & Qz
                                            ) = 0;   
      /**
       * Get number of ellipsoids in the geometry.
       *
       * @param geometry_index      The index of the geometry.
       * @return                  The number of ellipsoids in the geometry.
       */
      virtual size_t get_number_of_ellipsoids( size_t const & geometry_index ) = 0;
      
      /**
       * Get Ellipsoid Shape Parameters.
       *
       * @param geoemtry_index       The index of the geometry.
       * @param ellipsoid_number   The number of the ellipsoid.
       * @param scale_x            Upon return this argument holds the value of the scaling along the x-axis.
       * @param scale_y            Upon return this argument holds the value of the scaling along the y-axis.
       * @param scale_z            Upon return this argument holds the value of the scaling along the z-axis.
       */
      virtual void get_ellipsoid_shape( size_t const & geometry_index
                                       , size_t const & ellipsoid_number
                                       , float & scale_x
                                       , float & scale_y
                                       , float & scale_z
                                       ) = 0;
      
      /**
       * Get Ellipsoid Shape Position.
       * 
       * @param geoemtry_index       The index of the geometry.
       * @param ellipsoid_number   The number of the ellipsoid.
       * @param x                  Upon return this argument holds the x-coordinate.
       * @param y                  Upon return this argument holds the y-coordinate.
       * @param z                  Upon return this argument holds the z-coordinate.
       */
      virtual void get_ellipsoid_position(size_t const & geometry_index
                                          , size_t const & ellipsoid_number
                                          , float & x
                                          , float & y
                                          , float & z
                                          ) = 0;
      
      /**
       * Get Ellipsoid Shape Orientation.
       * 
       * @param geoemtry_index     The index of the geometry.
       * @param ellipsoid_number The number of the ellipsoid.
       * @param Qs                  Upon return this argument holds the real part of the rotation quaternion.
       * @param Qx                  Upon return this argument holds the first imaginary part of the rotation quaternion.
       * @param Qy                  Upon return this argument holds the second imaginary part of the rotation quaternion.
       * @param Qz                  Upon return this argument holds the third imaginary part of the rotation quaternion.
       */
      virtual void get_ellipsoid_orientation(  size_t const & geometry_index
                                             , size_t const & ellipsoid_number
                                             , float & Qs
                                             , float & Qx
                                             , float & Qy
                                             , float & Qz
                                             ) = 0; 
      
      /**
       * Get number of spheres in the geometry.
       *
       * @param geometry_index      The index of the geometry.
       * @return                  The number of spheres in the geometry.
       */
      virtual size_t get_number_of_spheres( size_t const & geometry_index ) = 0;
      
      /**
       * Get Sphere Shape Parameters.
       *
       * @param geoemtry_index       The index of the geometry.
       * @param sphere_number      The number of the sphere.
       * @param radius             Upon return this argument holds the value of the radius of the sphere.
       */
      virtual void get_sphere_shape( size_t const & geometry_index
                                    , size_t const & sphere_number
                                    , float & radius
                                    ) = 0;
      
      /**
       * Get Sphere Shape Position.
       * 
       * @param geoemtry_index       The index of the geometry.
       * @param sphere_number      The number of the sphere.
       * @param x                  Upon return this argument holds the x-coordinate.
       * @param y                  Upon return this argument holds the y-coordinate.
       * @param z                  Upon return this argument holds the z-coordinate.
       */
      virtual void get_sphere_position(size_t const & geometry_index
                                       , size_t const & sphere_number
                                       , float & x
                                       , float & y
                                       , float & z
                                       ) = 0;
      
      /**
       * Get Sphere Shape Orientation.
       * 
       * @param geoemtry_index    The index of the geometry.
       * @param sphere_number   The number of the sphere.
       * @param Qs                  Upon return this argument holds the real part of the rotation quaternion.
       * @param Qx                  Upon return this argument holds the first imaginary part of the rotation quaternion.
       * @param Qy                  Upon return this argument holds the second imaginary part of the rotation quaternion.
       * @param Qz                  Upon return this argument holds the third imaginary part of the rotation quaternion.
       */
      virtual void get_sphere_orientation(  size_t const & geometry_index
                                          , size_t const & sphere_number
                                          , float & Qs
                                          , float & Qx
                                          , float & Qy
                                          , float & Qz
                                          ) = 0; 
      
      
      /**
       * Get number of convexes in the geometry.
       *
       * @param geometry_index      The index of the geometry.
       * @return                  The number of convexes in the geometry.
       */
      virtual size_t get_number_of_convexes( size_t const & geometry_index ) = 0;
      
      /**
       * Get Convex Shape Parameters.
       *
       * @param geoemtry_index       The index of the geometry.
       * @param convex_number      The number of the convex.
       * @param no_points          The number of points in the convex hull.
       */      
      virtual void get_convex_shape( size_t const & geometry_index
                                    , size_t const & convex_number
                                    , size_t & no_points
                                    ) = 0;
      
      /**
       * Get Convex Shape Parameters.
       *
       * @param geoemtry_index       The index of the geometry.
       * @param convex_number      The number of the convex.
       * @param coordinates        Upon return this argument holds the coordinates of all the points. The array must be preallocated by the caller and contain atleast 3 times the 'get_convex_shape(...,no_points)' entries.
       */      
      virtual void get_convex_shape( size_t const & geometry_index
                                    , size_t const & convex_number
                                    , float * coordinates
                                    ) = 0;
      
      /**
       * Get Convex Shape Position.
       * 
       * @param geoemtry_index       The index of the geometry.
       * @param convex_number      The number of the convex.
       * @param x                  Upon return this argument holds the x-coordinate.
       * @param y                  Upon return this argument holds the y-coordinate.
       * @param z                  Upon return this argument holds the z-coordinate.
       */
      virtual void get_convex_position(  size_t const & geometry_index
                                       , size_t const & convex_number
                                       , float & x
                                       , float & y
                                       , float & z
                                       ) = 0;
      
      /**
       * Get Convex Shape Orientation.
       * 
       * @param geoemtry_index    The index of the geometry.
       * @param convex_number   The number of the convex.
       * @param Qs                  Upon return this argument holds the real part of the rotation quaternion.
       * @param Qx                  Upon return this argument holds the first imaginary part of the rotation quaternion.
       * @param Qy                  Upon return this argument holds the second imaginary part of the rotation quaternion.
       * @param Qz                  Upon return this argument holds the third imaginary part of the rotation quaternion.
       */
      virtual void get_convex_orientation(  size_t const & geometry_index
                                          , size_t const & convex_number
                                          , float & Qs
                                          , float & Qx
                                          , float & Qy
                                          , float & Qz
                                          ) = 0;         


      /**
       * Get number of tetrameshes in the geometry.
       *
       * @param geometry_index    The index of the geometry.
       * @return                  The number of tetrameshes in the geometry. This can only be 1 or 0 as a
       *                          geometry only can have one tetramesh shape.
       */
      virtual size_t get_number_of_tetrameshes( size_t const & geometry_index ) = 0;

      /**
       * Get Tetramesh Size Parameters.
       * Observe that a shape can only have one tetramesh and tetrameshes live in body space by definition
       * (hence no body space transformations).
       *
       * @param geometry_idx    The index representing the geometry.
       * @param N               The number of vertices in the tetrahedra mesh.
       * @param K               The number of tetrahedra in the tetrahedra mesh.
       */
      virtual void get_tetramesh_shape(  size_t const & geometry_idx
                                       , size_t       & N
                                       , size_t      & K
                                       ) = 0;

      /**
       * Get Tetramesh Shape Parameters.
       * Observe that a shape can only have one tetramesh and tetrameshes live in body space by definition
       * (hence no body space transformations).
       *
       * @param geometry_idx    The index representing the geometry.
       * @param vertices        An array of vertex indices must be allocated to be N long.
       * @param tetrahedra      An array of tetrahedra indices must be allocated to be 4*K long.
       *                        Each consequtive four tuple defines the four vertex indices of a
       *                        single tetrahedron
       * @param coordinates     Array of coordinates in order x1 y1 z1... xN yN zN. Must be allocated
       *                        to be 3 N long.
       */
      virtual void get_tetramesh_shape(  size_t const & geometry_idx
                                       , size_t * vertices
                                       , size_t * tetrahedra
                                       , float  * coordinates
                                       ) = 0;
      
      /**
       * Get Number of Rigid Bodies.
       *
       * @return The number of bodies.
       */
      virtual size_t get_number_of_rigid_bodies() = 0;
      
      /**
       * Get Object Indices of Rigid Bodies.
       *
       * @param index_array    Upon return this array will contain the object indices of
       *                       all the rigid bodies. The array must be allocated by the
       *                       caller and should contain at least 'get_number_of_rigid_bodies()'
       *                       entries.
       */
      virtual void get_rigid_body_indices( size_t * index_array ) = 0;
      
      /**
       * Get Rigid Body Name.
       *
       * @param body_index    The rigid body object index.
       * @return            A human readable name for the rigid body or the empty string if the name is unknown.
       */
      virtual std::string get_rigid_body_name( size_t const & body_index ) = 0;
      
      /**
       * Get Rigid Body Position.
       * @param body_index    The rigid body object index.
       * @param x           Upon return this argument holds the x-coordinate.
       * @param y           Upon return this argument holds the y-coordinate.
       * @param z           Upon return this argument holds the z-coordinate.
       */
      virtual void get_rigid_body_position( size_t const & body_index
                                           , float & x
                                           , float & y
                                           , float & z 
                                           ) = 0; 
      
      
      /**
       * Get Rigid Body Orientation.
       * @param body_index  The rigid body object index.
       * @param qs          Upon return this argument holds the real-part of a quaternion representation of the orientation.
       * @param qx          Upon return this argument holds the x-component of the imaginary part of a quaternion representation of the orientation.
       * @param qy          Upon return this argument holds the y-component of the imaginary part of a quaternion representation of the orientation.
       * @param qz          Upon return this argument holds the z-component of the imaginary part of a quaternion representation of the orientation.
       */
      virtual void get_rigid_body_orientation( size_t const & body_index
                                              , float & qs
                                              , float & qx
                                              , float & qy 
                                              , float & qz 
                                              ) = 0; 
      
      /**
       * Get Rigid Body Velocity.
       * @param body_index    The rigid body object index.
       * @param x           Upon return this argument holds the x-coordinate of the velocity.
       * @param y           Upon return this argument holds the y-coordinate of the velocity.
       * @param z           Upon return this argument holds the z-coordinate of the velocity.
       */
      virtual void get_rigid_body_velocity( size_t const & body_index
                                           , float & x
                                           , float & y
                                           , float & z 
                                           ) = 0; 
      
      /**
       * Get Rigid Body Spin.
       * @param body_index    The rigid body object index.
       * @param x           Upon return this argument holds the x-coordinate of the angular velocity.
       * @param y           Upon return this argument holds the y-coordinate of the angular velocity.
       * @param z           Upon return this argument holds the z-coordinate of the angular velocity.
       */
      virtual void get_rigid_body_spin( size_t const & body_index
                                       , float & x
                                       , float & y
                                       , float & z 
                                       ) = 0; 
      
      /**
       * Get Rigid Body Inertia Tensor (in the body frame).
       * @param body_index    The rigid body object index.
       * @param xx          Upon return this argument holds the xx moment of the inertia tensor.
       * @param yy          Upon return this argument holds the yy moment of the inertia tensor.
       * @param zz          Upon return this argument holds the zz moment of the inertia tensor.
       */
      virtual void get_rigid_body_inertia( size_t const & body_index
                                          , float & xx
                                          , float & yy
                                          , float & zz 
                                          ) = 0; 
      
      /**
       * Get Mass of Rigid Body.
       * @param     The rigid body object index.
       * @return    The mass of the rigid body.
       */
      virtual float get_rigid_body_mass( size_t const & body_index ) = 0;
      
      /**
       * Get fixed state of Rigid Body.
       * @param     The rigid body object index.
       * @return    Boolean flag indicating whether the rigid body is fixed.
       */
      virtual bool get_rigid_body_fixed( size_t const & body_index ) = 0;
      
      /**
       * Get active state of Rigid Body.
       * @param     The rigid body object index.
       * @return    Boolean flag indicating whether the rigid body is active.
       */
      virtual bool get_rigid_body_active( size_t const & body_index ) = 0;
      
      
      /**
       * Get Material of Rigid Body.
       * @param body_index     The rigid body object index.
       * @return             The material index of the rigid body.
       */
      virtual size_t get_rigid_body_material( size_t const & body_index ) = 0;
      
      /**
       * Get Collision Geometry of Rigid Body.
       * @param body_index    The rigid body object index.
       * @return            The index value of the collsion geometry of the rigid body.
       */
      virtual size_t get_rigid_body_collision_geometry(size_t const & body_index) = 0;
      
      /**
       * Get Number of connected external forces.
       * @param body_index    The rigid body object index.
       * @return            The number of external forces.
       */
      virtual size_t get_number_of_connected_forces(size_t const & body_index) = 0;
      
      /**
       * Get indices of connected forces..
       * @param body_index     The rigid body object index.
       * @param index_array  Upon return this array holds the force indices of all
       *                     external forces connected to this rigid body. The array
       *                     must be allocated by the caller and should contain at
       *                     least 'get_number_of_connected_forces()' entries.
       */
      virtual void get_connected_force_indices(  size_t const & body_index
                                               , size_t * index_array 
                                               ) = 0;

      /**
       * Get current axis aligned bounding box of specified rigid body.
       */
      virtual void get_rigid_body_bounding_box(
                                       size_t const & body_index
                                       , float & min_x
                                       , float & min_y
                                       , float & min_z
                                       , float & max_x
                                       , float & max_y
                                       , float & max_z
                                       ) = 0;





      /**
       * Get the number of pin forces.
       *
       * @return   The number of pin forces.
       */
      virtual size_t get_number_of_pin_forces() const = 0;

      /**
       * Get indices of pin forces..
       *
       * @param index_array  Upon return this array holds the pin force indices. The array
       *                     must be allocated by the caller and should contain at
       *                     least 'get_number_of_pin_forces()' entries.
       */
      virtual void get_pin_force_indices( size_t * index_array ) = 0;

      /**
       * Get characteristic time for pin force.
       *
       * @param tau   The value of the characteristic time.
       */
      virtual void get_pin_tau(size_t const & force_idx, float & tau ) = 0;

      /**
       * Get target position for pin force.
       *
       * @param x  The x-coordinate.
       * @param y  The y-coordinate.
       * @param z  The z-coordinate.
       */
      virtual void get_pin_target(size_t const & force_idx, float  & x, float  & y, float  & z ) = 0;

      /**
       * Get anchor position for pin force.
       *
       * @param x  The x-coordinate.
       * @param y  The y-coordinate.
       * @param z  The z-coordinate.
       */
      virtual void get_pin_anchor(size_t const & force_idx, float  & x, float  & y, float  & z ) = 0;

      /**
       * Test if rigid body is a scripted body.
       * A scripted body is not a free moving rigid body, instead its motion is completely determined by some given close form solution.
       *
       * @param body_index    The index of body
       *
       * @return              If the rigid body is scripted then the return value is true otherwise it is fale.
       *
       */
      virtual bool has_scripted_motion( size_t const & body_index ) = 0;

      /**
       * Get scripted motion index of scripted body.
       *
       * @param body_index    The index of the body whos connected scripted motion is wanted.
       *
       * @return              A unique index represented the scripted motion
       */
      virtual size_t get_scripted_motion( size_t const & body_index ) = 0;

      /**
       * Get total number of scripted motions.
       *
       * @return The number of scripted motions.
       */
      virtual size_t get_number_of_scripted_motions() = 0;

      /**
       * Get scripted motion indices.
       *
       * @param index_array  Upon return this array holds the scripted motion indices. The array
       *                     must be allocated by the caller and should contain at
       *                     least 'get_number_of_scripted_motion()' entries.
       */
      virtual void get_scripted_motion_indices( size_t * index_array ) = 0;

      /**
       * Test the type of scripted motion.
       *
       * @param motion_index    The unique identifier of the scripted motion.
       *
       * @return                True if the motion was a key frame scripted motion.
       */
      virtual bool is_scripted_motion_keyframe( size_t const & motion_index ) = 0;

      /**
       * Test the type of scripted motion.
       *
       * @param motion_index    The unique identifier of the scripted motion.
       *
       * @return                True if the motion was a oscilation scripted motion.
       */
      virtual bool is_scripted_motion_oscilation( size_t const & motion_index ) = 0;

      /**
       * Get motion parameters of harmonic oscilator scripted motion.
       *
       * @param motion_index    The unique motion index corresponding to the motion (must be of type oscilation and not say key frame)
       * @param amplitude       The amplitude of the oscilation
       * @param frequency       The frequency of the oscilation
       * @param phase           The phase of the oscilation
       * @param dir_x           The x-component of the direction vector giving the direction of oscilation.
       * @param dir_y           The y-component of the direction vector giving the direction of oscilation.
       * @param dir_z           The z-component of the direction vector giving the direction of oscilation.
       *
       */
      virtual void get_scripted_oscilation_paramters(
                                                     size_t const & motion_index
                                                     , float & amplitude
                                                     , float & frequency
                                                     , float & phase
                                                     , float & dir_x
                                                     , float & dir_y
                                                     , float & dir_z
                                                     ) = 0;


      /**
       * Get the number of key frame positions of the scripted motion.
       *
       * @param motion_index     A unique index of a keyframe scripted motion.
       *
       * @return                 The number of key ticks on positions.
       */
      virtual size_t get_number_of_key_frame_positions( size_t const & motion_index) = 0;

      /**
       * Get the number of key frame orientations of the scripted motion.
       *
       * @param motion_index     A unique index of a keyframe scripted motion.
       *
       * @return                 The number of key ticks on orientations.
       */
      virtual size_t get_number_of_key_frame_orientations( size_t const & motion_index) = 0;

      /**
       * Get all key frame positions values.
       * All arrays must be allocated by the caller and should contain at
       * least 'get_number_of_key_frame_positions()' entries.
       *
       * @param motion_index     A unique index of a keyframe scripted motion.
       * @param time_array       Upon return this array holds time values of all key ticks.
       * @param x_array          Upon return this array holds x values of all key ticks.
       * @param y_array          Upon return this array holds y values of all key ticks.
       * @param z_array          Upon return this array holds z values of all key ticks.
       */
      virtual void get_key_frame_positions(
                                           size_t const & motion_index
                                           , float * time_array
                                           , float * x_array
                                           , float * y_array
                                           , float * z_array
                                           ) = 0;

      /**
       * Get all key frame orientation values.
       * All arrays must be allocated by the caller and should contain at
       * least 'get_number_of_key_frame_positions()' entries.
       *
       * @param motion_index     A unique index of a keyframe scripted motion.
       * @param time_array       Upon return this array holds time values of all key ticks.
       * @param qs_array         Upon return this array holds qs values of all key ticks.
       * @param qx_array         Upon return this array holds qx values of all key ticks.
       * @param qy_array         Upon return this array holds qy values of all key ticks.
       * @param qz_array         Upon return this array holds qz values of all key ticks.
       */
      virtual void get_key_frame_orientations(
                                              size_t const & motion_index
                                              , float * time_array
                                              , float * qs_array
                                              , float * qx_array
                                              , float * qy_array
                                              , float * qz_array
                                              ) = 0;

    };
  
}// namespace content

// CONTENT_API_OUTPUT_H
#endif
