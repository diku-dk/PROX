#ifndef PROCEDURAL_H
#define PROCEDURAL_H

#include <procedural_noise.h>
#include <procedural_types.h>
#include <procedural_factory.h>
#include <procedural_make_obj_packing.h>
#include <procedural_make_grain_packing.h>
#include <procedural_make_sliding_point.h>
#include <procedural_make_two_points.h>
#include <procedural_make_point_in_crack.h>
#include <procedural_make_cliff_edge.h>
#include <procedural_make_internal_edge.h>

#include <tiny_math_types.h>

#include <content.h>
#include <mesh_array.h>

#include <vector>

namespace procedural
{

template <typename T>
void make_arch(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
               const T& r_outer, const T& r_inner, const T& pillar_height, const T& stone_depth,
               size_t const& arch_slices, size_t const& pillar_segments, MaterialInfo<T> mat_info);

template <typename T>
void make_cannonball(content::API* engine, const T& radius, const EigenVector3<T>& position,
                     const EigenQuaternion<T>& orientation, const EigenVector3<T>& direction, MaterialInfo<T> mat_info);

template <typename T> GeometryHandleEigen<T> make_cannonball_geometry(content::API* physics, const T& radius);

template <typename T>
size_t make_cannonball_rigid_body(content::API* physics, GeometryHandleEigen<T> const& ball,
                                  const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
                                  const EigenVector3<T>& direction, MaterialInfo<T> mat_info);

template <typename T>
void make_colosseum(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
                    const T& r_outer, const T& r_inner, size_t const& slices, size_t const& segments,
                    MaterialInfo<T> mat_info);

template <typename T>
void make_slide(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
                const T& degree, MaterialInfo<T> mat_info);

template <typename T>
void make_temple(content::API* engine, const EigenVector3<T>& positionNew, const EigenQuaternion<T>& orientationNew,
                 const T& temple_height, const T& pillar_width, size_t const& num_pillars_x,
                 size_t const& num_pillars_z, MaterialInfo<T> mat_info);

template <typename T>
void make_dome(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
               const T& r_outer, const T& r_inner, size_t const& slices, size_t const& segments,
               MaterialInfo<T> mat_info);

template <typename T>
void make_ground(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
                 MaterialInfo<T> mat_info, const T& width, const T& height, const T& depth);

template <typename T>
void make_pantheon(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
                   const T& r_outer, const T& r_inner, const T& height, size_t const& slices, size_t const& segments,
                   MaterialInfo<T> mat_info);

template <typename T>
void make_tower(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
                const T& r_outer, const T& r_inner, const T& height, size_t const& slices, size_t const& segments,
                MaterialInfo<T> mat_info, bool const& use_cubes);

template <typename T>
void make_wall(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
               const T& width, const T& height, const T& depth, size_t const& layers, size_t const& span,
               MaterialInfo<T> mat_info);

template <typename T>
void make_greek_pillar(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
                       const T& pillar_width, const T& pillar_height, const T& pillar_depth,
                       size_t const& pillar_segments, size_t const& pillar_slices, MaterialInfo<T> mat_info,
                       mesh_array::TetGenSettings tetset = mesh_array::tetgen_default_settings());

template <typename T>
void make_sphere_cube(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
                      const T& sphere_radius, size_t const& spheres, MaterialInfo<T> mat_info);

template <typename T>
void make_tetrahedral_stack(content::API* engine, const EigenVector3<T>& position,
                            const EigenQuaternion<T>& orientation, const T& cube_width, const T& cube_height,
                            const T& cube_depth, size_t const& sub_divisions, MaterialInfo<T> mat_info);
template <typename T>
void make_stack(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
                const T& stone_dim, size_t const& layers, MaterialInfo<T> mat_info);

template <typename T>
void make_heavy_sphere_light_sphere(content::API* engine, const EigenVector3<T>& position,
                                    const EigenQuaternion<T>& orientation, const T& sphere_radius, size_t const& scale,
                                    MaterialInfo<T> mat_info);

template <typename T>
void make_dropping_spheres(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
                           const T& sphere_radius, size_t const& w, size_t const& h, size_t const& d,
                           MaterialInfo<T> mat_info);

template <typename T>
void make_spheres(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
                  const T& radius, MaterialInfo<T> mat_info);

template <typename T>
void make_tetrahedron(content::API* engine, const EigenVector3<T>& one, const EigenVector3<T>& two,
                      const EigenVector3<T>& three, const EigenVector3<T>& four, const EigenVector3<T>& position,
                      const EigenQuaternion<T>& orientation, MaterialInfo<T> mat_info, bool const fixed = false);

template <typename T>
void make_tetrahedron(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
                      MaterialInfo<T> mat_info, bool const fixed = false);

template <typename T>
void make_obj(content::API* engine, std::string const& name, const T& scale, const EigenVector3<T>& position,
              const EigenQuaternion<T>& orientation, MaterialInfo<T> mat_info, bool const fixed = false,
              bool const blind_copy = false, std::string const& material = "Stone",
              mesh_array::TetGenSettings = mesh_array::tetgen_default_settings());

template <typename T>
void make_box(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
              const T& width, const T& height, const T& depth, MaterialInfo<T> mat_info, bool const fixed = false);

template <typename T>
void make_twisted_stack(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
                        const T& stone_dim, size_t const& layers, const T& degree, MaterialInfo<T> mat_info);
template <typename T>
void make_sphere_layer(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
                       const T& sphere_radius, size_t const& spheres_width, size_t const& spheres_length,
                       MaterialInfo<T> mat_info);

template <typename T>
void make_box_container(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
                        const T& width, const T& height, const T& depth, const T& wall_thickness,
                        MaterialInfo<T> mat_info);

template <typename T>
void make_sphere_packing(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
                         const T& min_radius, const T& max_radius, const T& width, const T& height, const T& depth,
                         size_t const& number_of_spheres, MaterialInfo<T> mat_info);

template <typename T>
void make_sliding_point(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
                        const T& scene_size, MaterialInfo<T> const& mat_info, bool const& use_spike,
                        bool const& use_wedge, mesh_array::TetGenSettings tetset);

template <typename T>
void make_two_points(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
                     const T& scene_size, MaterialInfo<T> const& mat_info, bool const& use_spike, bool const& use_wedge,
                     bool const& use_spike_and_wedge, mesh_array::TetGenSettings tetset);

template <typename MT>
void make_point_in_crack(content::API* engine, typename MT::vector3_type const& position,
                         typename MT::quaternion_type const& orientation, typename MT::real_type const& scene_size,
                         MaterialInfo<typename MT::real_type> const& mat_info, mesh_array::TetGenSettings tetset);

template <typename MT>
void make_cliff_edge(content::API* engine, typename MT::vector3_type const& position,
                     typename MT::quaternion_type const& orientation, typename MT::real_type const& scene_size,
                     MaterialInfo<typename MT::real_type> const& mat_info, mesh_array::TetGenSettings tetset);

template <typename MT>
void make_internal_edge(content::API* engine, typename MT::vector3_type const& position,
                        typename MT::quaternion_type const& orientation, typename MT::real_type const& scene_size,
                        MaterialInfo<typename MT::real_type> const& mat_info, mesh_array::TetGenSettings tetset);

template <typename MT> MaterialInfo<typename MT::real_type> create_material_info(content::API* engine);


} //namespace procedural

#include <procedural_make_scene.h>

// PROCEDURAL_H
#endif
