#ifndef PROCEDURAL_FACTORY_H
#define PROCEDURAL_FACTORY_H

#include <procedural_types.h>
#include <procedural_create_geometry_handle_obj.h>
#include <procedural_create_geometry_handle_mesh.h>

#include <content.h>
#include <mesh_array.h>

#include <string>

namespace procedural
{

  template<typename MT>
  size_t create_rigid_body(  content::API *  engine
                           , typename MT::vector3_type const & Tb2w
                           , typename MT::quaternion_type const & Qb2w
                           , GeometryHandle<MT>  const & geometry
                           , size_t const & mid
                           , typename MT::real_type const & density
                           , bool const fixed = false
                           , std::string const visual_material = "Visualizer/stone"
                           );

  template <typename T>
  size_t create_rigid_body(content::API* engine, const EigenVector3<T>& Tb2w, const EigenQuaternion<T>& Qb2w,
                           GeometryHandleEigen<T> const& geometry, size_t const& mid, const T& density,
                           bool const fixed = false, std::string const material_name = "Visualizer/stone");

  template<typename T>
  GeometryHandleEigen<T> create_geometry_handle_convex(
                                                   content::API * engine
                                                   , std::vector<EigenVector3<T>> const &  vertices
                                                   );

  template<typename MT>
  GeometryHandle<MT> create_geometry_handle_convex_old(
        content::API * engine
        , std::vector<typename MT::vector3_type> const &  vertices
        );

  template <typename T> GeometryHandleEigen<T> create_geometry_handle_sphere(content::API* engine, const T& radius);

  template<typename MT>
  GeometryHandle<MT> create_geometry_handle_capsule(
                                                    content::API * engine
                                                    , typename MT::real_type const & radius
                                                    , typename MT::real_type const & height
                                                    );

  template<typename MT>
  GeometryHandle<MT> create_geometry_handle_box(
                                                content::API * engine
                                                , typename MT::real_type  const & width
                                                , typename MT::real_type  const & height
                                                , typename MT::real_type  const & depth
                                                );

  template <typename T>
  GeometryHandleEigen<T> create_geometry_handle_box(content::API* engine, const T& width, const T& height,
                                                    const T& depth);

  template<typename T>
  GeometryHandleEigen<T> create_geometry_handle_cuboid(
                                                   content::API * engine
                                                   , EigenVector3<T>* vertices
                                                   );

  template<typename T>
  GeometryHandleEigen<T> create_geometry_handle_tetrahedron(
      content::API * engine
      , EigenVector3<T> one
      , EigenVector3<T> two
      , EigenVector3<T> three
      , EigenVector3<T> four
      );

  template<typename T>
  GeometryHandleEigen<T> create_geometry_handle_pillar_segment( content::API * engine
                                                               , const T bottom_radius
                                                               , const T top_radius
                                                               , const T height
                                                               , size_t const & slices
                                                               , mesh_array::TetGenSettings tetset
                                                               );

  template<typename MT>
  void compute_arch_stone_vertices(
                                   typename MT::real_type const & theta
                                   , typename MT::real_type const & depth
                                   , typename MT::real_type const & r_outer
                                   , typename MT::real_type const & r_inner
                                   , typename MT::vector3_type * vertices
                                   );

  template<typename T>
  void compute_arch_stone_vertices_eigen(
      const T& theta
      , const T& depth
      , const T& r_outer
      , const T& r_inner
      , EigenVector3<T>* vertices
      );
  template<typename MT>
  void compute_body_to_world_transform(
                                   typename MT::vector3_type      const & Tb2m
                                   , typename MT::quaternion_type const & Qb2m
                                   , typename MT::vector3_type    const & Tm2l
                                   , typename MT::quaternion_type const & Qm2l
                                   , typename MT::vector3_type    const & Tl2w
                                   , typename MT::quaternion_type const & Ql2w
                                   , typename MT::vector3_type    & Tb2w
                                   , typename MT::quaternion_type & Qb2w
                                   );

  template <typename T>
  void compute_body_to_world_transform(const EigenVector3<T>& Tb2m, const EigenQuaternion<T>& Qb2m,
                                       const EigenVector3<T>& Tm2l, const EigenQuaternion<T>& Qm2l,
                                       const EigenVector3<T>& Tl2w, const EigenQuaternion<T>& Ql2w,
                                       EigenVector3<T>& Tb2w, EigenQuaternion<T>& Qb2w);

  template <typename MT> size_t get_material_id(MaterialInfo<typename MT::real_type> info, std::string const material);

  template<typename MT>
  typename MT::real_type get_material_density(
                                              MaterialInfo<typename MT::real_type> info
                                              , std::string const material
                                              );

  template <typename T> size_t get_material_id_eigen(MaterialInfo<T> info, std::string const material);

  template <typename T> T get_material_density_eigen(MaterialInfo<T> info, std::string const material);

} // end of namespace procedural


// PROCEDURAL_FACTORY_H
#endif
