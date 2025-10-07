#include <procedural_factory.h>
#include <procedural_factory_types.h>
#include <procedural_types.h>
#include <procedural_noise.h>

#include <content.h>

#include <mesh_array.h>

#include <eigenhelperall.h>
#include <mass.h>
#include <util_string_helper.h>

#include <string>

namespace procedural
{

template <typename T>
GeometryHandleEigen<T> create_geometry_handle_convex(content::API* engine, std::vector<EigenVector3<T>> const& vertices)
{

    static size_t counter = 0u;

    std::string const geom_name = "convex_" + util::to_string(counter++);

    detail::MeshData<T> data;

    //fix vertex sequence
    mesh_array::make_convex<T>(vertices, data.m_mesh, data.m_X, data.m_Y,
                               data.m_Z);

    mass::Properties<T> props_mf = mass::compute_mesh<T>(1, data.m_mesh.triangle_size(), &data);
    mass::Properties<T> props_bf = mass::translate_to_body_frame(props_mf);
    mass::Properties<T> props = mass::rotate_to_body_frame(props_bf);

    // Change geometry from model space to body space
    const EigenVector3<T> d = -EigenVector3<T>(props_mf.m_x, props_mf.m_y, props_mf.m_z);
    const EigenQuaternion<T> R = EigenQuaternion<T>(props.m_Qs, props.m_Qx, props.m_Qy, props.m_Qz);
    mesh_array::translate<T>(d, data.m_mesh, data.m_X, data.m_Y, data.m_Z);
    mesh_array::rotate<T>((R).conjugate(), data.m_mesh, data.m_X, data.m_Y, data.m_Z);

    size_t const gid = engine->create_collision_geometry(geom_name);
    size_t const sid = engine->create_convex_shape(gid);
    size_t const N = data.m_mesh.vertex_size();

    std::vector<T> coords;
    coords.resize(3u * N);

    for (size_t i = 0; i < N; ++i)
    {
        mesh_array::Vertex const v = data.m_mesh.vertex(i);
        coords[3 * i] = data.m_X(v);
        coords[3 * i + 1] = data.m_Y(v);
        coords[3 * i + 2] = data.m_Z(v);
    }

    engine->set_convex_shape(gid, sid, N, &coords[0]);

    return GeometryHandleEigen<T>(props.m_m, props.m_Ixx, props.m_Iyy, props.m_Izz,
                                  EigenVector3<T>(props_mf.m_x, props_mf.m_y, props_mf.m_z),
                                  EigenQuaternion<T>(props.m_Qs, props.m_Qx, props.m_Qy, props.m_Qz), gid);
}

template GeometryHandleEigen<float> create_geometry_handle_convex(content::API* engine,
                                                                  std::vector<EigenVector3<float>> const& vertices);

template <typename T> GeometryHandleEigen<T> create_geometry_handle_sphere(content::API* engine, const T& radius)
{
    static size_t counter = 0u;

    std::string const geom_name = "sphere_" + util::to_string(counter++);

    mass::Properties<T> props = mass::compute_sphere<T>(1, radius);

    size_t const gid = engine->create_collision_geometry(geom_name);
    size_t const sid = engine->create_sphere_shape(gid);

    engine->set_sphere_shape(gid, sid, radius);

    return GeometryHandleEigen<T>(props.m_m, props.m_Ixx, props.m_Iyy, props.m_Izz, EigenVector3<T>(0, 0, 0),
                                  EigenQuaternion<T>::Identity(), gid);
}

template GeometryHandleEigen<float> create_geometry_handle_sphere<float>(content::API* engine, const float& radius);

template <typename T>
GeometryHandleEigen<T> create_geometry_handle_capsule(content::API* engine, const float& radius, const float& height)
{

    static size_t counter = 0u;

    std::string const geom_name = "capsule_" + util::to_string(counter++);

    mass::Properties<T> props = mass::compute_capsule<T>(1, radius, 0.5f * height);

    size_t const gid = engine->create_collision_geometry(geom_name);
    size_t const sid = engine->create_capsule_shape(gid);

    engine->set_capsule_shape(gid, sid, radius, height);

    return GeometryHandleEigen<T>(props.m_m, props.m_Ixx, props.m_Iyy, props.m_Izz, EigenVector3<T>(0, 0, 0),
                                  EigenQuaternion<T>::Identity(), gid);
}

template GeometryHandleEigen<float> create_geometry_handle_capsule<float>(content::API* engine, const float& radius,
                                                                          const float& height);

template <typename T>
GeometryHandleEigen<T> create_geometry_handle_box(content::API* engine, const T& width, const T& height, const T& depth)
{

    static size_t counter = 0u;

    std::string const geom_name = "box_" + util::to_string(counter++);

    mass::Properties<T> props = mass::compute_box<T>(1, width / 2, height / 2, depth / 2);

    size_t const gid = engine->create_collision_geometry(geom_name);
    size_t const sid = engine->create_box_shape(gid);

    engine->set_box_shape(gid, sid, width, height, depth);

    return GeometryHandleEigen<T>(props.m_m, props.m_Ixx, props.m_Iyy, props.m_Izz, EigenVector3<T>(0, 0, 0),
                                  EigenQuaternion<T>::Identity(), gid);
}

template GeometryHandleEigen<float> create_geometry_handle_box<float>(content::API* engine, const float& width,
                                                                      const float& height, const float& depth);

template <typename T>
GeometryHandleEigen<T> create_geometry_handle_tetrahedron(content::API* engine, EigenVector3<T> one,
                                                          EigenVector3<T> two, EigenVector3<T> three,
                                                          EigenVector3<T> four)
{

    static size_t counter = 0u;

    std::string const geom_name = "tetra_" + util::to_string(counter++);

    detail::MeshData<T> data;

    mesh_array::make_tetrahedron<T>(one, two, three, four, data.m_mesh, data.m_X, data.m_Y, data.m_Z);

    mass::Properties<T> props_mf = mass::compute_mesh<T>(1, data.m_mesh.triangle_size(), &data);
    mass::Properties<T> props_bf = mass::translate_to_body_frame(props_mf);
    mass::Properties<T> props = mass::rotate_to_body_frame(props_bf);

      // Change geometry from model space to body space
    const EigenVector3<T> d = -EigenVector3<T>(props_mf.m_x, props_mf.m_y, props_mf.m_z);
    const EigenQuaternion<T> R = EigenQuaternion<T>(props.m_Qs, props.m_Qx, props.m_Qy, props.m_Qz);
    mesh_array::translate<T>(d, data.m_mesh, data.m_X, data.m_Y, data.m_Z);
    mesh_array::rotate<T>((R).conjugate(), data.m_mesh, data.m_X, data.m_Y, data.m_Z);

    size_t const gid = engine->create_collision_geometry(geom_name);
    size_t const sid = engine->create_convex_shape(gid);
    size_t const N = data.m_mesh.vertex_size();

    std::vector<T> coords;
    coords.resize(3u * N);

    for (size_t i = 0; i < N; ++i)
    {
        mesh_array::Vertex const v = data.m_mesh.vertex(i);
        coords[3 * i] = data.m_X(v);
        coords[3 * i + 1] = data.m_Y(v);
        coords[3 * i + 2] = data.m_Z(v);
    }

    engine->set_convex_shape(gid, sid, N, &coords[0]);

    return GeometryHandleEigen<T>(props.m_m, props.m_Ixx, props.m_Iyy, props.m_Izz,
                                  EigenVector3<T>(props_mf.m_x, props_mf.m_y, props_mf.m_z),
                                  EigenQuaternion<T>(props.m_Qs, props.m_Qx, props.m_Qy, props.m_Qz), gid);
}

template GeometryHandleEigen<float>
create_geometry_handle_tetrahedron<float>(content::API* engine, EigenVector3<float> one, EigenVector3<float> two,
                                          EigenVector3<float> three, EigenVector3<float> four);

template <typename T>
GeometryHandleEigen<T> create_geometry_handle_cuboid(content::API* engine, EigenVector3<T>* vertices)
{
    static size_t counter = 0u;

    std::string const geom_name = "cuboid_" + util::to_string(counter++);

    detail::MeshData<T> data;

    //fix vertex sequence
    mesh_array::make_cuboid<T>(vertices[0], vertices[1], vertices[2], vertices[3], vertices[4], vertices[5],
                               vertices[6], vertices[7], data.m_mesh, data.m_X, data.m_Y, data.m_Z);

    mass::Properties<T> props_mf = mass::compute_mesh<T>(1, data.m_mesh.triangle_size(), &data);
    mass::Properties<T> props_bf = mass::translate_to_body_frame(props_mf);
    mass::Properties<T> props = mass::rotate_to_body_frame(props_bf);

      // Change geometry from model space to body space
    const EigenVector3<T> d = -EigenVector3<T>(props_mf.m_x, props_mf.m_y, props_mf.m_z);
    const EigenQuaternion<T> R = EigenQuaternion<T>(props.m_Qs, props.m_Qx, props.m_Qy, props.m_Qz);
    mesh_array::translate<T>(d, data.m_mesh, data.m_X, data.m_Y, data.m_Z);
    mesh_array::rotate<T>((R).conjugate(), data.m_mesh, data.m_X, data.m_Y, data.m_Z);

    size_t const gid = engine->create_collision_geometry(geom_name);
    size_t const sid = engine->create_convex_shape(gid);
    size_t const N = data.m_mesh.vertex_size();

    std::vector<T> coords;
    coords.resize(3u * N);

    for (size_t i = 0; i < N; ++i)
    {
        mesh_array::Vertex const v = data.m_mesh.vertex(i);
        coords[3 * i] = data.m_X(v);
        coords[3 * i + 1] = data.m_Y(v);
        coords[3 * i + 2] = data.m_Z(v);
    }

    engine->set_convex_shape(gid, sid, N, &coords[0]);

    return GeometryHandleEigen<T>(props.m_m, props.m_Ixx, props.m_Iyy, props.m_Izz,
                                  EigenVector3<T>(props_mf.m_x, props_mf.m_y, props_mf.m_z),
                                  EigenQuaternion<T>(props.m_Qs, props.m_Qx, props.m_Qy, props.m_Qz), gid);
}

template GeometryHandleEigen<float> create_geometry_handle_cuboid<float>(content::API* engine,
                                                                         EigenVector3<float>* vertices);

template <typename T>
GeometryHandleEigen<T> create_geometry_handle_pillar_segment(content::API* engine, const T bottom_radius,
                                                             const T top_radius, const T height, size_t const& slices,
                                                             mesh_array::TetGenSettings tetset)
{
    static size_t counter = 0u;

    std::string const geom_name = "pillar_segment_" + util::to_string(counter++);

    mesh_array::T3Mesh mesh;
    mesh_array::VertexAttribute<T, mesh_array::T3Mesh> X;
    mesh_array::VertexAttribute<T, mesh_array::T3Mesh> Y;
    mesh_array::VertexAttribute<T, mesh_array::T3Mesh> Z;

    mesh_array::make_conical<T>(bottom_radius, top_radius, height, slices, mesh, X, Y, Z);

    mass::Properties<T> props_mf = mass::compute_conical_solid<T>(1, bottom_radius, top_radius, height);
    mass::Properties<T> props_bf = mass::translate_to_body_frame(props_mf);
    mass::Properties<T> props = mass::rotate_to_body_frame(props_bf);

    // Change geometry from model space to body space
    const EigenVector3<T> d = -EigenVector3<T>(props_mf.m_x, props_mf.m_y, props_mf.m_z);
    const EigenQuaternion<T> R = EigenQuaternion<T>(props.m_Qs, props.m_Qx, props.m_Qy, props.m_Qz);
    mesh_array::translate<T>(d, mesh, X, Y, Z);
    mesh_array::rotate<T>((R).conjugate(), mesh, X, Y, Z);

    mesh_array::T4Mesh volume;
    mesh_array::VertexAttribute<T, mesh_array::T4Mesh> volX;
    mesh_array::VertexAttribute<T, mesh_array::T4Mesh> volY;
    mesh_array::VertexAttribute<T, mesh_array::T4Mesh> volZ;

    mesh_array::tetgen(mesh, X, Y, Z, volume, volX, volY, volZ, tetset);

    size_t const gid = engine->create_collision_geometry(geom_name);
    size_t const N = volume.vertex_size();
    size_t const K = volume.tetrahedron_size();

    std::vector<size_t> verts;
    verts.resize(N);

    for (size_t i = 0; i < N; ++i) { verts[i] = i; }

    std::vector<T> coords;
    coords.resize(3u * N);

    for (size_t i = 0; i < N; ++i)
    {
        mesh_array::Vertex const v = volume.vertex(i);
        coords[3 * i] = volX(v);
        coords[3 * i + 1] = volY(v);
        coords[3 * i + 2] = volZ(v);
    }

    std::vector<size_t> tets;
    tets.resize(4u * K);

    for (size_t i = 0; i < K; ++i)
    {
        mesh_array::Tetrahedron const t = volume.tetrahedron(i);
        tets[4 * i] = t.i();
        tets[4 * i + 1] = t.j();
        tets[4 * i + 2] = t.k();
        tets[4 * i + 3] = t.m();
    }
    engine->set_tetramesh_shape(gid, N, K, &verts[0], &tets[0], &coords[0]);

    return GeometryHandleEigen<T>(props.m_m, props.m_Ixx, props.m_Iyy, props.m_Izz,
                                  EigenVector3<T>(props_mf.m_x, props_mf.m_y, props_mf.m_z),
                                  EigenQuaternion<T>(props.m_Qs, props.m_Qx, props.m_Qy, props.m_Qz), gid);
}

template GeometryHandleEigen<float>
create_geometry_handle_pillar_segment<float>(content::API* engine, const float bottom_radius, const float top_radius,
                                             const float height, size_t const& slices,
                                             mesh_array::TetGenSettings tetset);

template <typename T>
size_t create_rigid_body(content::API* engine, const EigenVector3<T>& Tb2w, const EigenQuaternion<T>& Qb2w,
                         GeometryHandleEigen<T> const& geometry, size_t const& mid, const T& density, bool const fixed,
                         std::string const material_name)
{

    static size_t counter = 0u;

    std::string const body_name = "body_" + util::to_string(counter++);

    size_t const rid = engine->create_rigid_body(body_name);

      //2011-05-04 Mort: should these be passed as a parameter? Right now I have just hardcoded them to zero;
      //set_rigid_body spin?
    engine->set_rigid_body_velocity(rid, 0, 0, 0);

    if (Noise::on())
    {
        const EigenVector3<T> noise = randomEigen<T>(-Noise::scale(), Noise::scale());

          //engine->set_rigid_body_position( rid, Tb2w(0)+ noise(0), Tb2w(1)+ noise(1), Tb2w(2)+ noise(2) );
        engine->set_rigid_body_position(rid, Tb2w(0) + noise(0), Tb2w(1), Tb2w(2) + noise(2));
    }
    else { engine->set_rigid_body_position(rid, Tb2w(0), Tb2w(1), Tb2w(2)); }

    engine->set_rigid_body_orientation(rid, Qb2w.w(), Qb2w.x(), Qb2w.y(), Qb2w.z());

    engine->set_rigid_body_mass(rid, geometry.m_m * density);

    engine->set_rigid_body_inertia(rid, geometry.m_Ixx * density, geometry.m_Iyy * density, geometry.m_Izz * density);

    engine->connect_collision_geometry(rid, geometry.m_gid);

    engine->set_rigid_body_material(rid, mid);

    engine->set_rigid_body_fixed(rid, fixed);

    return rid;
}

template size_t create_rigid_body<float>(content::API* engine, const EigenVector3<float>& Tb2w,
                                         const EigenQuaternion<float>& Qb2w, GeometryHandleEigen<float> const& geometry,
                                         size_t const& mid, const float& density, bool const fixed,
                                         std::string const material_name);

template <typename T>
void compute_arch_stone_vertices_eigen(const T& theta, const T& depth, const T& r_outer, const T& r_inner,
                                       EigenVector3<T>* vertices)
{

    using std::cos;
    using std::sin;

    T const center_height = (r_outer + r_inner) * 0.5f;
    T const max_height = r_outer * cos(theta * 0.5f);
    T const min_height = r_inner * cos(theta * 0.5f);
    T const half_width_inner = r_inner * sin(theta * 0.5f);
    T const half_width_outer = r_outer * sin(theta * 0.5f);
    T const half_depth = depth * 0.5f;

    vertices[0] = EigenVector3<T>(-half_width_inner, min_height - center_height, half_depth);
    vertices[1] = EigenVector3<T>(half_width_inner, min_height - center_height, half_depth);
    vertices[2] = EigenVector3<T>(half_width_outer, max_height - center_height, half_depth);
    vertices[3] = EigenVector3<T>(-half_width_outer, max_height - center_height, half_depth);

    vertices[4] = vertices[0] - EigenVector3<T>(0, 0, depth);
    vertices[5] = vertices[1] - EigenVector3<T>(0, 0, depth);
    vertices[6] = vertices[2] - EigenVector3<T>(0, 0, depth);
    vertices[7] = vertices[3] - EigenVector3<T>(0, 0, depth);
}

template void compute_arch_stone_vertices_eigen<float>(const float& theta, const float& depth, const float& r_outer,
                                                       const float& r_inner, EigenVector3<float>* vertices);


template <typename T>
void compute_body_to_world_transform(const EigenVector3<T>& T_body2model, const EigenQuaternion<T>& Q_body2model,
                                     const EigenVector3<T>& T_model2local, const EigenQuaternion<T>& Q_model2local,
                                     const EigenVector3<T>& T_local2world, const EigenQuaternion<T>& Q_local2world,
                                     EigenVector3<T>& T_body2world, EigenQuaternion<T>& Q_body2world)
{

    const CoordSysEigen<T> Xb2m = CoordSysEigen<T>::make(T_body2model, Q_body2model);
    const CoordSysEigen<T> Xm2l = CoordSysEigen<T>::make(T_model2local, Q_model2local);
    const CoordSysEigen<T> Xl2w = CoordSysEigen<T>::make(T_local2world, Q_local2world);

    const CoordSysEigen<T> Xb2l = prod(Xm2l, Xb2m);
    const CoordSysEigen<T> Xb2w = prod(Xl2w, Xb2l);

    Q_body2world = Xb2w.Q();
    T_body2world = Xb2w.T();
}


template void compute_body_to_world_transform<float>(
    const EigenVector3<float>& T_body2model, const EigenQuaternion<float>& Q_body2model,
    const EigenVector3<float>& T_model2local, const EigenQuaternion<float>& Q_model2local,
    const EigenVector3<float>& T_local2world, const EigenQuaternion<float>& Q_local2world,
    EigenVector3<float>& T_body2world, EigenQuaternion<float>& Q_body2world);

template <typename T> size_t get_material_id_eigen(MaterialInfo<T> info, std::string const material)
{
    if (material.compare("Stone") == 0) return info.m_stone_mid;

    if (material.compare("Ground") == 0) return info.m_ground_mid;

    if (material.compare("Cannonball") == 0) return info.m_cannonball_mid;

    return 0u;
}

template size_t get_material_id_eigen<float>(MaterialInfo<float> info, std::string const material);

template <typename T> T get_material_density_eigen(MaterialInfo<T> info, std::string const material)
{
    if (material.compare("Stone") == 0) return info.m_stone_density;

    if (material.compare("Ground") == 0) return info.m_ground_density;

    if (material.compare("Cannonball") == 0) return info.m_cannonball_density;

    return 0;
}

template float get_material_density_eigen<float>(MaterialInfo<float> info, std::string const material);

} //namespace procedural
