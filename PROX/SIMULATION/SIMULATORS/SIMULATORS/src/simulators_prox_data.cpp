#include <simulators_prox_data.h>

#include <prox_compute_raycast.h>
#include <prox_update_inertia_tensor.h>

#include <util_profiling.h>
#include <util_string_helper.h>
#include <util_log.h>

#include <prox_rigid_body.h>

#include <cassert>

#include <steppers/prox_time_stepper.h>
#include <grid_helpers.h>

namespace simulators
{

ProxData::ProxData()
    : m_geometry_names()
    , m_materials()
    , m_bodies()
    , m_contacts()
    , m_broad()
    , m_narrow()
    , m_property_counter(0u)
    , m_time_step(0.01f)
    , m_time(0.0f)
    , m_params()
    , m_use_only_tetrameshes(false)
    , m_tetgen_settings(mesh_array::tetgen_quality_settings())
    , m_all_scripted_bodies()
{
    clear();
}

ProxData::~ProxData() { clear(); }

void ProxData::clear()
{
    m_bodies.clear();
    m_contacts.clear();
    m_broad.clear();
    m_narrow.clear();
    m_geometry_names.clear();
    m_materials.clear();
    m_property_counter = 0u;
    m_time = 0.0f;
    m_time_step = 0.01f;
    m_params = params_type();
    m_use_only_tetrameshes = false;

    m_properties.resize(m_number_of_materials);

    for (size_t i = 0u; i < m_number_of_materials; ++i)
    {
        m_properties[i].resize(m_number_of_materials);

        for (size_t j = 0u; j < m_number_of_materials; ++j) { m_exist_property[i][j] = false; }
    }

    m_gravity = {};
    m_damping = {};

    m_force_callbacks.clear();
    m_pin_forces.clear();

    m_all_scripted_bodies.clear();
    m_motion_callbacks.clear();
    m_keyframe_motions.clear();
    m_oscillation_motions.clear();

    util::Profiling::reset();
}

void ProxData::step_simulation(float const& dt)
{
    assert(dt > 0.0f || !"step_simulation(): invalid step size");
    assert(dt <= m_time_step || !"step_simulation(): invalid step size");

    prox::time_stepper(dt, m_bodies, m_properties, m_gravity, m_damping,
                       m_params, m_broad, m_narrow, m_contacts);

    T E_kinetic;
    T E_potential;
    get_total_energy(E_kinetic, E_potential);

    RECORD("dt", dt);
    RECORD("Ekin", E_kinetic);
    RECORD("Epot", E_potential);
}

bool ProxData::compute_raycast(float const& p_x, float const& p_y, float const& p_z, float const& ray_x,
                               float const& ray_y, float const& ray_z, size_t& body_index, float& hit_x, float& hit_y,
                               float& hit_z, float& hit_distance)
{
    auto ray = geometry::make_ray(EigenVector3<float>(p_x, p_y, p_z),
                                  EigenVector3<float>(ray_x, ray_y, ray_z));

    EigenVector3<T> point;
    T distance;

    bool const did_hit = prox::compute_raycast<float>(
        ray, m_bodies, m_narrow, body_index, point, distance);

    if (did_hit)
    {
        hit_x = point(0);
        hit_y = point(1);
        hit_z = point(2);
        hit_distance = distance;
    }

    return did_hit;
}

//Inconsistent call, but after doing setTetrameshShape, also generate a sdf?
void ProxData::make_tetramesh_geoemtry(geometry_type& geometry, mesh_array::T3Mesh const& surface_in,
                                       mesh_array::VertexAttribute<T, mesh_array::T3Mesh> const& surface_X_in,
                                       mesh_array::VertexAttribute<T, mesh_array::T3Mesh> const& surface_Y_in,
                                       mesh_array::VertexAttribute<T, mesh_array::T3Mesh> const& surface_Z_in)
{
    mesh_array::T3Mesh surface = surface_in;
    mesh_array::VertexAttribute<T, mesh_array::T3Mesh> surface_X = surface_X_in;
    mesh_array::VertexAttribute<T, mesh_array::T3Mesh> surface_Y = surface_Y_in;
    mesh_array::VertexAttribute<T, mesh_array::T3Mesh> surface_Z = surface_Z_in;

    if (kdop::SelectContactPointAlgorithm::is_using_closest_point())
    {
//      mesh_array::shrink<TT>(
//                             VT::numeric_cast(0.99)
//                             , surface
//                             , surface_X
//                             , surface_Y
//                             , surface_Z
//                             );
        std::vector<EigenVector3<T>> normals;

        mesh_array::compute_vertex_normals<T>(surface, surface_X, surface_Y, surface_Z, normals);

        T const distance = (-0.01);

        mesh_array::displace_vertices<T>(surface, surface_X, surface_Y, surface_Z, distance, normals);
    }

    mesh_array::T4Mesh volume;
    mesh_array::VertexAttribute<T, mesh_array::T4Mesh> volume_X;
    mesh_array::VertexAttribute<T, mesh_array::T4Mesh> volume_Y;
    mesh_array::VertexAttribute<T, mesh_array::T4Mesh> volume_Z;

    mesh_array::tetgen(surface, surface_X, surface_Y, surface_Z, volume, volume_X, volume_Y, volume_Z,
                       m_tetgen_settings);

    geometry.m_tetramesh.set_tetramesh_shape(volume, volume_X, volume_Y,
                                             volume_Z);
    grid::Grid<T, T> grid;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    grid::build_VF_from_T3Mesh(surface, surface_X, surface_Y, surface_Z, V, F);
    grid = grid::projectGridToSDF<T, T>(
        V, F, Eigen::Matrix<size_t, 3, 1>(32, 32, 32));
    geometry.m_signedDistanceMap.setSignedDistanceGrid(grid);
    grid::writeGridToFiles<T, T>(grid, 1);
/*    grid::Grid<T, T> grid;
    grid.create() geometry.m_signedDistanceMap =*/
//    geometry.add_shape()
}

void ProxData::get_total_energy(float& kinetic, float& potential)
{
    kinetic = 0.0f;
    potential = 0.0f;

    for (size_t i = 1u; i < m_bodies.size(); ++i)
    {
        rigid_body_type const& body = m_bodies[i];

        if (body.is_fixed()) continue;

      // 2016-06-21 Kenny code review: Hmm, should scripted bodies be
      //            ignore? They can move around and when colliding with
      //            free bodies they will transfer energy. Hence, they
      //            work like an "external" force on the system right now.
        if (body.is_scripted()) continue;

        float const m = body.get_mass();
        float const h = dot(m_gravity.up(), body.get_position());
        float const v = norm(body.get_velocity());

        auto w = body.get_spin();
        auto I_bf = body.get_inertia_bf();
        auto R = body.get_orientation().toRotationMatrix();

        EigenMatrix3<float> I;
        prox::detail::update_inertia_tensor(R, I_bf, I);

        T const wIw = (w).dot(I * w);

        kinetic += 0.5 * (m * v * v + wIw);
        potential += m_gravity.acceleration() * m * h;
    }
}

} //namespace simulators
