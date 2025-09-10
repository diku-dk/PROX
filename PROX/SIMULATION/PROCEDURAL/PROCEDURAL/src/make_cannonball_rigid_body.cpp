#include <procedural_factory.h>

#include <eigenhelperall.h>

#include <util_string_helper.h>

#include <string>

namespace procedural
{
template <typename T>
size_t make_cannonball_rigid_body(content::API* engine, GeometryHandleEigen<T> const& ball, const EigenVector3<T>& xbf,
                                  const EigenQuaternion<T>& Qbf, const EigenVector3<T>& vbf, MaterialInfo<T> mat_info)
{

    T const density = get_material_density_eigen<T>(mat_info, "Cannonball");
    size_t const mid = get_material_id_eigen<T>(mat_info, "Cannonball");

    static size_t counter = 0u;

    std::string const body_name = "ball_" + util::to_string(counter++);

    size_t const rid = engine->create_rigid_body(body_name);

    engine->set_rigid_body_position(rid, xbf(0), xbf(1), xbf(2));

    engine->set_rigid_body_orientation(rid, Qbf.w(), Qbf.x(), Qbf.y(), Qbf.z());

    engine->set_rigid_body_velocity(rid, vbf(0), vbf(1), vbf(2));

    engine->set_rigid_body_mass(rid, ball.m_m * density);

    engine->set_rigid_body_inertia(rid, ball.m_Ixx * density, ball.m_Iyy * density, ball.m_Izz * density);

    engine->connect_collision_geometry(rid, ball.m_gid);

    engine->set_rigid_body_material(rid, mid);

    engine->set_rigid_body_fixed(rid, false);

    return rid;
}

using MTf = tiny::MathTypes<float>;

template size_t make_cannonball_rigid_body<float>(content::API* engine, GeometryHandleEigen<float> const& ball,
                                                  const EigenVector3<float>& xbf, const EigenQuaternion<float>& Qbf,
                                                  const EigenVector3<float>& vbf, MaterialInfo<float> mat_info);

} //namespace procedural
