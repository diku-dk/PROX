#include <procedural.h>
#include <procedural_factory.h>

#include <tiny_math_types.h>

namespace procedural
{

template <typename T>
void make_cannonball(content::API* engine, const T& radius, const EigenVector3<T>& position,
                     const EigenQuaternion<T>& orientation, const EigenVector3<T>& direction, MaterialInfo<T> mat_info)

{
    GeometryHandleEigen<T> cannon_ball = make_cannonball_geometry<T>(engine, radius);
    make_cannonball_rigid_body(engine, cannon_ball, position, orientation, direction, mat_info);
}

using MTf = tiny::MathTypes<float>;

template void make_cannonball<float>(content::API* engine, const float& radius, const EigenVector3<float>& position,
                                     const EigenQuaternion<float>& orientation, const EigenVector3<float>& direction,
                                     MaterialInfo<float> mat_info);

} // namespace procedural
