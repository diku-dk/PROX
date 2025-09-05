#include <procedural.h>
#include <procedural_factory.h>

#include <tiny_math_types.h>

namespace procedural
{
template <typename T>
void make_slide(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
                const T& degree, MaterialInfo<T> mat_info)
{
    size_t const mid = get_material_id_eigen<T>(mat_info, "Stone");
    size_t const stone_density = get_material_density_eigen<T>(mat_info, "Stone");

    T const radians = degree * (std::numbers::pi_v<T> / 180);
    T const radius = 0.5f;
    const EigenVector3<T> hill_extents = EigenVector3<T>((6.0), (0.2), 4);
    const EigenVector3<T> box_extents = EigenVector3<T>(radius, radius, radius) * 2;

    GeometryHandleEigen<T> hill
        = create_geometry_handle_box<T>(engine, hill_extents(0), hill_extents(1), hill_extents(2));
    GeometryHandleEigen<T> stone
        = create_geometry_handle_box<T>(engine, box_extents(0), box_extents(1), box_extents(2));
    //  GeometryHandle<MT> sphere = create_geometry_handle_sphere<MT>( engine, radius);

    //BF or MF here?
    EigenVector3<T> Tm = rotate(
        Rotatey<T>(radians),
        EigenVector3<T>(
            hill_extents(0) * 0.5f, 0,
            0));  //V::make(  std::cos(radians)*hill_extents(0)*0.5f, std::sin(radians)*hill_extents(0)*0.5f, 0);
    EigenQuaternion<T> Qm = Rotatez(radians);

    // Body to Model transform
    EigenVector3<T> Tb = hill.Tb2m();
    EigenQuaternion<T> Qb = hill.Qb2m();

    // Body to World transform
    EigenVector3<T> Tw = rotate(Qm, Tb) + Tm;
    EigenQuaternion<T> Qw = Qm * Qb;

    // Apply any user transforms
    EigenVector3<T> Tu = rotate(orientation, Tw) + position;
    EigenQuaternion<T> Qu = orientation * Qw;

    create_rigid_body<T>(engine, Tu, Qu, hill, mid, stone_density, true);

    EigenVector3<T> TmB = rotate(
        Rotatey<T>(radians), EigenVector3<T>(hill_extents(0) - box_extents(0), (box_extents(1) + hill_extents(1)), 0));
    //V TmB = V::make(  std::cos(radians)*hill_extents(0)*0.5f-box_extents(0)*0.5f
      //              , std::sin(radians)*hill_extents(0)*0.5f+box_extents(1)*0.5f
      //              , radius*2);
   // TmB   = TmB + Tm;

    // Body to Model transform
    Tb = stone.Tb2m();
    Qb = stone.Qb2m();

    // Body to World transform
    Tw = rotate(Qm, Tb) + TmB;
    Qw = Qm * Qb;

    // Apply any user transforms
    Tu = rotate(orientation, Tw) + position;
    Qu = orientation * Qw;

    create_rigid_body<T>(engine, Tu, Qu, stone, mid, stone_density);
    /*
     V TmS = V::make(  std::cos(radians)*hill_extents(0)*0.5f-radius
     , std::sin(radians)*hill_extents(0)*0.5f+radius
     , -radius*2);
     TmS   = TmS + Tm;

     // Body to Model transform
     Tb = sphere.Tbf();
     Qb = sphere.Qbf();

     // Body to World transform
     Tw = rotate(Qm, Tb) + TmS;
     Qw = Qm*Qb;

     // Apply any user transforms
     Tu = rotate(orientation, Tw) + position;
     Qu = orientation*Qw;

     create_rigid_body<MT>(  engine
     , Tu
     , Qu
     , sphere
     , mid
     , stone_density
     );
     */
}

template void make_slide<float>(content::API* engine, const EigenVector3<float>& position,
                                const EigenQuaternion<float>& orientation, const float& degree,
                                MaterialInfo<float> mat_info);

} // namespace procedural

//namespace procedural
