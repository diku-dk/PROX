#ifndef RIGID_BODY_GUI_CHALK_WRITE_GRAINS_H
#define RIGID_BODY_GUI_CHALK_WRITE_GRAINS_H

#include <content.h>

#include <mesh_array_t3mesh.h>
#include <mesh_array_vertex_attribute.h>

#include <tiny_math_types.h>
#include <tiny_coordsys_functions.h>

#include <util_log.h>

#include <fstream>

namespace rigid_body
{
  namespace gui
  {
    namespace chalk
    {

      inline void write_grains(
                               std::string const & filename
                               , content::API * engine
                               )
      {
        typedef tiny::MathTypes<float> MT;
        typedef MT::vector3_type       V;
        typedef MT::quaternion_type    Q;
        typedef MT::value_traits       VT;
        typedef MT::coordsys_type      C;

        std::ofstream file( filename.c_str() );

        if(!file.is_open())
        {
          util::Log logging;


          logging << "Could not open file = " << filename << util::Log::newline();

          return;
        }

        float x  = VT::zero();
        float y  = VT::zero();
        float z  = VT::zero();
        float qs = VT::zero();
        float qx = VT::zero();
        float qy = VT::zero();
        float qz = VT::zero();

        std::vector<size_t> rids;
        rids.resize(engine->get_number_of_rigid_bodies());

        engine->get_rigid_body_indices( &rids[0] );

        for(size_t i = 0u; i < rids.size(); ++i)
        {
          size_t const rid = rids[i];

          engine->get_rigid_body_position(rid,x,y,z);
          engine->get_rigid_body_orientation(rid,qs,qx,qy,qz);

          V const T_b2w = V::make(x, y, z);
          Q const Q_b2w = Q(qs, qx, qy, qz);

          size_t const gid                = engine->get_rigid_body_collision_geometry(rid);
          size_t const number_of_convexes = engine->get_number_of_convexes(gid);

          if(number_of_convexes!=1)  // only support for one convex (=grain) shape per body
            continue;

          size_t const convex_number = 0u;

          engine->get_convex_position( gid, convex_number , x, y, z );
          engine->get_convex_orientation( gid, convex_number, qs, qx, qy, qz );

          V const T_s2b = V::make(x, y, z);
          Q const Q_s2b = Q(qs, qx, qy, qz);

          size_t no_points = 0u;

          engine->get_convex_shape( gid, convex_number , no_points );

          std::vector<float> coordinates;
          coordinates.resize(no_points*3u);

          engine->get_convex_shape( gid, convex_number, &coordinates[0] );

          C const X_s2b = C::make(T_s2b,Q_s2b);
          C const X_b2w = C::make(T_b2w,Q_b2w);
          C const X_s2w = tiny::prod(X_b2w, X_s2b);

          {
            util::Log logging;

            logging << "Grain " << i << util::Log::newline();
          }

          for(size_t k=0u;k<no_points;++k)
          {
            V const p = V::make( coordinates[3u*k], coordinates[3u*k+1u], coordinates[3u*k+2u]);
            V const q = tiny::xform_point(X_s2w, p);//changed from X_s2b
            file << "(" << q(0) << "," << q(1) << "," << q(2) << ")";
            
            if(k < no_points-1)
              file << ",";
          }
          file << "\n";
        }
        file.flush();
        file.close();
      }
      
      
    }//namespace chalk
  }//namespace gui
}//namespace rigid_body

// RIGID_BODY_GUI_CHALK_WRITE_GRAINS_H
#endif