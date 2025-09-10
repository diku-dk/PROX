#ifndef RIGID_BODY_GUI_CHALK_WRITE_GRAINS_H
#define RIGID_BODY_GUI_CHALK_WRITE_GRAINS_H

#include <content.h>

#include <mesh_array_t3mesh.h>
#include <mesh_array_vertex_attribute.h>

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
          using T = float;
          std::ofstream file(filename.c_str());

          if (!file.is_open())
          {
              util::Log logging;

              logging << "Could not open file = " << filename
                      << util::Log::newline();

              return;
          }

        float x  = 0;
        float y  = 0;
        float z  = 0;
        float qs = 0;
        float qx = 0;
        float qy = 0;
        float qz = 0;

        std::vector<size_t> rids;
        rids.resize(engine->get_number_of_rigid_bodies());

        engine->get_rigid_body_indices( &rids[0] );

        for(size_t i = 0u; i < rids.size(); ++i)
        {
          size_t const rid = rids[i];

          engine->get_rigid_body_position(rid,x,y,z);
          engine->get_rigid_body_orientation(rid,qs,qx,qy,qz);

          const EigenVector3<T> T_b2w = EigenVector3<T>(x, y, z);
          const EigenQuaternion<T> Q_b2w = EigenQuaternion<T>(qs, qx, qy, qz);

          size_t const gid                = engine->get_rigid_body_collision_geometry(rid);
          size_t const number_of_convexes = engine->get_number_of_convexes(gid);

          if(number_of_convexes!=1)  // only support for one convex (=grain) shape per body
            continue;

          size_t const convex_number = 0u;

          engine->get_convex_position( gid, convex_number , x, y, z );
          engine->get_convex_orientation( gid, convex_number, qs, qx, qy, qz );

          const EigenVector3<T> T_s2b = EigenVector3<T>(x, y, z);
          const EigenQuaternion<T> Q_s2b = EigenQuaternion<T>(qs, qx, qy, qz);

          size_t no_points = 0u;

          engine->get_convex_shape( gid, convex_number , no_points );

          std::vector<float> coordinates;
          coordinates.resize(no_points*3u);

          engine->get_convex_shape( gid, convex_number, &coordinates[0] );

          CoordSysEigen<T> const X_s2b = CoordSysEigen<T>::make(T_s2b, Q_s2b);
          CoordSysEigen<T> const X_b2w = CoordSysEigen<T>::make(T_b2w, Q_b2w);
          CoordSysEigen<T> const X_s2w = prod(X_b2w, X_s2b);

          {
            util::Log logging;

            logging << "Grain " << i << util::Log::newline();
          }

          for(size_t k=0u;k<no_points;++k)
          {
              const EigenVector3<T> p = EigenVector3<T>(
                  coordinates[3u * k], coordinates[3u * k + 1u],
                  coordinates[3u * k + 2u]);
              const EigenVector3<T> q
                  = xform_point(X_s2w, p); //changed from X_s2b
              file << "(" << q(0) << "," << q(1) << "," << q(2) << ")";

              if (k < no_points - 1) file << ",";
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
