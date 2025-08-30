#include <convex_shapes.h>

#include <tiny_math_types.h>
#include <tiny_is_finite.h>

#include <cassert>
#include <cmath>     // needed for std::min and std::sqrt

namespace convex
{
  template<typename M>
  typename M::real_type const & Cone<M>::half_height() const { return this->m_half_height; }

  template<typename M>
  typename M::real_type       & Cone<M>::half_height()       { return this->m_half_height; }

  template<typename M>
  typename M::real_type const & Cone<M>::base_radius() const { return this->m_base_radius; }

  template<typename M>
    typename M::real_type       & Cone<M>::base_radius()       { return this->m_base_radius; }

  template<typename M>
  Cone<M>::Cone()
  : m_half_height( 1 )
  , m_base_radius( 1 )
  {}

  template<typename M>
  typename M::vector3_type Cone<M>::get_support_point(typename M::vector3_type const& v) const
  {
    using std::sqrt;

    typedef typename M::real_type     T;
    typedef typename M::vector3_type  V;
    typedef typename M::value_traits VT;

    T const & vx = v(0);
    T const & vy = v(1);
    T const & vz = v(2);

    assert( is_number(vx) || !"NAN encountered");
    assert( is_number(vy) || !"NAN encountered");
    assert( is_number(vz) || !"NAN encountered");
    assert( is_finite(vx) || !"INF encountered");
    assert( is_finite(vy) || !"INF encountered");
    assert( is_finite(vz) || !"INF encountered");

    assert( is_number(this->m_half_height)    || !"NAN encountered");
    assert( is_finite(this->m_half_height)    || !"INF encountered");
    assert( this->m_half_height >= 0 || !"Negative half height");
    assert( is_number(this->m_base_radius)    || !"NAN encountered");
    assert( is_finite(this->m_base_radius)    || !"INF encountered");
    assert( this->m_base_radius >= 0 || !"Negative base radius");

    T const vv = vx*vx + vy*vy + vz*vz;
    assert( is_number(vv)    || !"NAN encountered");
    assert( is_finite(vv)    || !"INF encountered");

    T px = 0;
    T py = 0;
    T pz = 0;

    // Test if we have a valid search direction. If not we simply return some boundary point of the cone
    if (vv <= 0 )
    {
      px = 0;
      py = 0;
      pz = -this->m_half_height;

      assert( is_number(px) || !"NAN encountered");
      assert( is_number(py) || !"NAN encountered");
      assert( is_number(pz) || !"NAN encountered");
      assert( is_finite(px) || !"INF encountered");
      assert( is_finite(py) || !"INF encountered");
      assert( is_finite(pz) || !"INF encountered");

      return V::make(px,py,pz);
    }

    /*
    // Now we know that we have a valid search direction. We proceed by
    // performing a case-by-case analysis. The theoretical back-ground
    // of the analysis is outlined below.

    // Since a Cone is radially symmetric around the z-axis we choose to
    // perform the analysis in the radial half-plane orthogonal to the
    // search direction.
    //
    // In this radial half plane the cone is simply a triangle. Looking like this
    //
    //                                  //
    //  z=+h --+        apex            //
    //         |\                       //
    //         | \                      //
    //  z=0 ---+  \                     //
    //         |   \                    //
    //         |    \                   //
    //  z=-h --+-----+---->             //
    //         |     |
    //        r=0   r=R
    //          base
    //
    // The height can be computed as
    //
    //   H = 2*h
    //
    // The angle, alpha at the apex, can be found from the
    // trigonometric relation
    //
    //   sin(alpha) = \frac{R}{H}
    //
    // If the appex is the support point then it must mean
    // that the search direction, s, is contained in the voronoi
    // region of the appex. By geometric relations we observe
    // that the voroni region is defined by the alpha-angle.
    //
    //                                _                             //
    //         |    /                 /|  search direction          //
    //         |  /                  /                              //
    //         |/  alpha            /                               //
    //  apex   +----------         / gamma                          //
    //         |\                 +-------                          //
    //         | \                                                  //
    //  z=0 ---+  \                                                 //
    //         |   \                                                //
    //         |    \                                               //
    //         +-----+---->                                         //
    //
    // So we need to test if angle, gamma, of the search direction is
    // larger than the apex angle. That is
    //
    //   gamma >= alpha   (*1)
    //
    // From geometry we find
    //
    // sine(gamma) =  \frac{s_z}{\norm{s}}
    //
    // Fortunately the sine is a monotome function of the
    // interval [0..pi/2] so (*1) can stated as
    //
    //    sine(gamma) >= sine(alpha)
    //
    // By substitution of trigonometric terms this can be rewritten as
    //
    //   s_z >= \norm{s} \frac{B}{H}
    //
    // Observe that we have included the voronoi region of the inclined
    // cone side into the voronoi region of the appex. Therefor the next
    // voroni-region to consider is the base corner.
    //
    // If the first test fails then it suffices to test if the search
    // direction has any component in the radial direction. If this is
    // the case then the base corner must be a support point.
    //
    // If the two previous test-case both fails it means we have
    // a downward search direction and any point on the base of the
    // cone is a support point.
    */

    T const sine_alpha = this->m_base_radius / (2*this->m_half_height);
    assert( is_number(sine_alpha)    || !"NAN encountered");
    assert( is_finite(sine_alpha)    || !"INF encountered");
    assert( sine_alpha >= 0 || !"sine(alpha) can not be negative");
    assert( sine_alpha <= 1  || !"sine(alpha) can not be larger than one");

    T const norm_v     = sqrt(vv);
    assert( is_number(norm_v)    || !"NAN encountered");
    assert( is_finite(norm_v)    || !"INF encountered");
    assert( norm_v >= 0 || !"Norm can not be negative");

    // Test if search direction is in Apex voronoi region
    if( vz >= norm_v*sine_alpha)
    {
      px = 0;
      py = 0;
      pz = this->m_half_height;

      assert( is_number(px) || !"NAN encountered");
      assert( is_number(py) || !"NAN encountered");
      assert( is_number(pz) || !"NAN encountered");
      assert( is_finite(px) || !"INF encountered");
      assert( is_finite(py) || !"INF encountered");
      assert( is_finite(pz) || !"INF encountered");

      return V::make(px,py,pz);
    }

    T const norm_sigma     = sqrt( vx*vx + vy*vy );
    assert( is_number(norm_sigma)    || !"NAN encountered");
    assert( is_finite(norm_sigma)    || !"INF encountered");
    assert( norm_sigma >= 0 || !"Cone::operator(): Norm can not be negative");

    // Test if search direction has any radial component
    if(norm_sigma > 0 )
    {
      px = this->m_base_radius*vx/norm_sigma;
      py = this->m_base_radius*vy/norm_sigma;
      pz = -(this->m_half_height);

      assert( is_number(px) || !"NAN encountered");
      assert( is_number(py) || !"NAN encountered");
      assert( is_number(pz) || !"NAN encountered");
      assert( is_finite(px) || !"INF encountered");
      assert( is_finite(py) || !"INF encountered");
      assert( is_finite(pz) || !"INF encountered");

      return V::make(px,py,pz);
    }

    // Search direction must be straight down
    px = 0;
    py = 0;
    pz = -this->m_half_height;

    assert( is_number(px) || !"NAN encountered");
    assert( is_number(py) || !"NAN encountered");
    assert( is_number(pz) || !"NAN encountered");
    assert( is_finite(px) || !"INF encountered");
    assert( is_finite(py) || !"INF encountered");
    assert( is_finite(pz) || !"INF encountered");

    return V::make(px,py,pz);
  }
  template<typename M>
  auto Cone<M>::get_support_point(EigenVector3<T> dir) const -> EigenVector3<T>
  {
      using std::sqrt;

      typedef typename M::real_type     T;
      typedef typename M::vector3_type  V;
      typedef typename M::value_traits VT;

      T const & vx = dir(0);
      T const & vy = dir(1);
      T const & vz = dir(2);

      assert( is_number(vx) || !"NAN encountered");
      assert( is_number(vy) || !"NAN encountered");
      assert( is_number(vz) || !"NAN encountered");
      assert( is_finite(vx) || !"INF encountered");
      assert( is_finite(vy) || !"INF encountered");
      assert( is_finite(vz) || !"INF encountered");

      assert( is_number(this->m_half_height)    || !"NAN encountered");
      assert( is_finite(this->m_half_height)    || !"INF encountered");
      assert( this->m_half_height >= 0 || !"Negative half height");
      assert( is_number(this->m_base_radius)    || !"NAN encountered");
      assert( is_finite(this->m_base_radius)    || !"INF encountered");
      assert( this->m_base_radius >= 0 || !"Negative base radius");

      T const vv = vx*vx + vy*vy + vz*vz;
      assert( is_number(vv)    || !"NAN encountered");
      assert( is_finite(vv)    || !"INF encountered");

      T px = 0;
      T py = 0;
      T pz = 0;

      // Test if we have a valid search direction. If not we simply return some boundary point of the cone
      if (vv <= 0 )
      {
          px = 0;
          py = 0;
          pz = -this->m_half_height;

          assert( is_number(px) || !"NAN encountered");
          assert( is_number(py) || !"NAN encountered");
          assert( is_number(pz) || !"NAN encountered");
          assert( is_finite(px) || !"INF encountered");
          assert( is_finite(py) || !"INF encountered");
          assert( is_finite(pz) || !"INF encountered");

          return {px,py,pz};
      }

      /*
    // Now we know that we have a valid search direction. We proceed by
    // performing a case-by-case analysis. The theoretical back-ground
    // of the analysis is outlined below.

    // Since a Cone is radially symmetric around the z-axis we choose to
    // perform the analysis in the radial half-plane orthogonal to the
    // search direction.
    //
    // In this radial half plane the cone is simply a triangle. Looking like this
    //
    //                                  //
    //  z=+h --+        apex            //
    //         |\                       //
    //         | \                      //
    //  z=0 ---+  \                     //
    //         |   \                    //
    //         |    \                   //
    //  z=-h --+-----+---->             //
    //         |     |
    //        r=0   r=R
    //          base
    //
    // The height can be computed as
    //
    //   H = 2*h
    //
    // The angle, alpha at the apex, can be found from the
    // trigonometric relation
    //
    //   sin(alpha) = \frac{R}{H}
    //
    // If the appex is the support point then it must mean
    // that the search direction, s, is contained in the voronoi
    // region of the appex. By geometric relations we observe
    // that the voroni region is defined by the alpha-angle.
    //
    //                                _                             //
    //         |    /                 /|  search direction          //
    //         |  /                  /                              //
    //         |/  alpha            /                               //
    //  apex   +----------         / gamma                          //
    //         |\                 +-------                          //
    //         | \                                                  //
    //  z=0 ---+  \                                                 //
    //         |   \                                                //
    //         |    \                                               //
    //         +-----+---->                                         //
    //
    // So we need to test if angle, gamma, of the search direction is
    // larger than the apex angle. That is
    //
    //   gamma >= alpha   (*1)
    //
    // From geometry we find
    //
    // sine(gamma) =  \frac{s_z}{\norm{s}}
    //
    // Fortunately the sine is a monotome function of the
    // interval [0..pi/2] so (*1) can stated as
    //
    //    sine(gamma) >= sine(alpha)
    //
    // By substitution of trigonometric terms this can be rewritten as
    //
    //   s_z >= \norm{s} \frac{B}{H}
    //
    // Observe that we have included the voronoi region of the inclined
    // cone side into the voronoi region of the appex. Therefor the next
    // voroni-region to consider is the base corner.
    //
    // If the first test fails then it suffices to test if the search
    // direction has any component in the radial direction. If this is
    // the case then the base corner must be a support point.
    //
    // If the two previous test-case both fails it means we have
    // a downward search direction and any point on the base of the
    // cone is a support point.
    */

      T const sine_alpha = this->m_base_radius / (2*this->m_half_height);
      assert( is_number(sine_alpha)    || !"NAN encountered");
      assert( is_finite(sine_alpha)    || !"INF encountered");
      assert( sine_alpha >= 0 || !"sine(alpha) can not be negative");
      assert( sine_alpha <= 1  || !"sine(alpha) can not be larger than one");

      T const norm_v     = sqrt(vv);
      assert( is_number(norm_v)    || !"NAN encountered");
      assert( is_finite(norm_v)    || !"INF encountered");
      assert( norm_v >= 0 || !"Norm can not be negative");

      // Test if search direction is in Apex voronoi region
      if( vz >= norm_v*sine_alpha)
      {
          px = 0;
          py = 0;
          pz = this->m_half_height;

          assert( is_number(px) || !"NAN encountered");
          assert( is_number(py) || !"NAN encountered");
          assert( is_number(pz) || !"NAN encountered");
          assert( is_finite(px) || !"INF encountered");
          assert( is_finite(py) || !"INF encountered");
          assert( is_finite(pz) || !"INF encountered");

          return {px,py,pz};
      }

      T const norm_sigma     = sqrt( vx*vx + vy*vy );
      assert( is_number(norm_sigma)    || !"NAN encountered");
      assert( is_finite(norm_sigma)    || !"INF encountered");
      assert( norm_sigma >= 0 || !"Cone::operator(): Norm can not be negative");

      // Test if search direction has any radial component
      if(norm_sigma > 0 )
      {
          px = this->m_base_radius*vx/norm_sigma;
          py = this->m_base_radius*vy/norm_sigma;
          pz = -(this->m_half_height);

          assert( is_number(px) || !"NAN encountered");
          assert( is_number(py) || !"NAN encountered");
          assert( is_number(pz) || !"NAN encountered");
          assert( is_finite(px) || !"INF encountered");
          assert( is_finite(py) || !"INF encountered");
          assert( is_finite(pz) || !"INF encountered");

          return {px,py,pz};
      }

      // Search direction must be straight down
      px = 0;
      py = 0;
      pz = -this->m_half_height;

      assert( is_number(px) || !"NAN encountered");
      assert( is_number(py) || !"NAN encountered");
      assert( is_number(pz) || !"NAN encountered");
      assert( is_finite(px) || !"INF encountered");
      assert( is_finite(py) || !"INF encountered");
      assert( is_finite(pz) || !"INF encountered");

      return {px,py,pz};
  }

  template<typename M>
  typename M::real_type Cone<M>::get_scale() const
  {
    using std::min;

    typedef typename M::real_type     T;
    typedef typename M::value_traits VT;

    assert( is_number(this->m_half_height)    || !"NAN encountered");
    assert( is_finite(this->m_half_height)    || !"INF encountered");
    assert( this->m_half_height >= 0 || !"Negative half height");
    assert( is_number(this->m_base_radius)    || !"NAN encountered");
    assert( is_finite(this->m_base_radius)    || !"INF encountered");
    assert( this->m_base_radius >= 0 || !"Negative base radius");

    T const d = 2 * ((this->m_base_radius > 0) ? this->m_base_radius : std::numeric_limits<T>::max());
    T const h = 2 * ((this->m_half_height > 0) ? this->m_half_height : std::numeric_limits<T>::max());

    return min(h, d);
  }

  using Mf = tiny::MathTypes<float>;
  using Md = tiny::MathTypes<double>;

  template class Cone<Mf>;
  template class Cone<Md>;

} // namespace convex
