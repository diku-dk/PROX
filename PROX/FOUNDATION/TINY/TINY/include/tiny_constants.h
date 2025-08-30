#ifndef TINY_CONSTANTS_H
#define TINY_CONSTANTS_H

#include <boost/numeric/conversion/bounds.hpp>  // for  boost::numeric::bounds<T>::lowest() and  boost::numeric::bounds<T>::highest()
#include <boost/cast.hpp>                       // for boost::numeric_cast

#include <cmath>

#ifdef _MSC_VER
#  ifndef M_PI
#    define M_PI 3.14159265358979323846 /* pi */
#  endif
#  ifndef M_PI_2
#    define M_PI_2 1.57079632679489661923 /* pi/2 */
#  endif
#  ifndef M_PI_4
#    define M_PI_4 0.78539816339744830962  /* pi/4 */
#  endif
#endif

namespace tiny
{

    namespace detail
    {

      template<typename T>
      inline T pi() { return boost::numeric_cast<T>(M_PI); }

      template<typename T>
      inline T pi_half() {  return boost::numeric_cast<T>(M_PI_2); }

      template<typename T>
      inline T pi_quarter() {  return boost::numeric_cast<T>(M_PI_4); }

      template <typename T>
      inline T highest() { return boost::numeric::bounds<T>::highest(); }

      template <typename T>
      inline T lowest() { return boost::numeric::bounds<T>::lowest(); }

      template <typename T>
      inline T infinity(){ return highest<T>(); }

      template<typename T>
      inline T radians_per_degree() {  return boost::numeric_cast<T>(0.017453292519943295769236907684886); }

      template<typename T>
      inline T degrees_per_radian() {  return boost::numeric_cast<T>(57.295779513082320876798154814105); }

      template<typename R, typename T>
      inline R numeric_cast(T const & value) {  return boost::numeric_cast<R>(value); }

    } // namespace detail

}  // namespace tiny

// TINY_CONSTANTS_H
#endif
