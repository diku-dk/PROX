#ifndef HYPER_SCRIPTED_MOTION_H
#define HYPER_SCRIPTED_MOTION_H

#include <tiny.h>

#include <cmath>
#include <cassert>

namespace hyper
{

  template<typename MT>
  class Body;

  template<typename MT>
  class ScriptedMotion
  {
  public:

    virtual void compute( Body<MT> & body, typename MT::real_type const & time) = 0;

  };


} // namespace hyper

// HYPER_SCRIPTED_MOTION_H
#endif 
