#include <big_return_codes.h>

namespace big
{
  
  std::string get_error_message(size_t const & error_code)
  {
    std::string msg;
    switch(error_code)
    {
      case OK:                               msg = "no error"; break;
      case NON_DESCENT_DIRECTION:            msg = "Non descent direction was encountered"; break;
      case BACKTRACKING_FAILED:              msg = "Back-tracking failure during line-search"; break;
      case STAGNATION:                       msg = "Stagnation test passed"; break;
      case RELATIVE_CONVERGENCE:             msg = "Relative convergence test passed"; break;
      case ABSOLUTE_CONVERGENCE:             msg = "Absolute convergence test passed"; break;
      case ITERATING:                        msg = "Halted while iterating or did not converge with maximum number of iterations"; break;
      case DESCEND_DIRECTION_IN_NORMAL_CONE: msg = "Descend Direction is in Normal Cone"; break;
      case LOCAL_MINIMA:                     msg = "Local minimia"; break;
      case MAX_LIMIT:                        msg = "MAximum limit reached"; break;
      default:                               msg = "unrecognised error"; break;
    };
    return msg;
  }
  
} // namespace big

