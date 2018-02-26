#ifndef BIG_RETURN_CODES_H
#define BIG_RETURN_CODES_H

#include <string>

namespace big
{
  
  /**
   * @file
   * This file contains declarations of some readable constants that is
   * used for setting proper error code/status before exiting any given
   * numerical method in the optimization library. The intention is that
   * this ensures that all numerical error codes have the same global meaning
   * regardless of the methods/functions that is used.
   *
   * @author Kenny Erleben
   */
  
  /**
   * No error occured.
   */
  size_t const OK = 0;
  
  /**
   * A non descent direction was encountered. This may occur
   * during a line-search method. In which case one expects
   * to perform a line-search along a descent direction.
   */
  size_t const NON_DESCENT_DIRECTION = 1;
  
  /**
   * A back-tracking operation failure was detected during some
   * line-search method. This often indicates that the step-length
   * have become too small and no decrease in (merit) function
   * value was obtained.
   */
  size_t const BACKTRACKING_FAILED = 2;
  
  /**
   * Stagnation occured.
   * Stagnation means that the maximum difference between the
   * components of a new iterate and the old iterate dropped
   * below some small threshold value. Basically it means
   * that the two iterates are nearly the same and no progress
   * have been made by the numerical method used.
   */
  size_t const STAGNATION = 3;
  
  /**
   * Relative Convergence Test Succeded.
   * This means that the relative improvement in function value has dropped below
   * a given threshold value. Given the current iterate \f$x^{n+1}\f$ and the
   * previous iterate\f$x^n\f$ then the test is
   *
   * \f[ \frac{| f(x^{n+1}) - f(x^n) |}{|f(x^n)|} < \varepsilon,  \f]
   *
   *where \f$\varepsilon \eq 0\f$ is a user specified test-threshold.
   */
  size_t const RELATIVE_CONVERGENCE = 4;
  
  /**
   * Absolute Convergence Test Succeded.
   * This means that the absolute function value has dropped below
   * a given threshold value. Given the current iterate \f$x^{n}\f$ then the test is
   *
   * \f[ f(x^n) < \varepsilon,  \f]
   *
   *where \f$\varepsilon \eq 0\f$ is a user specified test-threshold. In case
   * of minimization the test could also be for a stationary point. The test
   * would then be
   *
   * \f[ | \nabla f(x^n) | < \varepsilon,  \f]
   *
   */
  size_t const ABSOLUTE_CONVERGENCE = 5;
  
  /**
   * Halted while iterating.
   * This means that somehow for whatever unknown reason the method/function has
   * halted while iterating. In some cases this indicates some internal error. In
   * some cases this may occur if a method/function is being feed a bad (ill-posed
   * or inconsistent) problem to solve.
   *
   * However in most cases it simply means that the method did not converge within
   * some given maximum number of iterations.
   */
  size_t const ITERATING  = 6;
  
  /**
   * Descend Direction is in Normal Cone.
   * This means that the current iterate must be placed on the boundary of
   * the feasible region and that the descend direction is point away from
   * the feasible region in the normal direction. Thus one can not even
   * slide along the boundary of the feasible region in order to minimize
   * ones function value.
   *
   * In many cases this means that the current iterate is as good as its
   * get (when considering local information only).
   */
  size_t const DESCEND_DIRECTION_IN_NORMAL_CONE  = 7;

  
  /**
   *
   */
  size_t const LOCAL_MINIMA  = 8;
  
  /**
   *
   */
  size_t const MAX_LIMIT  = 9;
  
  /**
   * Get Error Message.
   * This function decodes an given error code value
   * into a user friendly and human readable text string. The
   * text string may be convenient for displaying error
   * messages in a log or on screen for an end-user.
   *
   * @param error_code   The value of an error code.
   *
   * @return             A textual and human readable error message string.
   */
  std::string get_error_message(size_t const & error_code);
  
  /**
   * In Non-Active Constraint Set Bitmask.
   */
  size_t const IN_NON_ACTIVE = 0u;
  
  /**
   * In Active Constraint Set Bitmask.
   */
  size_t const IN_ACTIVE = 1u;
  
  
} // namespace big

// BIG_RETURN_CODES_H
#endif
