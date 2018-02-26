#ifndef BIG_ABSOLUTE_CONVERGENCE_H
#define BIG_ABSOLUTE_CONVERGENCE_H

namespace big
{
  
  /**
   * Absolute Convergence Test Function.
   *
   * @param f            The function value at the current iterate.
   * @param tolerance    This argument holds the value used in the absolute
   *                     stopping criteria.
   *                     Setting the value to zero will make the test in-effective.
   * @return             If the absolute convergence test passes then the return
   *                     value is true otherwise it is false.
   */
  template < typename T >
  bool absolute_convergence(
                            T const & f
                            , T const & tolerance
                            );
  
} // namespace big

// BIG_ABSOLUTE_CONVERGENCE_H
#endif
