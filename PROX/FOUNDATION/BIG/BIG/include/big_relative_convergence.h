#ifndef BIG_RELATIVE_CONVERGENCE_H
#define BIG_RELATIVE_CONVERGENCE_H

namespace big
{
  
  /**
   * Relative Convergence Test Function.
   *
   * @param f_old        The function value at the previous iterate.
   * @param f            The function value at the current iterate.
   * @param tolerance    This argument holds the value used in the relative stopping criteria.
   *                     Setting the value to zero will make the test in-effective.
   * @return             If the relative convergence test passes then the return value is true otherwise it is false.
   */
  template < typename T >
  bool relative_convergence(
                            T const & f_old
                            , T const & f
                            , T const & tolerance
                            );

} // namespace big

// BIG_RELATIVE_CONVERGENCE_H
#endif
