#ifndef BIG_COMPUTE_INDEX_REORDERING_H
#define BIG_COMPUTE_INDEX_REORDERING_H

#include <big_types.h>

namespace big
{
  
  /**
   * Compute Reordering of Constraints
   *
   * @param bitmask       A bitmask vector indicating set membership of the variables.
   *
   * @param old2new       Upon return this vector holds an index permutation of old indices
   *                      to new indices. The permutation is such that the first sub-block of
   *                      the new vector corresponds to active constraints. The second corresponds
   *                      to lower constraints and the last sub-block corresponds to the upper
   *                      constraints set.
   * @oaram new2old       Upon return this vector holds the reversible permuation of old2new.
   *
   */
  void  compute_index_reordering(
                                 ublas::vector<size_t> const & bitmask
                                 , ublas::vector<size_t>    & old2new
                                 , ublas::vector<size_t>    & new2old
                                 );
} // namespace big

// BIG_COMPUTE_INDEX_REORDERING_H
#endif
