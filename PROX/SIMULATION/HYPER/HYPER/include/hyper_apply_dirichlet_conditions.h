#ifndef HYPER_APPLY_DIRICHLET_CONDITIONS_H
#define HYPER_APPLY_DIRICHLET_CONDITIONS_H

#include <hyper_dirichlet_info.h>
#include <hyper_math_policy.h>

#include <vector>

namespace hyper
{

  /**
   * Apply Dirichlet Conditions to linear system A x = b where A is general
   * sparse matrix.
   */
  template<typename MT>
  inline void apply_dirichlet_conditions(
                                           std::vector<DirichletInfo<MT> > const & dirichlet_conditions
                                         , typename MT::compressed_block3x3_type & A
                                         , typename MT::vector_block3x1_type & b
                                         )
  {
    typedef typename MT::matrix3x3_type            M;
    typedef typename MT::vector3_type              V;
//    typedef typename MT::real_type                 T;
//    typedef typename MT::value_traits              VT;

    typedef typename MT::vector_block3x1_type      vector_block3x1_type;
//    typedef typename MT::diagonal_block3x3_type    diagonal_block3x3_type;
    typedef typename MT::compressed_block3x3_type  compressed_block3x3_type;
//    typedef typename MT::element_matrices_type     element_matrices_type;


    typedef typename compressed_block3x3_type::row_iterator  row_iterator;

    //--- Adjust right hand side -----------------------------------------------
    {
      vector_block3x1_type dirichlet_values;

      dirichlet_values.resize( b.size() );
      dirichlet_values.clear_data();

      for(unsigned int i = 0u; i < dirichlet_conditions.size(); ++i)
      {
        unsigned int const & idx   = dirichlet_conditions[i].idx();
        V            const & value = dirichlet_conditions[i].value();

        dirichlet_values(idx) = MT::convert( value );
      }

      vector_block3x1_type tmp;

      sparse::prod(A, dirichlet_values, tmp, true);
      sparse::sub(tmp ,b);
    }
    //--- Make bit-mask for fast testing ---------------------------------------
    std::vector<bool> is_dirichlet(b.size(), false);

    for(unsigned int i = 0u; i < dirichlet_conditions.size(); ++i)
    {
      unsigned int const & idx = dirichlet_conditions[i].idx();
      is_dirichlet[idx]        = true;
    }

    //--- Adjust all columns of A ----------------------------------------------
    for(size_t r = 0u; r < A.nrows(); ++r)
    {
      row_iterator row_iter = A.row_begin(r);
      row_iterator row_end  = A.row_end(r);

      for(;row_iter != row_end; ++row_iter)
      {
        unsigned int const c = col(row_iter);

        *row_iter = is_dirichlet[c] ? MT::convert( M::zero() ) : *row_iter;
      }
    }

    //--- Adjust all Dirichlet rows of A ---------------------------------------
    for(unsigned int i = 0u; i < dirichlet_conditions.size(); ++i)
    {
      unsigned int const & idx   = dirichlet_conditions[i].idx();
      V            const & value = dirichlet_conditions[i].value();

      row_iterator row_iter = A.row_begin(idx);
      row_iterator row_end  = A.row_end(idx);

      for(;row_iter != row_end; ++row_iter)
      {
        *row_iter = MT::convert( M::zero() );
      }

      A(idx,idx) = MT::convert( M::identity() );
      b(idx)     = MT::convert( value         );
    }
  }

  /**
   * Apply Dirichlet Conditions to linear system A x = b where A is a
   * diagonal matrix.
   */
  template<typename MT>
  inline void apply_dirichlet_conditions(
                                           std::vector<DirichletInfo<MT> > const & dirichlet_conditions
                                         , typename MT::diagonal_block3x3_type & A
                                         , typename MT::vector_block3x1_type & b
                                         )
  {
    typedef typename MT::matrix3x3_type            M;
    typedef typename MT::vector3_type              V;
//    typedef typename MT::real_type                 T;
//    typedef typename MT::value_traits              VT;

//    typedef typename MT::vector_block3x1_type      vector_block3x1_type;
//    typedef typename MT::diagonal_block3x3_type    diagonal_block3x3_type;
//    typedef typename MT::compressed_block3x3_type  compressed_block3x3_type;
//    typedef typename MT::element_matrices_type     element_matrices_type;

    for(unsigned int i = 0u; i < dirichlet_conditions.size(); ++i)
    {
      unsigned int const & idx   = dirichlet_conditions[i].idx();
      V            const & value = dirichlet_conditions[i].value();

      b(idx)     = MT::convert( value );
      A(idx,idx) = MT::convert( M::identity() );
    }
  }
  
}// namespace hyper

// HYPER_APPLY_DIRICHLET_CONDITIONS_H
#endif
