#ifndef HYPER_ASSEMBLE_MATRIX_H
#define HYPER_ASSEMBLE_MATRIX_H

#include <hyper_math_policy.h>

#include <mesh_array.h>

#include <vector>
#include <cassert>

namespace hyper
{

  namespace details
  {
    /**
     * Add an unique Index to an Array.
     * It is assumed that the array has been pre-allocated and
     * is sufficiently large to hold what ever data needed.
     *
     * @param idx    The index value to be added.
     * @param C      An array that keeps all unique indices stored in
     *               incremental order. Upon return C is updated.
     * @param K      The number of unique indices currently stored in C. Upon
     *               return K will be updated.
     */
    inline void add_unique_index_to_array(
                                          unsigned int const & idx
                                          , std::vector<unsigned int> & C
                                          , unsigned int & K
                                          )
    {
      if (K==0u) //add to empty array
      {
        C[0] = idx;
        K = 1u;
      }
      else if ( C[K-1u] < idx )// add to back of array
      {
        C[K] = idx;
        K = K + 1u;
      }
      else if ( C[0] > idx )// add to front of array
      {
        for (unsigned int i = K; i>0;--i)
        {
          C[i] = C[i-1];
        }
        C[0] = idx;
        K = K + 1u;
      }else{  //add inside array
        unsigned int p = 0u;
        while (C[p] < idx)
        {
          p = p + 1u;
        }
        if ( C[p] != idx ) // only insert unique
        {
          for (unsigned int i = K; i>p;--i)
          {
            C[i] = C[i-1];
          }
          C[p] = idx;
          K = K + 1u;
        }
      }
    }

  }// details

  template<typename MT>
  inline void assemble_matrix(
                              mesh_array::T4Mesh  const & mesh
                              , mesh_array::Neighborhoods const & neighbors
                              , std::vector< typename MT::element_matrices_type > const & Ae
                              , typename MT::compressed_block3x3_type & A
                              )
    {
      typedef typename MT::matrix3x3_type            M;
//      typedef typename MT::vector3_type              V;
//      typedef typename MT::real_type                 T;
//      typedef typename MT::value_traits              VT;

//      typedef typename MT::vector_block3x1_type      vector_block3x1_type;
//      typedef typename MT::diagonal_block3x3_type    diagonal_block3x3_type;
//      typedef typename MT::compressed_block3x3_type  compressed_block3x3_type;
      typedef typename MT::element_matrices_type     element_matrices_type;

      unsigned int const UNASSIGNED = 0xFFFFFFFF;

      unsigned int const N          = neighbors.m_offset.size() - 1;

      assert(N == mesh.vertex_size() || !"assemble_diagonal(): mesh did not match neighborhood data");

      //--- Step 1: Compute how much space is needed  ----------------------------
    std::vector<unsigned int> row_sizes;
    row_sizes.resize( N );

    for(size_t r = 0u; r < N; ++r)
    {
      unsigned int const Kupper = (neighbors.m_offset[r+1] - neighbors.m_offset[r])*4u;

      std::vector<unsigned int> C;
      C.resize(Kupper,UNASSIGNED);

      unsigned int K = 0u;

      for (unsigned int entry = neighbors.m_offset[r]; entry < neighbors.m_offset[r+1];++entry)
      {
        unsigned int            const   idx         = neighbors.m_V2T[ entry ].second;
        mesh_array::Tetrahedron const & tetrahedron = mesh.tetrahedron(idx);

        unsigned int            const   global_i    = tetrahedron.i();
        unsigned int            const   global_j    = tetrahedron.j();
        unsigned int            const   global_k    = tetrahedron.k();
        unsigned int            const   global_m    = tetrahedron.m();

        details::add_unique_index_to_array(global_i, C, K);
        details::add_unique_index_to_array(global_j, C, K);
        details::add_unique_index_to_array(global_k, C, K);
        details::add_unique_index_to_array(global_m, C, K);
      }

      row_sizes[r] = K;
    }

    //--- Step 2: Allocate space  ----------------------------------------------
    unsigned int nnz = 0u;

    for(size_t i = 0; i < N; ++i)
    {
      nnz += row_sizes[i];
    }

    A.resize(N,N,nnz,false);

    //--- Step 3: Initialize indexing data and value data to zero  -------------
    for(size_t r = 0u; r < N; ++r)
    {
      unsigned int const Kupper = (neighbors.m_offset[r+1] - neighbors.m_offset[r])*4u;

      std::vector<unsigned int> C;

      C.resize(Kupper,UNASSIGNED);

      unsigned int K = 0u;

      for (unsigned int entry = neighbors.m_offset[r]; entry < neighbors.m_offset[r+1];++entry)
      {
        unsigned int            const   idx         = neighbors.m_V2T[ entry ].second;
        mesh_array::Tetrahedron const & tetrahedron = mesh.tetrahedron(idx);

        unsigned int            const   global_i    = tetrahedron.i();
        unsigned int            const   global_j    = tetrahedron.j();
        unsigned int            const   global_k    = tetrahedron.k();
        unsigned int            const   global_m    = tetrahedron.m();

        details::add_unique_index_to_array(global_i, C, K);
        details::add_unique_index_to_array(global_j, C, K);
        details::add_unique_index_to_array(global_k, C, K);
        details::add_unique_index_to_array(global_m, C, K);
      }

      for( unsigned int k = 0; k < K; ++k)
      {
        if (C[k] == UNASSIGNED)
          break;

        A(r,C[k]) = MT::convert( M::zero() );
      }

    }

    //--- Step 4: Fill in data  ------------------------------------------------
    for(size_t r = 0u; r < N; ++r)
    {
      for (unsigned int entry = neighbors.m_offset[r]; entry < neighbors.m_offset[r+1];++entry)
      {
        unsigned int            const   idx         = neighbors.m_V2T[ entry ].second;
        mesh_array::Tetrahedron const & tetrahedron = mesh.tetrahedron(idx);
        unsigned int            const   local_r     = tetrahedron.get_local_index(r);

        element_matrices_type const & E = Ae[idx];

        unsigned int const global_i = tetrahedron.i();
        unsigned int const global_j = tetrahedron.j();
        unsigned int const global_k = tetrahedron.k();
        unsigned int const global_m = tetrahedron.m();

        unsigned int const local_i  = 0u;
        unsigned int const local_j  = 1u;
        unsigned int const local_k  = 2u;
        unsigned int const local_m  = 3u;

        M const A_ri = MT::convert( A(r,global_i) ) + E(local_r, local_i);
        M const A_rj = MT::convert( A(r,global_j) ) + E(local_r, local_j);
        M const A_rk = MT::convert( A(r,global_k) ) + E(local_r, local_k);
        M const A_rm = MT::convert( A(r,global_m) ) + E(local_r, local_m);

        A(r, global_i) = MT::convert( A_ri );
        A(r, global_j) = MT::convert( A_rj );
        A(r, global_k) = MT::convert( A_rk );
        A(r, global_m) = MT::convert( A_rm );
      }
    }
  }
  
}// namespace hyper

// HYPER_ASSEMBLE_MATRIX_H
#endif
