#ifndef HYPER_ASSEMBLE_DIAGONAL_H
#define HYPER_ASSEMBLE_DIAGONAL_H

#include <hyper_math_policy.h>

#include <mesh_array.h>

#include <vector>
#include <cassert>

namespace hyper
{

  template<typename MT>
  inline void assemble_diagonal(
                                mesh_array::T4Mesh  const & mesh
                                , mesh_array::Neighborhoods const & neighbors
                                , std::vector< typename MT::element_matrices_type > const & Ae
                                , typename MT::diagonal_block3x3_type & D
                                )
  {
    typedef typename MT::matrix3x3_type            M;
    typedef typename MT::vector3_type              V;
    typedef typename MT::real_type                 T;
    typedef typename MT::value_traits              VT;

    typedef typename MT::vector_block3x1_type      vector_block3x1_type;
    typedef typename MT::diagonal_block3x3_type    diagonal_block3x3_type;
    typedef typename MT::compressed_block3x3_type  compressed_block3x3_type;
    typedef typename MT::element_matrices_type     element_matrices_type;

    unsigned int const N = neighbors.m_offset.size() - 1;

    assert(N == mesh.vertex_size() || !"assemble_diagonal(): mesh did not match neighborhood data");

    D.resize( N );

    for(size_t r = 0u; r < N; ++r)
    {
      M sum = M::zero();

      for (unsigned int entry = neighbors.m_offset[r]; entry < neighbors.m_offset[r+1];++entry)
      {
        unsigned int            const   idx         = neighbors.m_V2T[ entry ].second;
        mesh_array::Tetrahedron const & tetrahedron = mesh.tetrahedron(idx);
        unsigned int            const   local_r     = tetrahedron.get_local_index(r);

        element_matrices_type const & A = Ae[idx];

        assert(local_r != mesh_array::UNASSIGNED() || !"assemble_diagonal(): unassigned index, internal error");

        M const & A_rr = A(local_r,local_r);

        sum += A_rr;
      }

      D[r] = MT::convert( sum );
    }

  }
  
}// namespace hyper

// HYPER_ASSEMBLE_DIAGONAL_H
#endif
