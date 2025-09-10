#ifndef PROX_MATH_H
#define PROX_MATH_H

#include <cstddef>

#include <sparse.h>
#include <tiny_math_types.h>

#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>

namespace prox
{

template <size_t N, size_t M, typename T>
using MatrixBlock = sparse::Block<N, M, T>;

template <size_t N, typename T>
using DiagonalMatrix = sparse::DiagonalMatrix<MatrixBlock<N, N, T>>;

template <size_t N, size_t M, typename T>
using CRMatrix = sparse::CompressedRowMatrix<MatrixBlock<N, M, T>>;

template <size_t N, typename T>
using NeedlesslyComplicatedVector = sparse::Vector<MatrixBlock<N, 1, T>>;

template <typename T> using NCVec4 = NeedlesslyComplicatedVector<4, T>;
template <typename T> using NCVec6 = NeedlesslyComplicatedVector<6, T>;
template <typename T> using NCVec7 = NeedlesslyComplicatedVector<7, T>;

template <typename T> using EigenDenseVec7 = Eigen::Matrix<T, 7, 1>;
template <typename T> using EigenDenseVec6 = Eigen::Matrix<T, 6, 1>;
template <typename T> using EigenDenseVec4 = Eigen::Matrix<T, 4, 1>;

template <typename T> using EigenSparseVec7 = std::vector<EigenDenseVec7<T>>;
template <typename T> using EigenSparseVec6 = std::vector<EigenDenseVec6<T>>;
template <typename T> using EigenSparseVec4 = std::vector<EigenDenseVec4<T>>;

template <typename T> NCVec7<T> convertEigenToSparse(EigenSparseVec7<T>& input)
{
    NCVec7<T> h;
    h.resize(input.size());
    for (size_t i = 0; i < input.size(); ++i)
    {
        auto& b = h(i);

        b(0) = input[i](0);
        b(1) = input[i](1);
        b(2) = input[i](2);
        b(3) = input[i](3);
        b(4) = input[i](4);
        b(5) = input[i](5);
        b(6) = input[i](6);
    }
}

template <typename T>
sparse::Vector<MatrixBlock<6, 1, T>>
convertEigenToSparse(EigenSparseVec6<T>& input)
{
    NCVec6<T> h;
    h.resize(input.size());
    for (size_t i = 0; i < input.size(); ++i)
    {
        MatrixBlock<6, 1, T>& b = h(i);

        b(0) = input[i](0);
        b(1) = input[i](1);
        b(2) = input[i](2);
        b(3) = input[i](3);
        b(4) = input[i](4);
        b(5) = input[i](5);
    }
    return h;
}

template <typename T> NCVec4<T> convertEigenToSparse(EigenSparseVec4<T>& input)
{
    NCVec4<T> h;
    h.resize(input.size());
    for (size_t i = 0; i < input.size(); ++i)
    {
        auto& b = h(i);

        b(0) = input[i](0);
        b(1) = input[i](1);
        b(2) = input[i](2);
        b(3) = input[i](3);
    }
}

//compute z = x - R(J W J^T x  + b)  = x - R( A x  + b)
template <typename T>
void computeZ(const NCVec4<T>& x, const DiagonalMatrix<4, T>& R,
              const CRMatrix<4, 6, T>& J, const CRMatrix<6, 4, T>& WJT,
              const NCVec4<T>& b, NCVec4<T>& z)
{
    NCVec6<T> t1(WJT.nrows());
    NCVec4<T> t2(x.nrows()), t3(x.nrows());
    z.resize(x.size());
    sparse::prod(WJT, x, t1);       // t1 = WJT*x
    sparse::prod(
        J, t1, t2,
        false); // t2 = J*t1    // 2009-08-04 Kenny: We should probably extend SPARSE with axpy_prod for compressed row matrices
    sparse::add(b, t2, t2); // t2 += b
    sparse::prod(R, t2, t3); // t3 = R*t2
    sparse::sub(x, t3, z); // z  = x-t3
}

template <typename T> T computeInfNorm(const NCVec4<T>& x)
{
    // 2009-08-04 Kenny: Argh, this is ugly ugly ugly code, if block4_type were a tiny::Vector<4,T> type then this could be made much more pretty
    // 2025-09-06: Yup, agreement all around.

    assert(x.size() > 0 || !"compute_norm_inf() zero sized vector");
    T max = 0;
    for (const auto& block : x)
        for (size_t i = 0; i < 4; ++i)
        {
            if (auto cand = std::abs(block(i)); cand > max) max = cand;
        }
    return max;
}

template <typename T>
void computeZk(const MatrixBlock<4, 1, T>& x_k, const NCVec6<T>& w,
               const MatrixBlock<4, 4, T>& R_k, const CRMatrix<4, 6, T>& J,
               const MatrixBlock<4, 1, T>& b_k, MatrixBlock<4, 1, T>& z_k,
               const size_t& k)
{
    // z_k = x_k - R_kk ( J w + b_k )
    z_k.clear_data();
    MatrixBlock<4, 1, T> temp(0);
    sparse::row_prod(J, w, z_k, k); // z_k += J_k w
    sparse::add(b_k, z_k); // z_k += b_k
    sparse::prod(R_k, z_k, temp); // temp += R_kk z_k
    sparse::sub(x_k, temp, z_k); // z_k = x_k - temp
}

// 2009-08-13 Kenny: We should add unit-test for this method, it is a core method and we need to make sure it works correctly.
// compute b = (I+E)J u + J W (dt h)
template <typename T>
void computeB(CRMatrix<4, 6, T> const& J, NCVec6<T> const& Wdth,
              NCVec6<T> const& u, NCVec4<T> const& e, NCVec4<T> const& g,
              NCVec4<T>& b)
{
    b.resize(g.size());

    NCVec4<T> v_minus(g.size());
    NCVec4<T> v_plus(g.size());

    sparse::prod(J, u, v_minus); // v_minus = J*u
    v_plus = v_minus;
    sparse::prod(J, Wdth, v_plus, false); // v_plus  += J*Wdth
    sparse::element_prod(e, v_minus, b); // b = (e#v_minus)
    sparse::add(g, b); // b += g : correction/stabilization term
    sparse::add(v_plus, b); // b += v_plus
}


// 2009-08-13 Kenny: Optimization replace with diagonal_mass_type, maybe wait to optimize until all it working
template <typename T>
void computeWJT(DiagonalMatrix<6, T> const& W, CRMatrix<4, 6, T> const& J,
                CRMatrix<6, 4, T>& WJT)
{
    CRMatrix<6, 4, T> JT;

    sparse::transpose(J, JT);
    WJT.resize(W.nrows(), JT.ncols(), JT.size());
    sparse::prod(W, JT, WJT);
}

} // namespace prox

#endif // PROX_MATH_H
