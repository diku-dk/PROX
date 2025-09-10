#ifndef EIGENSPARSEHELPERS_H
#define EIGENSPARSEHELPERS_H

#include <cstddef>

#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>

namespace prox
{

template <typename T> using EigenDenseVec7 = Eigen::Matrix<T, 7, 1>;
template <typename T> using EigenDenseVec6 = Eigen::Matrix<T, 6, 1>;
template <typename T> using EigenDenseVec4 = Eigen::Matrix<T, 4, 1>;

template <typename T> using EigenSparseVec7 = std::vector<EigenDenseVec7<T>>;
template <typename T> using EigenSparseVec6 = std::vector<EigenDenseVec6<T>>;
template <typename T> using EigenSparseVec4 = std::vector<EigenDenseVec4<T>>;

template <typename T>
void computeB_Eigen(Eigen::SparseMatrix<T> const& J,
                    Eigen::VectorX<T> const& Wdth, Eigen::VectorX<T> const& u,
                    Eigen::VectorX<T> const& e, Eigen::VectorX<T> const& g,
                    Eigen::VectorX<T>& b)
{
    b.resize(g.size());

    Eigen::VectorX<T> v_minus = J * u;
    Eigen::VectorX<T> v_plus = v_minus;
    v_plus += J * Wdth;

    // Perform operations in the exact order of the original
    b = e.cwiseProduct(v_minus); // Element-wise product
    b += g; // Add stabilization term
    b += v_plus; // Add v_plus
}

/*template <typename T>
void computeB_Eigen(Eigen::SparseMatrix<T> const& J,
                    Eigen::VectorX<T> const& Wdth, Eigen::VectorX<T> const& u,
                    Eigen::VectorX<T> const& e, Eigen::VectorX<T> const& g,
                    Eigen::VectorX<T>& b)
{
    b.resize(g.size());

    Eigen::VectorX<T> v_minus(g.size());
    Eigen::VectorX<T> v_plus(g.size());

    // v_minus = J*u
    v_minus = J * u;

    // v_plus = v_minus + J*Wdth
    v_plus = v_minus + J * Wdth;

    // b = (e .* v_minus) + g + v_plus (element-wise product and addition)
    b = e.array() * v_minus.array() + g.array() + v_plus.array();
}*/

template <typename T>
void computeWJT_Eigen(Eigen::SparseMatrix<T> const& W,
                      Eigen::SparseMatrix<T> const& J,
                      Eigen::SparseMatrix<T>& WJT)
{
    // Compute JT = J.transpose() explicitly
    Eigen::SparseMatrix<T> JT = J.transpose();

    // Compute WJT = W * JT
    WJT = W * JT;
}

} // namespace prox

#endif //EIGENSPARSEHELPERS_H
