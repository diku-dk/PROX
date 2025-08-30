#ifndef PROX_NONNEGATIVE_H
#define PROX_NONNEGATIVE_H

#include <algorithm>

namespace prox::detail
{

template <typename T>
inline static void nonnegative(const T& z_n, T& lambda_n)
{
    lambda_n = std::max<T>(0, z_n);
}

} // namespace prox::detail

// PROX_NONNEGATIVE_H
#endif
