#ifndef PROX_NORMAL_SUB_SOLVER_H
#define PROX_NORMAL_SUB_SOLVER_H

#include "prox_enums.h"
#include <algorithm>
#include <stdexcept>

namespace prox
{

template <typename T>
void normalSolver(normal_sub_solver_type type, const T& z_n, T& lambda_n)
{
    switch (type)
    {
    case nonnegative:     lambda_n = std::max<T>(0, z_n); return;
    case normal_origin:   lambda_n = 0; return;
    case normal_infinity: return;
    default:              throw std::logic_error("Code is broken.");
    }
}

} //namespace prox

// PROX_NORMAL_SUB_SOLVER_H
#endif
