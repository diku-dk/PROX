#ifndef PROX_VELOCITY_UPDATE_H
#define PROX_VELOCITY_UPDATE_H

#include "eigenhelperfunctions.h"

namespace prox
{


template <typename T>
inline void
velocity_update_eigen(const Eigen::VectorX<T>& u, const Eigen::VectorX<T>& Wdth,
                      const Eigen::VectorX<T>& fc, Eigen::VectorX<T>& unew)
{
    unew = u + Wdth + fc;
}

template <typename T>
inline void velocity_update_eigen(const Eigen::VectorX<T>& u,
                                  const Eigen::VectorX<T>& Wdth,
                                  Eigen::VectorX<T>& unew)
{
    unew = u + Wdth;
}
} // namespace prox

// PROX_VELOCITY_UPDATE_H
#endif
