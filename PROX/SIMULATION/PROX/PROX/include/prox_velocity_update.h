#ifndef PROX_VELOCITY_UPDATE_H
#define PROX_VELOCITY_UPDATE_H

#include "prox_math.h"

namespace prox
{

template <typename T>
inline void velocity_update(const NCVec6<T>& u, const NCVec6<T>& Wdth,
                            const NCVec6<T>& fc, NCVec6<T>& unew)
{
    unew.resize(u.nrows());
    sparse::add(u, Wdth, fc, unew);
}

template <typename T>
inline void velocity_update(const NCVec6<T>& u, const NCVec6<T>& Wdth,
                            NCVec6<T>& unew)
{
    unew.resize(u.nrows());
    sparse::add(u, Wdth, unew);
}

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
