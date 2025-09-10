#ifndef EIGENHELPERALL_H
#define EIGENHELPERALL_H

#include <eigen3/Eigen/Dense>
#include <numbers>
#include <random>

template <typename Number> using EigenVector3 = Eigen::Vector<Number, 3>;

template <typename Number> using EigenVector4 = Eigen::Vector<Number, 4>;

template <typename Number> using EigenMatrix3 = Eigen::Matrix<Number, 3, 3>;

template <typename Number> using EigenQuaternion = Eigen::Quaternion<Number>;

template <typename Number> using EigenQuaternion = Eigen::Quaternion<Number>;

template <typename Number> using EigenMatrix3 = Eigen::Matrix<Number, 3, 3>;

template <typename Vector> inline auto dot(Vector a, Vector b)
{
    return a.dot(b);
}

template <typename Vector> inline auto abs(Vector a) { return a.cwiseAbs(); }

template <typename Vector> inline auto cross(Vector a, Vector b)
{
    return a.cross(b);
}

inline EigenQuaternion<float> prodEigen(const EigenQuaternion<float>& a,
                                        const EigenQuaternion<float>& b)
{
    return a * b;
}

inline EigenQuaternion<double> prodEigen(const EigenQuaternion<double>& a,
                                         const EigenQuaternion<double>& b)
{
    return a * b;
}

template <typename Vector> inline auto norm(Vector a) { return a.norm(); }

template <typename Vector> inline auto unit(Vector a) { return a.normalized(); }

template <typename Number> class CoordSysEigen
{
public:
protected:
    EigenVector3<Number> m_T;      ///< The Position.
    EigenQuaternion<Number> m_Q;      ///< The orientation in Quaternion form.

public:
    EigenVector3<Number>& T() { return m_T; }

    const EigenVector3<Number>& T() const { return m_T; }

    EigenQuaternion<Number>& Q() { return m_Q; }

    const EigenQuaternion<Number>& Q() const { return m_Q; }

public:
    CoordSysEigen()
        : m_T({0, 0, 0})
        , m_Q({1, 0, 0, 0})
    {
    }

    CoordSysEigen(const CoordSysEigen& X)
        : m_T(X.m_T)
        , m_Q(X.m_Q)
    {
    }

    CoordSysEigen(EigenVector3<Number> const& T_val,
                  const EigenQuaternion<Number>& Q_val)
    {
        m_T = T_val;
        m_Q = (Q_val).normalized();
    }

    CoordSysEigen(EigenVector3<Number> const& T_val,
                  const EigenMatrix3<Number>& R_val)
    {
        m_T = T_val;
        EigenQuaternion<Number> Quat(R_val);
        m_Q = Quat;
    }

    CoordSysEigen& operator=(CoordSysEigen const& C)
    {
        m_T = C.m_T;
        m_Q = C.m_Q;
        return *this;
    }

    bool operator==(CoordSysEigen const& C) const
    {
        return m_T == C.m_T && m_Q == C.m_Q;
    }

public:
    static CoordSysEigen identity()
    {
        CoordSysEigen X;
        X.T() = {0, 0, 0};
        X.Q() = EigenQuaternion<Number>::Identity();
        return X;
    }

    static CoordSysEigen make(const EigenVector3<Number>& T,
                              const EigenQuaternion<Number>& Q)
    {
        CoordSysEigen X;
        X.T() = T;
        X.Q() = Q;
        return X;
    }
};

template <typename Number>
inline EigenVector3<Number> rotate(const EigenQuaternion<Number>& Quat,
                                   const EigenVector3<Number>& rVec)
{
    //return prod(  prod(q , r)  , conj(q) ).imag();
    return Quat * rVec;
}

template <typename Number>
static EigenQuaternion<Number> Rotateu(const Number rad,
                                       const EigenVector3<Number>& axis)
{
    using std::cos;
    using std::sin;
    const Number theta = rad / 2;
    const Number ctheta = (cos(theta));
    const Number stheta = (sin(theta));
    EigenVector3<Number> vec = (axis).normalized() * stheta;
    return EigenQuaternion<Number>(ctheta, vec.x(), vec.y(), vec.z());
}

template <typename Number>
static EigenQuaternion<Number> Rotatex(const Number rad)
{
    using std::cos;
    using std::sin;
    const Number theta = rad / 2;
    const Number ctheta = (cos(theta));
    const Number stheta = (sin(theta));
    return EigenQuaternion<Number>(ctheta, stheta, 0, 0);
}

template <typename Number>
static EigenQuaternion<Number> Rotatey(const Number rad)
{
    using std::cos;
    using std::sin;
    const Number theta = rad / 2;
    const Number ctheta = (cos(theta));
    const Number stheta = (sin(theta));
    return EigenQuaternion<Number>(ctheta, 0, stheta, 0);
}

template <typename Number>
static EigenQuaternion<Number> Rotatez(const Number rad)
{
    using std::cos;
    using std::sin;
    const Number theta = rad / 2;
    const Number ctheta = (cos(theta));
    const Number stheta = (sin(theta));
    return EigenQuaternion<Number>(ctheta, 0, 0, stheta);
}

/*template<typename Number>
inline EigenVector3<Number> rotate(const EigenQuaternion<Number>& Quat, const Eigen::Matrix<Number, 3,1>::NeagtiveReturnType& rVec)
{
    //return prod(  prod(q , r)  , conj(q) ).imag();
    return Quat*rVec;
}*/

template <typename T>
inline void getAxisAngle(const EigenQuaternion<T>& Q, EigenVector3<T>& axis,
                         T& theta)
{
    using std::atan2;

    //
    // By definition a unit quaternion Q can be written as
    //
    //    Q = [s,v] = [cos(theta/2), n sin(theta/2)]
    //
    // where n is a unit vector. This is the same as a rotation of
    // theta radian around the axis n.
    //
    //
    // Rotations are difficult to work with for several reasons.
    //
    // Firstly both Q and -Q represent the same rotation. This is
    // easily proven, rotate a arbitary vector r by Q then we have
    //
    //   r^\prime = Q r Q^*
    //
    // Now rotate the same vector by -Q
    //
    //   r^\prime = (-Q) r (-Q)^* = Q r Q^*
    //
    // because -Q = [-s,-v] and (-Q)^* = [-s , v] = - [s,-v]^* = - Q^*.
    //
    // Thus the quaternion representation of a single rotation is not unique.
    //
    // Secondly the rotation it self is not well-posed. A rotation of theta
    // radians around the unit axis n could equally well be done as a rotation
    // of -theta radians around the negative unit axis n.
    //
    // This is seen by straightforward substitution
    //
    //  [ cos(-theta/2), sin(-theta/2) (-n) ] = [ cos(theta/2), sin(theta/2) n ]
    //
    // Thus we get the same quaternion regardless of whether we
    // use (+theta,+n) or (-theta,-n).
    //
    //
    // From the Quaternion we see that
    //
    //   \frac{v}{\norm{v}}  = \frac{ sin(theta/2) n }{| sin(theta/2) | } = sign(sin(theta/2)) n
    //
    // Thus we can easily get the rotation axis. However, we can not immediately
    // determine the positive rotation axis direction. The problem boils down to the
    // fact that we can not see the sign of the sinus-factor.
    //
    // Let us proceed by setting
    //
    //   x =    cos(theta/2)   =  s
    //   y =  | sin(theta/2) | =  \norm{v}
    //
    // Then we basically have two possibilities for finding theta
    //
    //  theta_1 = 2 atan2( y, x)        equivalent to      sign(sin(theta/2)) = 1
    //
    // or
    //
    //  theta_2 = 2 atan2( -y, x)       equivalent to      sign(sin(theta/2)) = -1
    //
    // If theta_1 is the solution we have
    //
    //  n = \frac{v}{\norm{v}}
    //
    // If theta_2 is the solution we must have
    //
    //  n = - \frac{v}{\norm{v}}
    //
    // Observe that we always have theta_2 = 2 pi - theta_1. Therefore theta_1 < theta_2.
    //
    // Let us imagine that we always choose $theta_1$ as the solution then
    // the correspoding quaternion for that solution would be
    //
    //         Q_1 = [cos(theta_1/2),  sin(theta_1/2)   \frac{v}{\norm{v}}]
    //             = [s ,  \norm{v}   \frac{v}{\norm{v}}]
    //             = Q
    //
    // Now if we choose theta_2 as the solution we would have
    //
    //         Q_2 = [cos(theta_2/2),  sin(theta_2/2)   -\frac{v}{\norm{v}}]
    //             = [s ,  -\norm{v}   -\frac{v}{\norm{v}}]
    //             = [s ,  \norm{v}   \frac{v}{\norm{v}}]
    //             = Q
    //
    // Thus we observe that regardless of which solution we pick we always have Q = Q_1 = Q_2.
    //
    // At this point one may be confused. However, it should be clear that theta_2 is equivalent
    // to the theta_1 rotation. The difference is simply that theta_2 corresponds to flipping the
    // rotation axis of the theta_1 case.
    //
    const T ct2 = Q.w(); //---   cos(theta/2)
    const T st2
        = (EigenVector3<T>(Q.x(), Q.y(), Q.z())).norm(); //---  |sin(theta/2)|

    theta = 2 * atan2(st2, ct2);

    assert(st2 >= 0
           || !"get_axis_angle(): |sin(theta/2)| must be non-negative");
    assert(theta >= 0 || !"get_axis_angle(): theta must be non-negative");

    axis = st2 > 0 ? EigenVector3<T>(Q.x(), Q.y(), Q.z()) / st2
                   : EigenVector3<T>(0, 0, 0);
}

template <typename T>
inline CoordSysEigen<T>
model_update(const EigenVector3<T>& TA, const EigenQuaternion<T>& QA,
             const EigenVector3<T>& TB, const EigenQuaternion<T>& QB)
{
    //---
    //---  p' = RA p + TA         (*1)  from A->WCS
    //---
    //---  p = RB^T (p' - TB)     (*2)  from WCS-B
    //---
    //--- Insert (*1) into (*2)  A -> B
    //---
    //---   p = RB^T ( RA p + TA - TB)
    //---     =  RB^T  RA p + RB^T (TA - TB)
    //--- So
    //---   R = RB^T  RA
    //---   T = RB^T (TA - TB)
    //---
    EigenQuaternion<T> q;

    if (fabs(1 - (QA.dot(QB))) < std::numeric_limits<T>::epsilon() * 10)
    {
        q = EigenQuaternion<T>::Identity();
    }
    else { q = (((QB.conjugate()) * QA).normalized()); }
    EigenVector3<T> vecA = (TA - TB);
    return CoordSysEigen<T>(rotate((QB).conjugate(), vecA), q);
}

template <typename T>
inline CoordSysEigen<T> model_update(CoordSysEigen<T> const& A,
                                     CoordSysEigen<T> const& B)
{
    return model_update(A.T(), A.Q(), B.T(), B.Q());
}

template <typename T>
inline EigenVector3<T> xform_point(CoordSysEigen<T> const& X,
                                   EigenVector3<T> const& p)
{
    return rotate(X.Q(), p) + X.T();
}

template <typename T>
inline EigenVector3<T> xform_vector(CoordSysEigen<T> const& X,
                                    EigenVector3<T> const& v)
{
    return rotate(X.Q(), v);
}

template <typename T>
inline CoordSysEigen<T> prod(CoordSysEigen<T> const& L,
                             CoordSysEigen<T> const& R)
{
    return CoordSysEigen<T>(rotate(L.Q(), R.T()) + L.T(),
                            (prodEigen(L.Q(), R.Q())).normalized());
}

template <typename T> static T convert_to_degrees_eigen(T const& radians)
{
    return radians * T(57.295779513082320876798154814105);
}

template <typename T> static T convert_to_radians_eigen(T const& degrees)
{
    return degrees
         * T(0.017453292519943295769236907684886);
}

template <typename T>
inline EigenVector3<T> randomEigen(const T& lower, const T& upper)
{
    std::random_device rd;
    std::uniform_real_distribution<T> dist(lower, upper);

    EigenVector3<T> v(dist(rd), dist(rd), dist(rd));
    return v;
}

template <typename T> inline EigenVector3<T> randomEigen()
{
    return randomEigen<T>(0, 1);
}

namespace Eigen
{

template <std::size_t Index, typename T>
T& get(EigenVector3<T>& vec) noexcept
requires(Index < 3)
{
    return vec.data()[Index];
}

template <std::size_t Index, typename T>
T&& get(EigenVector3<T>&& vec) noexcept
requires(Index < 3)
{
    return std::move(vec.data()[Index]);
}

template <std::size_t Index, typename T>
const T& get(const EigenVector3<T>& vec) noexcept
requires(Index < 3)
{
    return vec.data()[Index];
}

template <std::size_t Index, typename T>
const T&& get(const EigenVector3<T>&& vec) noexcept
requires(Index < 3)
{
    return std::move(vec.data()[Index]);
}

} // namespace Eigen

template <typename T>
struct std::tuple_size<EigenVector3<T>>
    : public std::integral_constant<std::size_t, 3>
{
};

template <std::size_t Index, typename T>
struct std::tuple_element<Index, EigenVector3<T>>
    : public std::integral_constant<std::size_t, 3>
{
    using type = T;
};

template <typename T>
void orthonormal_vectors(EigenVector3<T>& i, EigenVector3<T>& j,
                         const EigenVector3<T>& k)
{
    EigenVector3<T> absK(std::abs(k.x()), std::abs(k.y()), std::abs(k.z()));

    if (absK(0) > absK(1))
    {
        if (absK(0) > absK(2))
            i = {0, 1, 0};
        else
            i = {1, 0, 0};
    }
    else
    {
        if (absK(1) > absK(2))
            i = {0, 0, 1};
        else
            i = {1, 0, 0};
    }
    j = cross(k, i).normalized();
    i = cross(j, k);
}

template <typename T>
inline EigenQuaternion<T> slerp(EigenQuaternion<T> const& A,
                                EigenQuaternion<T> const& B, const T& w)
{
    using std::acos;
    using std::sin;

    assert(w >= 0 || !"slerp(): w must not be less than 0");
    assert(w <= 1 || !"slerp(): w must not be larger than 1");

    T const q_tiny = (10e-7);
    T norm = dot(A, B);

    bool flip = false;
    if (norm < 0)
    {
        norm = -norm;
        flip = true;
    }
    T weight = w;
    T inv_weight;
    if (1 - norm < q_tiny) { inv_weight = 1 - weight; }
    else
    {
        T const theta = T(acos(norm));
        T const s_val = T(1 / sin(theta));
        inv_weight = T(sin((1 - weight) * theta) * s_val);
        weight = T(sin(weight * theta) * s_val);
    }
    if (flip) { weight = -weight; }
    Eigen::Quaternion<T> res = Eigen::Quaternion<T>(
        (Eigen::Quaternion<T>(inv_weight * A.coeffs())).coeffs()
        + (Eigen::Quaternion<T>(weight * B.coeffs())).coeffs());
    return res;
}

#ifdef WIN32
#define is_number(val) (_isnan(val) == 0) ///< Is a number test
#else
#if (__APPLE__)
#define is_number(val) (std::isnan(val) == 0) ///< Is a number test
#else
#ifdef __linux__
#define is_number(val) (std::isnan(val) == false) ///< Is a number test
#else
#define is_number(val) (isnan(val) == 0) ///< Is a number test
#endif
#endif
#endif

#ifdef WIN32
#define is_finite(val) (_finite(val) != 0) ///< Is finite number test
#else
#if (__APPLE__)
#define is_finite(val) std::isfinite(val) ///<Is finite number test
#else
#define is_finite(val) (finite(val) != 0) ///< Is finite number test
#endif
#endif

template <typename T>
inline void
eigen(const Eigen::Matrix<T, 3, 3>& A,
      Eigen::Matrix<T, 3, 3>& V, // output: eigenvectors as columns
      Eigen::Matrix<T, 3, 1>& diag) // output: eigenvalues (3x1 vector)
{
    // Ensure symmetric input (the original algorithm assumes symmetric matrix).
    Eigen::Matrix<T, 3, 3> sym = (A + A.transpose()) * T(0.5);

    // Solve
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T, 3, 3>> es;
    es.compute(sym); // cheap for 3x3

    // Copy outputs
    diag = es.eigenvalues(); // ascending order by default
    V = es.eigenvectors(); // columns are eigenvectors
}
// EIGENHELPERS_H
#endif
