#include <iostream>
#include <mass.h>

#include <eigenhelperall.h>

#include <cassert>
#include <cmath>
#include <vector>

namespace mass
{
namespace details
{

template <typename T> class FaceInfo
{
protected:
    std::vector<T> m_coord[3];

    T m_n[3] = {T(0), T(1), T(0)};

    size_t m_A;
    size_t m_B;
    size_t m_C;

public:
    T const* get_a_ptr() const { return &(m_coord[m_A][0]); }

    T const* get_b_ptr() const { return &(m_coord[m_B][0]); }

    T const& get_nA() const { return m_n[m_A]; }

    T const& get_nB() const { return m_n[m_B]; }

    T const& get_nC() const { return m_n[m_C]; }

    T const& get_nx() const { return m_n[0]; }

    T const& get_ny() const { return m_n[1]; }

    T const& get_nz() const { return m_n[2]; }

    T compute_w() const { return -(m_n[0] * m_coord[0][0] + m_n[1] * m_coord[1][0] + m_n[2] * m_coord[2][0]); }

    size_t const& get_A() const { return m_A; }

    size_t const& get_B() const { return m_B; }

    size_t const& get_C() const { return m_C; }

    size_t get_number_of_vertices() const { return m_coord[0].size(); }

public:
    FaceInfo() {}

    void init(size_t const& N, T const* x, T const* y, T const* z)
    {
        using std::fabs;
        using std::sqrt;


        assert(N >= 3u || !"N was less than 3");
        assert(x || !"x coord pointer was null");
        assert(y || !"y coord pointer was null");
        assert(z || !"z coord pointer was null");

          //--- Copy data into internal data container
        m_coord[0].resize(N);
        m_coord[1].resize(N);
        m_coord[2].resize(N);

        for (size_t i = 0u; i < N; ++i)
        {
            assert(is_number(x[i]) || !"x[i] must be a number");
            assert(is_finite(x[i]) || !"x[i] must be a finite number");
            assert(is_number(y[i]) || !"y[i] must be a number");
            assert(is_finite(y[i]) || !"y[i] must be a finite number");
            assert(is_number(z[i]) || !"z[i] must be a number");
            assert(is_finite(z[i]) || !"z[i] must be a finite number");

            m_coord[0][i] = x[i];
            m_coord[1][i] = y[i];
            m_coord[2][i] = z[i];
        }

          //--- Find face normal from the cross-product of the two edges that span the largest area
        T max_area = 0;
        m_n[0] = 0;
        m_n[1] = 1;
        m_n[2] = 0;

        for (size_t i = 0u; i < N; ++i)
        {
            size_t const j = (i + 1) % N;
            size_t const k = (j + 1) % N;

            T const dx_1 = x[j] - x[i];
            T const dy_1 = y[j] - y[i];
            T const dz_1 = z[j] - z[i];
            /*std::cerr << "i, j: (" << i << ", " << j << ")  and x[j], x[i] = ("
                      << x[j] << " , " << x[i] << ")\n";
            std::cerr << "i, j: (" << i << ", " << j << ")  and y[j], y[i] = ("
                      << y[j] << " , " << y[i] << ")\n";
            std::cerr << "i, j: (" << i << ", " << j << ")  and z[j], z[i] = ("
                      << z[j] << " , " << z[i] << ")\n";*/

            T const dx_2 = x[k] - x[j];
            T const dy_2 = y[k] - y[j];
            T const dz_2 = z[k] - z[j];
            /*std::cerr << "k, j: (" << i << ", " << j << ")  and x[k], x[j] = ("
                      << x[k] << " , " << x[j] << ")\n";
            std::cerr << "k, j: (" << i << ", " << j << ")  and y[k], y[j] = ("
                      << y[k] << " , " << y[j] << ")\n";
            std::cerr << "k, j: (" << i << ", " << j << ")  and z[k], z[j] = ("
                      << z[k] << " , " << z[j] << ")\n";*/

            T const nx = dy_1 * dz_2 - dy_2 * dz_1;
            T const ny = dx_2 * dz_1 - dx_1 * dz_2;
            T const nz = dx_1 * dy_2 - dx_2 * dy_1;

            T const area = nx * nx + ny * ny + nz * nz;

            if (area > max_area)
            {
                max_area = area;
                m_n[0] = nx;
                m_n[1] = ny;
                m_n[2] = nz;
            }
            /*if (max_area < 0)
            {
                std::cerr << "MAX AREA TOO SMALL: " << max_area << "\n";
            }*/
        }

          //--- Normalize the best triangle normal and use it as face normal
        T const l = sqrt(max_area);
        m_n[0] /= l;
        m_n[1] /= l;
        m_n[2] /= l;

        //--- Find the best coordinate plane to project the face onto -- that is the coordinate plane that maximizes the projection area
        T const abs_nx = fabs(m_n[0]);
        T const abs_ny = fabs(m_n[1]);
        T const abs_nz = fabs(m_n[2]);

        enum
        {
            X = 0,
            Y = 1,
            Z = 2
        }; // Just some convenience constants

        m_C = (abs_ny > abs_nz) ? Y : Z;

        if (abs_nx > abs_ny && abs_nx > abs_nz) m_C = X;

        m_A = (m_C + 1) % 3;
        m_B = (m_C + 2) % 3;
    }
};

template <typename T> class ProjectionIntegralInfo
{

public:
    T m_P1 = T(0.0);
    T m_Pa = T(0.0);
    T m_Pb = T(0.0);
    T m_Paa = T(0.0);
    T m_Pab = T(0.0);
    T m_Pbb = T(0.0);
    T m_Paaa = T(0.0);
    T m_Paab = T(0.0);
    T m_Pabb = T(0.0);
    T m_Pbbb = T(0.0);

public:
    ProjectionIntegralInfo()
        : m_P1(0)
        , m_Pa(0)
        , m_Pb(0)
        , m_Paa(0)
        , m_Pab(0)
        , m_Pbb(0)
        , m_Paaa(0)
        , m_Paab(0)
        , m_Pabb(0)
        , m_Pbbb(0)
    {
    }
};

template <typename T> ProjectionIntegralInfo<T> compute_projected_face_integral(FaceInfo<T> const& face)
{
    using std::fabs;

    ProjectionIntegralInfo<T> result;

    T const* a = face.get_a_ptr();
    T const* b = face.get_b_ptr();
    size_t const N = face.get_number_of_vertices();

    for (size_t i = 0u; i < N; ++i)
    {
        T const a0 = a[i];
        T const b0 = b[i];
        T const a1 = a[(i + 1) % N];
        T const b1 = b[(i + 1) % N];

        T const da = a1 - a0;
        T const db = b1 - b0;
        T const a0_2 = a0 * a0;
        T const a0_3 = a0_2 * a0;
        T const a0_4 = a0_3 * a0;
        T const b0_2 = b0 * b0;
        T const b0_3 = b0_2 * b0;
        T const b0_4 = b0_3 * b0;
        T const a1_2 = a1 * a1;
        T const a1_3 = a1_2 * a1;
        T const b1_2 = b1 * b1;
        T const b1_3 = b1_2 * b1;
        T const C1 = a1 + a0;
        T const Ca = a1 * C1 + a0_2;
        T const Caa = a1 * Ca + a0_3;
        T const Caaa = a1 * Caa + a0_4;
        T const Cb = b1 * (b1 + b0) + b0_2;
        T const Cbb = b1 * Cb + b0_3;
        T const Cbbb = b1 * Cbb + b0_4;
        T const Cab = 3 * a1_2 + 2 * a1 * a0 + a0_2;
        T const Kab = a1_2 + 2 * a1 * a0 + 3 * a0_2;
        T const Caab = a0 * Cab + 4 * a1_3;
        T const Kaab = a1 * Kab + 4 * a0_3;
        T const Cabb = 4 * b1_3 + 3 * b1_2 * b0 + 2 * b1 * b0_2 + b0_3;
        T const Kabb = b1_3 + 2 * b1_2 * b0 + 3 * b1 * b0_2 + 4 * b0_3;

        result.m_P1 += db * C1;
        result.m_Pa += db * Ca;
        result.m_Paa += db * Caa;
        result.m_Paaa += db * Caaa;
        result.m_Pb += da * Cb;
        result.m_Pbb += da * Cbb;
        result.m_Pbbb += da * Cbbb;
        result.m_Pab += db * (b1 * Cab + b0 * Kab);
        result.m_Paab += db * (b1 * Caab + b0 * Kaab);
        result.m_Pabb += da * (a1 * Cabb + a0 * Kabb);
    }

    result.m_P1 /= 2;
    result.m_Pa /= T(6.0);
    result.m_Paa /= T(12.0);
    result.m_Paaa /= T(20.0);
    result.m_Pb /= -T(6.0);
    result.m_Pbb /= -T(12.0);
    result.m_Pbbb /= -T(20.0);
    result.m_Pab /= T(24.0);
    result.m_Paab /= T(60.0);
    result.m_Pabb /= -T(60.0);

    return result;
}

template <typename T> class FaceIntegralInfo
{

public:
    T m_Fa;
    T m_Fb;
    T m_Fc;
    T m_Faa;
    T m_Fbb;
    T m_Fcc;
    T m_Faaa;
    T m_Fbbb;
    T m_Fccc;
    T m_Faab;
    T m_Fbbc;
    T m_Fcca;

public:
    FaceIntegralInfo()
        : m_Fa(0)
        , m_Fb(0)
        , m_Fc(0)
        , m_Faa(0)
        , m_Fbb(0)
        , m_Fcc(0)
        , m_Faaa(0)
        , m_Fbbb(0)
        , m_Fccc(0)
        , m_Faab(0)
        , m_Fbbc(0)
        , m_Fcca(0)
    {
    }
};

template <typename T> FaceIntegralInfo<T> compute_face_integral(FaceInfo<T> const& face)
{
    using std::fabs;


    FaceIntegralInfo<T> result;

    T const w = face.compute_w();
    T const& nA = face.get_nA();
    T const& nB = face.get_nB();
    T const& nC = face.get_nC();
    T const k1 = 1 / nC;
    T const k2 = k1 * k1;
    T const k3 = k2 * k1;
    T const k4 = k3 * k1;

    ProjectionIntegralInfo<T> term = compute_projected_face_integral(face);

    result.m_Fa = k1 * term.m_Pa;
    result.m_Fb = k1 * term.m_Pb;
    result.m_Fc = -k2 * (nA * term.m_Pa + nB * term.m_Pb + w * term.m_P1);

    result.m_Faa = k1 * term.m_Paa;
    result.m_Fbb = k1 * term.m_Pbb;
    result.m_Fcc = k3
                 * (nA * nA * term.m_Paa + 2 * nA * nB * term.m_Pab + nB * nB * term.m_Pbb
                    + w * (2 * (nA * term.m_Pa + nB * term.m_Pb) + w * term.m_P1));

    result.m_Faaa = k1 * term.m_Paaa;
    result.m_Fbbb = k1 * term.m_Pbbb;
    result.m_Fccc = -k4
                  * (nA * nA * nA * term.m_Paaa + 3 * nA * nA * nB * term.m_Paab + 3 * nA * nB * nB * term.m_Pabb
                     + nB * nB * nB * term.m_Pbbb
                     + 3 * w * (nA * nA * term.m_Paa + 2 * nA * nB * term.m_Pab + nB * nB * term.m_Pbb)
                     + w * w * (3 * (nA * term.m_Pa + nB * term.m_Pb) + w * term.m_P1));

    result.m_Faab = k1 * term.m_Paab;
    result.m_Fbbc = -k2 * (nA * term.m_Pabb + nB * term.m_Pbbb + w * term.m_Pbb);
    result.m_Fcca = k3
                  * (nA * nA * term.m_Paaa + 2 * nA * nB * term.m_Paab + nB * nB * term.m_Pabb
                     + w * (2 * (nA * term.m_Paa + nB * term.m_Pab) + w * term.m_Pa));

    return result;
}

template <typename T> class VolumeIntegralInfo
{
public:
    T m_T0;
    T m_T1[3];
    T m_T2[3];
    T m_TP[3];

public:
    VolumeIntegralInfo()
        : m_T0(0)
    {
        m_T1[0] = 0;
        m_T1[1] = 0;
        m_T1[2] = 0;
        m_T2[0] = 0;
        m_T2[1] = 0;
        m_T2[2] = 0;
        m_TP[0] = 0;
        m_TP[1] = 0;
        m_TP[2] = 0;
    }

public:
    void get_center_of_mass(T& x, T& y, T& z) const
    {
        x = m_T1[0] / m_T0;
        y = m_T1[1] / m_T0;
        z = m_T1[2] / m_T0;
    }

    void get_mass(T const& density, T& mass) const
    {
        assert(is_number(density) || !"density must be a number");
        assert(is_finite(density) || !"density must be a finite number");
        assert(density > 0 || !"density must be positive");

        mass = density * m_T0;
    }

    void get_inertia_tensor(T const& density, T& Ixx, T& Iyy, T& Izz, T& Ixy, T& Ixz, T& Iyz) const
    {
        assert(is_number(density) || !"density must be a number");
        assert(is_finite(density) || !"density must be a finite number");
        assert(density > 0 || !"density must be positive");

        Ixx = density * (m_T2[1] + m_T2[2]);
        Iyy = density * (m_T2[2] + m_T2[0]);
        Izz = density * (m_T2[0] + m_T2[1]);
        Ixy = -density * m_TP[0];
        Iyz = -density * m_TP[1];
        Ixz = -density * m_TP[2];
    }
};

template <typename T>
VolumeIntegralInfo<T> compute_volume_integral(size_t const K, FaceCallbackInterface<T> const* callback)
{
    using std::fabs;


    VolumeIntegralInfo<T> result;

    //--- Loop over all faces
    for (size_t k = 0u; k < K; ++k)
    {
        // Extract/Compute geometry information about the current face iterate
        FaceInfo<T> face;

        size_t const N = callback->get_number_of_vertices(k);

        std::vector<T> xcoords(N);
        std::vector<T> ycoords(N);
        std::vector<T> zcoords(N);

        callback->get_x_coords(k, &xcoords[0]);
        callback->get_y_coords(k, &ycoords[0]);
        callback->get_z_coords(k, &zcoords[0]);

        face.init(N, &xcoords[0], &ycoords[0], &zcoords[0]);

        size_t const& A = face.get_A();
        size_t const& B = face.get_B();
        size_t const& C = face.get_C();
        T const& nA = face.get_nA();
        T const& nB = face.get_nB();
        T const& nC = face.get_nC();
        T const& nx = face.get_nx();

        FaceIntegralInfo<T> part = compute_face_integral(face);

        result.m_T0 += nx * ((A == 0) ? part.m_Fa : ((B == 0) ? part.m_Fb : part.m_Fc));

        result.m_T1[A] += nA * part.m_Faa;
        result.m_T1[B] += nB * part.m_Fbb;
        result.m_T1[C] += nC * part.m_Fcc;
        result.m_T2[A] += nA * part.m_Faaa;
        result.m_T2[B] += nB * part.m_Fbbb;
        result.m_T2[C] += nC * part.m_Fccc;
        result.m_TP[A] += nA * part.m_Faab;
        result.m_TP[B] += nB * part.m_Fbbc;
        result.m_TP[C] += nC * part.m_Fcca;
    }

    result.m_T1[0] /= 2;
    result.m_T1[1] /= 2;
    result.m_T1[2] /= 2;
    result.m_T2[0] /= 3;
    result.m_T2[1] /= 3;
    result.m_T2[2] /= 3;
    result.m_TP[0] /= 2;
    result.m_TP[1] /= 2;
    result.m_TP[2] /= 2;

    return result;
}

} // namespace details

template <typename T>
Properties<T> compute_mesh(T const& density, size_t const K, FaceCallbackInterface<T> const* callback)
{

    assert(is_number(density) || !"density must be a number");
    assert(is_finite(density) || !"density must be a finite number");
    assert(density > 0 || !"density must be positive");
    assert(K > 0u || !"K must be positive");
    details::VolumeIntegralInfo<T> result = details::compute_volume_integral(K, callback);

    Properties<T> value;

    result.get_mass(density, value.m_m);
    result.get_center_of_mass(value.m_x, value.m_y, value.m_z);
    result.get_inertia_tensor(density, value.m_Ixx, value.m_Iyy, value.m_Izz, value.m_Ixy, value.m_Ixz, value.m_Iyz);

    return value;
}

template Properties<float> compute_mesh<float>(float const& density, size_t const K,
                                               FaceCallbackInterface<float> const* callback);

template Properties<double> compute_mesh<double>(double const& density, size_t const K,
                                                 FaceCallbackInterface<double> const* callback);

} // namespace mass
