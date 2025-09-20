#ifndef GRID_HELPERS_H
#define GRID_HELPERS_H

#include "grid_node_position.h"
#include "igl/signed_distance.h"
#include "types/geometry_triangle.h"
#include <filesystem>
#include <grid_grid.h>
#include <mesh_array_vertex_attribute.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <ostream>
#include <fstream>
#include <random>

namespace grid
{

template <typename T> void rotate_Zup_to_Yup(Eigen::MatrixXd& V)
{
    const double a = -M_PI / 2.0; // inverse
    const double c = std::cos(a);
    const double s = std::sin(a);
    Eigen::Matrix3d R;
    R << 1, 0, 0, 0, c, -s, 0, s, c;
    V = (R * V.transpose()).transpose();
}

template <typename T> void rotate_Yup_to_Zup(Eigen::MatrixXd& V)
{
    const double a = M_PI / 2.0; // +90deg
    const double c = std::cos(a);
    const double s = std::sin(a);
    Eigen::Matrix3d R;
    R << 1, 0, 0, 0, c, -s, 0, s, c;
    // V is #V x 3; apply R to each row:
    V = (R * V.transpose()).transpose();
}

template <typename T> void rot_Xup_to_Zup(Eigen::MatrixXd& V)
{
    // Rotate -90° about Y: X -> Z
    double a = -M_PI / 2.0;
    double c = std::cos(a), s = std::sin(a);
    Eigen::Matrix3d R;
    R << c, 0, s, 0, 1, 0, -s, 0, c;
    V = (R * V.transpose()).transpose();
}

template <typename T> void rot_Xup_to_Zup_and_rotX(Eigen::MatrixXd& V)
{
    // Step 1: Rotate -90° about Y: X -> Z
    double a = -M_PI / 2.0;
    double c = std::cos(a), s = std::sin(a);
    Eigen::Matrix3d R_y;
    R_y << c, 0, s, 0, 1, 0, -s, 0, c;

    // Step 2: Rotate +90° about X
    double b = M_PI / 2.0;
    double cx = std::cos(b), sx = std::sin(b);
    Eigen::Matrix3d R_x;
    R_x << 1, 0, 0, 0, cx, -sx, 0, sx, cx;

    // Combine rotations: R = R_x * R_y
    Eigen::Matrix3d R = R_x * R_y;

    // Apply to vertices
    V = (R * V.transpose()).transpose();
}

static inline double deg2rad(double deg) { return deg * M_PI / 180.0; }

// Ensure matrix has 3 columns (x,y,z)
template <typename T> static inline void ensure_is_xyz(const Eigen::MatrixXd& V)
{
    if (V.cols() != 3)
        throw std::invalid_argument(
            "Vertex matrix must have 3 columns (x,y,z)");
}

// Rotate in-place: N x 3 matrix, degrees can be any real number
template <typename T> inline void rotateX(Eigen::MatrixXd& V, double degrees)
{
//    ensure_is_xyz<T>(V);
    double a = deg2rad(degrees);
    double c = std::cos(a), s = std::sin(a);
    Eigen::Matrix3d R;
    R << 1, 0, 0, 0, c, -s, 0, s, c;
    V = (R * V.transpose()).transpose();
}

template <typename T> inline void rotateY(Eigen::MatrixXd& V, double degrees)
{
   // ensure_is_xyz(V);
    double a = deg2rad(degrees);
    double c = std::cos(a), s = std::sin(a);
    Eigen::Matrix3d R;
    R << c, 0, s, 0, 1, 0, -s, 0, c;
    V = (R * V.transpose()).transpose();
}

template <typename T> inline void rotateZ(Eigen::MatrixXd& V, double degrees)
{
  //  ensure_is_xyz(V);
    double a = deg2rad(degrees);
    double c = std::cos(a), s = std::sin(a);
    Eigen::Matrix3d R;
    R << c, -s, 0, s, c, 0, 0, 0, 1;
    V = (R * V.transpose()).transpose();
}

template <typename D, typename T>
Grid<D, T> projectGridToSDF(Eigen::MatrixXd verts, Eigen::MatrixXi indices,
                            Eigen::Matrix<size_t, 3, 1> res)
{
//    rotate_Yup_to_Zup<T>(verts);
rot_Xup_to_Zup<T>(verts);
//    rotate_Zup_to_Yup<T>(verts);
//rotateY<T>(verts, -90.0);
//    rot_Xup_to_Zup_and_rotX<T>(verts);
Eigen::RowVector3d minv = verts.colwise().minCoeff();
Eigen::RowVector3d maxv = verts.colwise().maxCoeff();
Eigen::RowVector3d diag = maxv - minv;
double longest = diag.maxCoeff();
    //10% padding to our bounding box!
double pad = 0.10 * longest;
Eigen::Matrix<T, 3, 1> gmin((T)(minv.x() - pad), (T)(minv.y() - pad),
                            (T)(minv.z() - pad));
Eigen::Matrix<T, 3, 1> gmax((T)(maxv.x() + pad), (T)(maxv.y() + pad),
                            (T)(maxv.z() + pad));

    //create grid (res^3 nodes)
Eigen::Matrix<size_t, 3, 1> nodes((size_t)res.x(), (size_t)res.y(),
                                  (size_t)res.z());
Grid<D, T> G;
G.create(gmin, gmax, nodes);
std::cerr << "GMIN " << gmin << "\n";
std::cerr << "GMAX " << gmax << "\n";
const size_t total = G.m_nodes.x() * G.m_nodes.y() * G.m_nodes.z();
std::cout << "Created grid: " << G.I() << " x " << G.J() << " x " << G.K()
          << "  (total nodes = " << total << ")\n";

//Build the query points matrix P (total x 3) in the same linear order used by grid
Eigen::MatrixXd P((Eigen::Index)total, 3);
size_t idx_lin = 0;
for (size_t k = 0; k < G.K(); ++k)
{
    for (size_t j = 0; j < G.J(); ++j)
    {
        for (size_t i = 0; i < G.I(); ++i)
        {
            Eigen::Matrix<size_t, 3, 1> idx(i, j, k);
            Eigen::Matrix<T, 3, 1> p;
            grid::node_position(G, idx, p);
            P((Eigen::Index)idx_lin, 0) = p.x();
            P((Eigen::Index)idx_lin, 1) = p.y();
            P((Eigen::Index)idx_lin, 2) = p.z();
            ++idx_lin;
        }
    }
}

    std::cout << "Computing signed distances (libigl::signed_distance)...\n";
    Eigen::VectorXd S; // signed distances
    Eigen::VectorXi I; // indices to closest triangles (unused here)
    Eigen::MatrixXd C; // closest points on triangles
    Eigen::MatrixXd N; // normals used for signing (pseudonormal)

    igl::SignedDistanceType signType = igl::SIGNED_DISTANCE_TYPE_WINDING_NUMBER;
    igl::signed_distance(P, verts, indices, signType, S, I, C, N);

    for (size_t s = 0; s < total; ++s)
    {
        G.data()[s] = static_cast<D>(S((Eigen::Index)s));
    }
    return G;
}

template <typename T>
void build_VF_from_T3Mesh(
    mesh_array::T3Mesh const& mesh,
    mesh_array::VertexAttribute<T, mesh_array::T3Mesh> const& X,
    mesh_array::VertexAttribute<T, mesh_array::T3Mesh> const& Y,
    mesh_array::VertexAttribute<T, mesh_array::T3Mesh> const& Z,
    Eigen::MatrixXd& V, // out
    Eigen::MatrixXi& F // out (libigl expects int indices)
)
{
    const size_t nV = mesh.vertex_size();
    const size_t nF = mesh.triangle_size();

    V.resize((int)nV, 3);
    F.resize((int)nF, 3);

    // Fill vertex positions (cast to double because libigl uses double)
    for (size_t i = 0; i < nV; ++i)
    {
        auto const& vert = mesh.vertex(i);
        // operator()(Vertex) returns T (or reference), cast to double
        V((int)i, 0) = static_cast<double>(X(vert));
        V((int)i, 1) = static_cast<double>(Y(vert));
        V((int)i, 2) = static_cast<double>(Z(vert));
    }

    // Fill triangle indices (int). We'll assume triangles expose vertex(k).idx().
    for (size_t t = 0; t < nF; ++t)
    {
        auto const& tri = mesh.triangle(t);

        // --- Primary (likely) pattern: triangle.vertex(j).idx()
        int a = static_cast<int>(tri.i());
        int b = static_cast<int>(tri.j());
        int c = static_cast<int>(tri.k());

        F((int)t, 0) = a;
        F((int)t, 1) = b;
        F((int)t, 2) = c;

        // ----- If the above doesn't compile, try one of these (pick and replace) -----
        // Variant: triangle provides integer indices directly: tri[0], tri[1], tri[2]
        // int a = static_cast<int>( tri[0] );
        // int b = static_cast<int>( tri[1] );
        // int c = static_cast<int>( tri[2] );

        // Variant: triangle exposes members e.g. v0, v1, v2 which are Vertex
        // int a = static_cast<int>( tri.v0.idx() );
        // int b = static_cast<int>( tri.v1.idx() );
        // int c = static_cast<int>( tri.v2.idx() );

        // Variant: triangle.operator[](j) returns Vertex
        // int a = static_cast<int>( tri[0].idx() );
        // int b = static_cast<int>( tri[1].idx() );
        // int c = static_cast<int>( tri[2].idx() );
    }
}

template <typename T>
std::vector<GridTriangle<T>> build_triangle_list_from_T3Mesh(
    mesh_array::T3Mesh const& mesh,
    mesh_array::VertexAttribute<T, mesh_array::T3Mesh> const& X,
    mesh_array::VertexAttribute<T, mesh_array::T3Mesh> const& Y,
    mesh_array::VertexAttribute<T, mesh_array::T3Mesh> const& Z)
{
    size_t nF = mesh.triangle_size();
    std::vector<GridTriangle<T>> tris;
    tris.reserve(nF);

    for (size_t t = 0; t < nF; ++t)
    {
        auto const& tri = mesh.triangle(t);

        // primary assumption: tri.i(), tri.j(), tri.k() give integer vertex indices
        int ia = static_cast<int>(tri.i());
        int ib = static_cast<int>(tri.j());
        int ic = static_cast<int>(tri.k());

        auto va = mesh.vertex((size_t)ia);
        auto vb = mesh.vertex((size_t)ib);
        auto vc = mesh.vertex((size_t)ic);

        GridTriangle<T> tri_xyz;
        tri_xyz.v0 << X(va), Y(va), Z(va);
        tri_xyz.v1 << X(vb), Y(vb), Z(vb);
        tri_xyz.v2 << X(vc), Y(vc), Z(vc);

        tris.push_back(tri_xyz);
    }

    return tris;
}

template <typename D, typename T>
void writeGridToFiles(const Grid<D, T>& gridType, size_t globalRigidBodyIdx)
{
    uint32_t a = 0;
    uint32_t b = 1000;
    std::random_device rd; // Will seed the random number generator
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint32_t> dis(a, b);
    uint32_t random_value = dis(gen);
    globalRigidBodyIdx = size_t(random_value);
    namespace fs = std::filesystem;
    fs::path source_dir = fs::current_path();
    std::cerr << source_dir;
    fs::path bunnySDFRelative
        = "../output/sdftest/sdf" + std::to_string(globalRigidBodyIdx) + ".csv";
    fs::path bunnySDFFull = source_dir / bunnySDFRelative;
    fs::path bunnySDFNormalized = bunnySDFFull.lexically_normal();

    const std::string outCsv = bunnySDFNormalized.string();

    std::ofstream csv(outCsv);
    if (!csv)
    {
        std::cerr << "ERROR: cannot open " << outCsv << " for writing\n";
    }
    csv << "x,y,z,sdf\n";

    D maxZ = 0.0;
    for (size_t k = 0; k < gridType.K(); ++k)
    {
        for (size_t j = 0; j < gridType.J(); ++j)
        {
            for (size_t i = 0; i < gridType.I(); ++i)
            {
                Eigen::Matrix<size_t, 3, 1> idx(i, j, k);
                Eigen::Matrix<T, 3, 1> p;
                grid::node_position(gridType, idx, p);
                csv << std::setprecision(9) << p.x() << "," << p.y() << ","
                    << p.z() << "," << gridType(idx) << "\n";
                maxZ = std::max<D>(p.z(), maxZ);
            }
        }
    }
    std::cerr << "MAXZ: " << maxZ << "\n";
    csv.close();
}

//void computeVerticesFacesFromMesh(mesh_array::T3Mesh mesh) {}
} // namespace grid

#endif // GRID_HELPERS_H
