#include <mesh_array_make_convex.h>

#if defined(__cplusplus)
extern "C"
{
#endif
#include <stdio.h>
#include <stdlib.h>
#include <qhull.h>
#include <mem.h>
#include <qset.h>
#include <geom.h>
#include <merge.h>
#include <poly.h>
#include <io.h>
#include <stat.h>
#if defined(__cplusplus)
}
#endif

#include <vector>
#include <cassert>

namespace mesh_array
{

template <typename T>
void make_convex_2(std::vector<EigenVector3<T>> const& vertices, T3Mesh& mesh,
                   VertexAttribute<T, T3Mesh>& X, VertexAttribute<T, T3Mesh>& Y,
                   VertexAttribute<T, T3Mesh>& Z)
{
    size_t const N = vertices.size();

    assert(N >= 3u
           || !"make_convex(): Must be at least 3 vertices for a triangle fan");

    mesh.clear();

    // Allocate space for vertices and triangles
    // For a triangle fan with N vertices, we'll have (N-2) triangles
    size_t num_triangles = (N >= 3) ? (N - 2) : 0;
    mesh.set_capacity(N, num_triangles);

    X.bind(mesh);
    Y.bind(mesh);
    Z.bind(mesh);

    // Add all vertices to the mesh
    for (size_t i = 0u; i < N; ++i)
    {
        mesh.push_vertex();
        Vertex v = mesh.vertex(i);
        X(v) = vertices[i](0);
        Y(v) = vertices[i](1);
        Z(v) = vertices[i](2);
    }

    // Create triangle fan: (0,1,2), (0,2,3), (0,3,4), etc.
    for (size_t i = 2u; i < N; ++i)
    {
        Vertex v0 = mesh.vertex(0);
        Vertex v1 = mesh.vertex(i - 1);
        Vertex v2 = mesh.vertex(i);
        mesh.push_triangle(v0, v1, v2);
    }
}

template <typename T>
void make_convex(std::vector<EigenVector3<T>> const& vertices, T3Mesh& mesh, VertexAttribute<T, T3Mesh>& X,
                 VertexAttribute<T, T3Mesh>& Y, VertexAttribute<T, T3Mesh>& Z)
{
    size_t const N = vertices.size();

    assert(N > 0u || !"make_convex(): Must be a positive number of vertices");

    std::vector<coordT> coords;
    coords.resize(3u * N);
    for (size_t i = 0u; i < N; ++i)
    {
        coords[i * 3u] = vertices[i](0);
        coords[i * 3u + 1u] = vertices[i](1);
        coords[i * 3u + 2u] = vertices[i](2);
    }

    //------------------------------------------------

    boolT ismalloc = False;
    char flags[] = "qhull Qx Qt"; //QJ maybe Qt
    FILE* outfile = stdout;    //FILE *outfile = 0;
    FILE* errfile = stderr;    //FILE *errfile = 0;

    int K = N;
    int dim = 3;
    int exitcode = qh_new_qhull(dim, K, &coords[0], ismalloc, flags, outfile, errfile);

    if (!exitcode)
    {
      //--- Determine how many triangle faces we have in hull
        size_t F = 0u;
        {
            facetT* facet;
            qh visit_id++;
            FORALLfacets
            {
                facet->visitid = qh visit_id;
                ++F;
            }
        }

      //--- Allocate space for the convex hull mesh
        mesh.clear();
        // We'll allocate conservatively: number of hull vertices <= N, faces = F
        mesh.set_capacity(N, F); // still ok as an upper bound

        X.bind(mesh);
        Y.bind(mesh);
        Z.bind(mesh);

        // ---- collect unique point ids actually used by facets ----
        std::vector<int> used_point_ids; // list of used qh point ids
        used_point_ids.reserve(N);
        std::vector<int> id_to_newindex(
            N, -1); // map original qh point id -> index in mesh (or -1)

        {
            facetT* facet;
            vertexT *vertex, **vertexp;
            qh visit_id++;
            FORALLfacets
            {
                facet->visitid = qh visit_id;
                FOREACHvertex_(facet->vertices)
                {
                    int orig_id = qh_pointid(vertex->point);
                    if (id_to_newindex[orig_id] == -1)
                    {
                        id_to_newindex[orig_id]
                            = static_cast<int>(used_point_ids.size());
                        used_point_ids.push_back(orig_id);
                    }
                }
            }
        }

        // --- push only the used hull vertices into the mesh, in the order of used_point_ids ---
        for (size_t k = 0; k < used_point_ids.size(); ++k)
        {
            int orig_id = used_point_ids[k];
            mesh.push_vertex();
            Vertex v = mesh.vertex(static_cast<size_t>(k));
            // coords array: coords[3*orig_id + 0..2]
            X(v) = coords[3u * orig_id + 0u];
            Y(v) = coords[3u * orig_id + 1u];
            Z(v) = coords[3u * orig_id + 2u];
        }

        // --- now push triangles, remapping qh_pointid -> new mesh index ---
        {
            facetT* facet;
            vertexT *vertex, **vertexp;
            qh visit_id++;
            FORALLfacets
            {
                facet->visitid = qh visit_id;

                int i = 0;
                Vertex idx[3];

                FOREACHvertex_(facet->vertices)
                {
                    int orig_id = qh_pointid(vertex->point);
                    int new_idx = id_to_newindex[orig_id];
                    idx[i++] = mesh.vertex(static_cast<size_t>(new_idx));
                }

                if (facet->toporient)
                {
                    mesh.push_triangle(idx[2], idx[1], idx[0]);
                }
                else { mesh.push_triangle(idx[0], idx[1], idx[2]); }
            }
        }
    }

    qh_freeqhull(!qh_ALL);

    int curlong;
    int totlong;
    qh_memfreeshort(&curlong, &totlong);

    if (curlong || totlong)
    {
        fprintf(errfile, "qhull internal warning (main): did not free %d bytes of long memory (%d pieces)\n", totlong,
                curlong);
    }
}

template void make_convex<float>(std::vector<EigenVector3<float>> const& vertices, T3Mesh& mesh,
                                 VertexAttribute<float, T3Mesh>& X, VertexAttribute<float, T3Mesh>& Y,
                                 VertexAttribute<float, T3Mesh>& Z);

template void make_convex<double>(std::vector<EigenVector3<double>> const& vertices, T3Mesh& mesh,
                                  VertexAttribute<double, T3Mesh>& X, VertexAttribute<double, T3Mesh>& Y,
                                  VertexAttribute<double, T3Mesh>& Z);

template void
make_convex_2<float>(std::vector<EigenVector3<float>> const& vertices,
                     T3Mesh& mesh, VertexAttribute<float, T3Mesh>& X,
                     VertexAttribute<float, T3Mesh>& Y,
                     VertexAttribute<float, T3Mesh>& Z);

template void
make_convex_2<double>(std::vector<EigenVector3<double>> const& vertices,
                      T3Mesh& mesh, VertexAttribute<double, T3Mesh>& X,
                      VertexAttribute<double, T3Mesh>& Y,
                      VertexAttribute<double, T3Mesh>& Z);

} //namespace mesh_array
