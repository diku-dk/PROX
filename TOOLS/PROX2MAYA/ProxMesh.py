import sys


class T4Mesh:
    def __init__(self, name, number_of_vertices, number_of_tetrahedra):
        self.name = name
        self.vertices = [(sys.float_info.max, sys.float_info.max, sys.float_info.max) for _ in xrange(number_of_vertices)]
        self.tetrahedra = [(-1, -1, -1, -1) for _ in xrange(number_of_tetrahedra)]


class T3Mesh:
    def __init__(self, name):
        self.name = name
        self.vertices = []
        self.triangles = []


def make_unique(indices):
    """
    Performs cyclic permutation of the index tuple(i, j, k), such that the first index is the one with the lowest value.

    :param indices:
    :return:
    """
    i = indices[0]
    j = indices[1]
    k = indices[2]

    if i < j and i < k:
        return i, j, k

    if j < i and j < k:
        return j, k, i

    if k < i and k < j:
        return k, i, j

    raise RuntimeError('Internal error, triangle indices was flawed?')


def make_twin(indices):
    """
    Swaps to interfaces to give the tuple with opposite orientation. That is given a triangle with orientation
     (i,j,k) this function return the same triangle with orientation (j,i,k)

    :param indices:
    :return:
    """
    return indices[1], indices[0], indices[2]


def make_surface_mesh(t4mesh):
    """
    This function extracts the surface mesh of a tetrahedron mesh.

    :param t4mesh:      The tetrahedron mesh.
    :return:            A triangle mesh of the surface of the given tetrahedron mesh.
    """
    # First we build a data structure for looking up adjacent tetrahedrons
    adjacent_info = [(-1, -1, -1, -1) for _ in t4mesh.tetrahedra]

    lookup = dict()
    tetrahedron_idx = 0
    for tetrahedron in t4mesh.tetrahedra:
        i = tetrahedron[0]
        j = tetrahedron[1]
        k = tetrahedron[2]
        m = tetrahedron[3]

        tri_i = make_unique((j, k, m))
        tri_j = make_unique((i, m, k))
        tri_k = make_unique((i, j, m))
        tri_m = make_unique((j, i, k))

        twin_i = make_unique(make_twin(tri_i))
        twin_j = make_unique(make_twin(tri_j))
        twin_k = make_unique(make_twin(tri_k))
        twin_m = make_unique(make_twin(tri_m))

        lookup[twin_i] = tetrahedron_idx
        lookup[twin_j] = tetrahedron_idx
        lookup[twin_k] = tetrahedron_idx
        lookup[twin_m] = tetrahedron_idx

        tetrahedron_idx += 1

    tetrahedron_idx = 0
    for tetrahedron in t4mesh.tetrahedra:
        i = tetrahedron[0]
        j = tetrahedron[1]
        k = tetrahedron[2]
        m = tetrahedron[3]

        tri_i = make_unique((j, k, m))
        tri_j = make_unique((i, m, k))
        tri_k = make_unique((i, j, m))
        tri_m = make_unique((j, i, k))

        twin_i_idx = None
        twin_j_idx = None
        twin_k_idx = None
        twin_m_idx = None

        if tri_i in lookup:
            twin_i_idx = lookup[tri_i]

        if tri_j in lookup:
            twin_j_idx = lookup[tri_j]

        if tri_k in lookup:
            twin_k_idx = lookup[tri_k]

        if tri_m in lookup:
            twin_m_idx = lookup[tri_m]

        adjacent_info[tetrahedron_idx] = (twin_i_idx, twin_j_idx, twin_k_idx, twin_m_idx)
        tetrahedron_idx += 1

    # Second we extract all triangles not having an adjacent tetrahedron
    surface_triangles = []
    tetrahedron_idx = 0
    for tetrahedron in t4mesh.tetrahedra:
        i = tetrahedron[0]
        j = tetrahedron[1]
        k = tetrahedron[2]
        m = tetrahedron[3]

        neighbors = adjacent_info[tetrahedron_idx]

        tri_i = make_unique((j, k, m))
        tri_j = make_unique((i, m, k))
        tri_k = make_unique((i, j, m))
        tri_m = make_unique((j, i, k))

        if neighbors[0] is None:
            surface_triangles.append(tri_i)

        if neighbors[1] is None:
            surface_triangles.append(tri_j)

        if neighbors[2] is None:
            surface_triangles.append(tri_k)

        if neighbors[3] is None:
            surface_triangles.append(tri_m)

        tetrahedron_idx += 1

    # Third we build a data structure to help us re-map vertex indices
    vertex_count = 0
    indices = dict()
    for triangle in surface_triangles:
        i = triangle[0]
        j = triangle[1]
        k = triangle[2]
        if i not in indices:
            indices[i] = vertex_count
            vertex_count += 1
        if j not in indices:
            indices[j] = vertex_count
            vertex_count += 1
        if k not in indices:
            indices[k] = vertex_count
            vertex_count += 1

    # Fourth we re-map triangles into new index space
    mesh = T3Mesh(t4mesh.name)
    for triangle in surface_triangles:
        i = indices[triangle[0]]
        j = indices[triangle[1]]
        k = indices[triangle[2]]
        mesh.triangles.append(make_unique((i, j, k)))

    # Fifth we fill in the vertex data
    mesh.vertices = [(sys.float_info.max, sys.float_info.max, sys.float_info.max) for _ in xrange(vertex_count)]
    for old_idx, new_idx in indices.iteritems():
        mesh.vertices[new_idx] = t4mesh.vertices[old_idx]

    return mesh


def make_box_mesh(name, width, height, depth):
    """
    This function creates a surface mesh of a box shape.

    :param name:      The name of the box
    :param width:
    :param height:
    :param depth:
    :return:            A triangle mesh of the surface of the given box shape mesh.
    """

    p0 = (-width/2.0, -height/2.0, -depth/2.0)
    p1 = (width/2.0, -height/2.0, -depth/2.0)
    p2 = (width/2.0, height/2.0, -depth/2.0)
    p3 = (-width/2.0, height/2.0, -depth/2.0)
    p4 = (-width/2.0, -height/2.0, depth/2.0)
    p5 = (width/2.0, -height/2.0, depth/2.0)
    p6 = (width/2.0, height/2.0, depth/2.0)
    p7 = (-width/2.0, height/2.0, depth/2.0)

    mesh = T3Mesh(name)
    mesh.vertices = [p0, p1, p2, p3, p4, p5, p6, p7]
    mesh.triangles.append(make_unique((0, 2, 1)))
    mesh.triangles.append(make_unique((0, 3, 2)))

    mesh.triangles.append(make_unique((1, 2, 6)))
    mesh.triangles.append(make_unique((1, 6, 5)))

    mesh.triangles.append(make_unique((0, 4, 3)))
    mesh.triangles.append(make_unique((3, 4, 7)))

    mesh.triangles.append(make_unique((0, 1, 5)))
    mesh.triangles.append(make_unique((0, 5, 4)))

    mesh.triangles.append(make_unique((3, 7, 6)))
    mesh.triangles.append(make_unique((2, 3, 6)))

    mesh.triangles.append(make_unique((4, 5, 6)))
    mesh.triangles.append(make_unique((4, 6, 7)))

    return mesh

