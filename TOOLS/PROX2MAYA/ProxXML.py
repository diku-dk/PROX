from ProxMesh import T4Mesh, make_surface_mesh, make_box_mesh
from ProxMotionChannel import MotionChannel
from ProxUtils import make_euler_angles
import xml.etree.ElementTree as ET


class RawData:
    def __init__(self):
        self.geometries = []
        self.channels = []
        self.bodies = []


def load_box_shape(geometry_xml_tag):
    name = geometry_xml_tag.attrib['name']
    box_tag = geometry_xml_tag.find('box')
    width = float(box_tag.attrib['width'])
    height = float(box_tag.attrib['height'])
    depth = float(box_tag.attrib['depth'])
    return name, width, height, depth


def load_tetrahedral_mesh(geometry_xml_tag):
    """
    This function creates a tetrahedron mesh corresponding to the xml data
    given by the xml <geometry> tag as defined by the PROX engine.

    :param geometry_xml_tag:
    :return:
    """
    poly_name = geometry_xml_tag.attrib['name']
    number_of_vertices = sum(1 for _ in geometry_xml_tag.iter('vertex'))
    number_of_tetrahedra = sum(1 for _ in geometry_xml_tag.iter('tetrahedron'))
    mesh = T4Mesh(poly_name, number_of_vertices, number_of_tetrahedra)
    for vertex in geometry_xml_tag.iter('vertex'):
        idx = int(vertex.attrib['idx'])
        x = float(vertex.attrib['x'])
        y = float(vertex.attrib['y'])
        z = float(vertex.attrib['z'])
        mesh.vertices[idx] = (x, y, z)
    tetrahedron_idx = 0
    for tetrahedron in geometry_xml_tag.iter('tetrahedron'):
        i = int(tetrahedron.attrib['i'])
        j = int(tetrahedron.attrib['j'])
        k = int(tetrahedron.attrib['k'])
        m = int(tetrahedron.attrib['m'])
        mesh.tetrahedra[tetrahedron_idx] = (i, j, k, m)
        tetrahedron_idx += 1
    return mesh


def load_motion_channel(channel_xml_tag, verbose=False):
    body_name = channel_xml_tag.attrib['name']
    channel = MotionChannel(body_name)
    for key_tag in channel_xml_tag.iter('key'):
        time = float(key_tag.attrib['time'])
        x = float(key_tag.attrib['x'])
        y = float(key_tag.attrib['y'])
        z = float(key_tag.attrib['z'])
        qs = float(key_tag.attrib['qs'])
        qx = float(key_tag.attrib['qx'])
        qy = float(key_tag.attrib['qy'])
        qz = float(key_tag.attrib['qz'])
        (rx, ry, rz) = make_euler_angles(qs, qx, qy, qz)
        channel.time.append(time)
        channel.tx.append(x)
        channel.ty.append(y)
        channel.tz.append(z)
        channel.rx.append(rx)
        channel.ry.append(ry)
        channel.rz.append(rz)
    #if verbose:
    #    print 'Done loading motion channel for', body_name
    return channel


def load_from_xml(xml_filename, data, verbose=False):
    if verbose:
        print 'Started to parse', xml_filename
        print '  '

    xml = ET.parse(xml_filename)
    root = xml.getroot()

    if verbose:
        print 'Done parsing', xml_filename
        print '  '

    geometry_specs = dict()

    if root.tag == 'physics':
        for geometry_xml_tag in root.iter('geometry'):
            if not geometry_xml_tag.attrib['type'] == 'collection':
                continue

            tetramesh_tag = geometry_xml_tag.find('tetramesh')
            box_tag = geometry_xml_tag.find('box')
            if tetramesh_tag is not None:
                t4mesh = load_tetrahedral_mesh(geometry_xml_tag)
                t3mesh = make_surface_mesh(t4mesh)
                #if verbose:
                #    print '  '
                #    print 'volume mesh info of ', t4mesh.name
                #    print '  #vertices   = ', len(t4mesh.vertices)
                #    print '  #tetrahera  = ', len(t4mesh.tetrahedra)
                #    print '  '
                #    print 'surface mesh info of ', t3mesh.name
                #    print '  #vertices   = ', len(t3mesh.vertices)
                #    print '  #triangles  = ', len(t3mesh.triangles)
                #    print '  '
                data.geometries.append(t3mesh)
                geometry_specs[t4mesh.name] = [len(t4mesh.vertices), len(t4mesh.tetrahedra)]
            if box_tag is not None:
                name, width, height, depth = load_box_shape(geometry_xml_tag)
                t3mesh = make_box_mesh(name, width, height, depth)
                data.geometries.append(t3mesh)
                geometry_specs[name] = [0, 0]

        total_verts = 0
        total_tets = 0

        for obj_xml_tag in root.iter('object'):
            geometry_ref = obj_xml_tag.find('shape')
            if not obj_xml_tag.attrib['type'] == 'rigid':
                continue
            poly_name = geometry_ref.attrib['ref']
            body_name = obj_xml_tag.attrib['name']
            data.bodies.append((poly_name, body_name))
            total_verts += geometry_specs[poly_name][0]
            total_tets += geometry_specs[poly_name][1]

        if verbose:
            print '  '
            print '#objects = ', len(data.bodies)
            print '#vertices = ', total_verts
            print '#tetrahedra = ', total_tets
            print '  '

    elif root.tag == 'channels':
        for channel_xml_tag in root.iter('channel'):
            channel = load_motion_channel(channel_xml_tag, verbose)
            data.channels.append(channel)
    else:
        raise RuntimeError('XML file format not recognized')
