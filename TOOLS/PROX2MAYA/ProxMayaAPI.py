from ProxMEL import maya_mel_hide, maya_mel_fix_euler_angles
import maya.cmds as cmds
import maya.api.OpenMaya       as OpenMaya        # Common classes
import maya.api.OpenMayaAnim   as OpenMayaAnim    # Animation classes


def maya_api_addkeys(plugName, times, values):
    sel = OpenMaya.MSelectionList()
    sel.add(plugName)
    plug = sel.getPlug(0)

    animFn = OpenMayaAnim.MFnAnimCurve()
    animFn.create( plug, OpenMayaAnim.MFnAnimCurve.kAnimCurveTL)

    timeArray = OpenMaya.MTimeArray()
    valueArray = OpenMaya.MDoubleArray()

    for i in range(len(times)):
        timeArray.append(OpenMaya.MTime(times[i], OpenMaya.MTime.kSeconds))
        valueArray.append(values[i])

    animFn.addKeys(timeArray, valueArray)


def maya_api_create_mesh(mesh):
    vertices = [OpenMaya.MPoint(v[0], v[1], v[2]) for v in mesh.vertices]

    vertices_per_polygon = [3 for _ in mesh.triangles]

    polygon_indices = []
    for t in mesh.triangles:
        polygon_indices.append(t[0])
        polygon_indices.append(t[1])
        polygon_indices.append(t[2])

    meshFn = OpenMaya.MFnMesh()

    parentNode = meshFn.create(vertices, vertices_per_polygon, polygon_indices)

    meshFn.setName(mesh.name + 'Shape')

    transFn = OpenMaya.MFnTransform()
    transFn.setObject(parentNode)
    transFn.setName(mesh.name)


def maya_api_make_instance_copy( from_name, to_name):
    sel = OpenMaya.MSelectionList()
    sel.add(from_name)
    depNode = sel.getDependNode(0)
    dagPath = sel.getDagPath(0)

    depFn = OpenMaya.MFnDependencyNode()
    dagFn = OpenMaya.MFnDagNode()

    depFn.setObject(depNode)
    dagFn.setObject(dagPath)

    node = dagFn.duplicate(instance=True, instanceLeaf=False)

    #new_dagPath = OpenMaya.MDagPath.getAPathTo(node)
    #node = new_dagPath.transform()
    transFn = OpenMaya.MFnTransform()
    transFn.setObject(node)
    transFn.setName(to_name)


def maya_api_create_scene(data):
    for geometry in data.geometries:
        maya_api_create_mesh(geometry)

    for body in data.bodies:
        maya_api_make_instance_copy(body[0], body[1])

    for geometry in data.geometries:
        maya_mel_hide(geometry.name)

    for channel in data.channels:
        maya_api_addkeys(channel.name + '.tx', channel.time, channel.tx)
        maya_api_addkeys(channel.name + '.ty', channel.time, channel.ty)
        maya_api_addkeys(channel.name + '.tz', channel.time, channel.tz)
        maya_api_addkeys(channel.name + '.rx', channel.time, channel.rx)
        maya_api_addkeys(channel.name + '.ry', channel.time, channel.ry)
        maya_api_addkeys(channel.name + '.rz', channel.time, channel.rz)
        maya_mel_fix_euler_angles(channel.name)

    cmds.select(clear=True)
    cmds.select(cmds.ls('body*'), visible=True)
    cmds.polySetToFaceNormal()
    cmds.group(name='World')
    cmds.sets(forceElement='initialShadingGroup')