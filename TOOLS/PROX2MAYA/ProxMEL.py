import maya.cmds as cmds


def maya_mel_make_instance_copy(source_name, target_name, verbose=False):
    maya_name = cmds.instance(source_name)
    cmds.select(maya_name, r=True)
    cmds.rename(target_name)
    if verbose:
        print 'Made an instance copy of', source_name, ' with new name', target_name


def maya_mel_hide(name, verbose=False):
    cmds.select(name, r=True)
    cmds.hide()


def maya_mel_make_poly_geometry(mesh, verbose=False):
    first_triangle = True
    used = [False for _ in mesh.vertices]
    for triangle in mesh.triangles:
        i = triangle[0]
        j = triangle[1]
        k = triangle[2]

        pi = mesh.vertices[i]
        pj = mesh.vertices[j]
        pk = mesh.vertices[k]

        pattern = int(used[k]) << 2 | int(used[j]) << 1 | int(used[i]) << 0

        if pattern == 0:
            if first_triangle:
                cmds.polyCreateFacet(constructionHistory=False, p=[pi, pj, pk], name=mesh.name)
                first_triangle = False
            else:
                cmds.polyAppendVertex(constructionHistory=False, a=[pi, pj, pk])
        elif pattern == 1:
            cmds.polyAppendVertex(constructionHistory=False, a=[i, pj, pk])
        elif pattern == 2:
            cmds.polyAppendVertex(constructionHistory=False, a=[pi, j, pk])
        elif pattern == 3:
            cmds.polyAppendVertex(constructionHistory=False, a=[i, j, pk])
        elif pattern == 4:
            cmds.polyAppendVertex(constructionHistory=False, a=[pi, pj, k])
        elif pattern == 5:
            cmds.polyAppendVertex(constructionHistory=False, a=[i, pj, k])
        elif pattern == 6:
            cmds.polyAppendVertex(constructionHistory=False, a=[pi, j, k])
        elif pattern == 7:
            cmds.polyAppendVertex(constructionHistory=False, a=[i, j, k])
        else:
            raise RuntimeError('Internal error, illegal pattern value detected')
        used[i] = True
        used[j] = True
        used[k] = True
    cmds.polySetToFaceNormal()
    if verbose:
        print 'Created maya poly mesh of', mesh.name


def maya_mel_create_keyframe(body_name, time, pos, euler):
    time_value = format(time, '.2f') + 'sec'
    cmds.move(pos[0], pos[1], pos[2], body_name, absolute=True)
    cmds.rotate(euler[0], euler[1], euler[2], body_name, absolute=True)
    cmds.setKeyframe(body_name, time=time_value)


def maya_mel_fix_euler_angles(body_name):
    cmds.filterCurve(body_name + '_rotateX', body_name + '_rotateY', body_name + '_rotateZ')


def maya_mel_make_motion_channel(channel, verbose=False):
    body_name = channel.name
    for i in xrange(len(channel.time)):
        time = channel.time[i]
        x = channel.tx[i]
        y = channel.ty[i]
        z = channel.tz[i]
        rx = channel.rx[i]
        ry = channel.ry[i]
        rz = channel.rz[i]
        pos = (x, y, z)
        euler = (rx, ry, rz)
        maya_mel_create_keyframe(body_name, time, pos, euler)
    maya_mel_fix_euler_angles(body_name)
    if verbose:
        print 'Created motion channel for', body_name


def maya_mel_create_scene(data, verbose=False):
    for geometry in data.geometries:
        maya_mel_make_poly_geometry(geometry, verbose)

    for body in data.bodies:
        maya_mel_make_instance_copy(body[0], body[1], verbose)

    for geometry in data.geometries:
        maya_mel_hide(geometry.name, verbose)

    for channel in data.channels:
        maya_mel_make_motion_channel(channel, verbose)
