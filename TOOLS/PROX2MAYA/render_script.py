import sys
import maya.standalone
import maya.cmds as cmds
import maya.utils


if __name__ == '__main__':
    sys.path.append('/Applications/solidangle/mtoa/2017/scripts/')

    maya.standalone.initialize()

    cmds.loadPlugin('mtoa')

    cmds.file(force=True, new=True)
    cmds.file(rename='render_test.ma')
    cmds.file(force=True, type='mayaAscii', save=True)

    cmds.polyCube()
    cmds.polyCube()
    cmds.polyCube()
    cmds.polyCube()

    # Make a world group
    cmds.select(clear=True)
    cmds.select(cmds.ls('pCube*'), visible=True)
    cmds.group(name='World')

    # Create a camera
    cmds.delete(cmds.ls('camera*'))
    name = cmds.camera(worldCenterOfInterest=[0, 1, 0], worldUp=[0, 1, 0], position=[1, 0, 1])
    cmds.lookThru(name)

    # Create area light
    light_name = mtoa.utils.createLocator('aiAreaLight', asLight=True)
    cmds.setAttr(light_name[0] + '.aiTranslator', 'quad', type='string')
    cmds.setAttr(light_name[0] + '.intensity', 300)
    cmds.setAttr(light_name[0] + '.aiSamples', 2)
    cmds.setAttr(light_name[1] + '.translateX', 10)
    cmds.setAttr(light_name[1] + '.translateY', 10)
    cmds.setAttr(light_name[1] + '.translateZ', 10)

    cmds.file(force=True, type='mayaAscii', save=True)

    maya.standalone.uninitialize()
