import os
#import maya.standalone
#import maya.cmds as cmds
#import maya.api.OpenMaya       as OpenMaya        # Common classes
#import maya.api.OpenMayaAnim   as OpenMayaAnim    # Animation classes

def import_multiple_object_files(path, prefix):
    files_to_import = []
    for file in os.listdir(path):
        if file.startswith(prefix) and file.endswith('.obj'):
            files_to_import.append(file)
            print(file)

    for file in files_to_import:
        absolute_path = path + file
        #cmds.file('%s' % absolute_path, i=True, type='OBJ', ignoreVersion=True, options='mo=1')
        print(absolute_path)

if __name__ == '__main__':

#    maya.standalone.initialize()
#    cmds.file(force=True, new=True)
#    cmds.file(rename='/Users/kennyerleben/Documents/PROX_TRUNK/TOOLS/PROX2MAYA/bin/cantilever_000001.ma')
#    cmds.file(force=True, type='mayaAscii', save=True)
#    import_multiple_object_files(  '000001' );
#    cmds.file(force=True, type='mayaAscii', save=True)
    import_multiple_object_files('./bin/cantilever/','000001')
#    maya.standalone.uninitialize()
