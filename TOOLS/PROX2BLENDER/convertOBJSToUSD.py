import bpy
import os

inputFolder = "out/"
outputFolder = "out/"

#We collect all obj files :)
objFiles = [f for f in os.listdir(inputFolder) if f.endswith(".obj")]

for objFile in objFiles:
    objFile_path = os.path.join(inputFolder, objFile)
    usdFile_path = os.path.join(outputFolder, objFile.replace(".obj", ".usd"))

    #Select all objects and delete them
    bpy.ops.object.select_all(action="SELECT")
    bpy.ops.object.delete(use_global=False)

    #Import OBJ
    bpy.ops.wm.obj_import(filepath=objFile_path)

    #Assuming only one object is imported, create a new UV map and unwrap
    obj = bpy.context.view_layer.objects.active  # Get the active object (should be the imported one)
    
    if obj.type == 'MESH':
        #Add a new UV map (if not already existing)
        if not obj.data.uv_layers:
            obj.data.uv_layers.new(name="UVMap")
        
        #Enter Edit mode for unwrapping
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_all(action='SELECT')
        bpy.ops.uv.unwrap(method='ANGLE_BASED', margin=0.001)
        bpy.ops.object.mode_set(mode='OBJECT')

    #Export to USD
    bpy.ops.wm.usd_export(filepath=usdFile_path)

    #Remove the original OBJ file
    os.remove(objFile_path)
    print(f"Converted {objFile_path} to {usdFile_path} and deleted the .obj file.")

print("All .obj files processed, UV unwrapped, and converted to .usd.")

