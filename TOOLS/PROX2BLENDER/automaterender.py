import bpy
import math
import mathutils
import os

#Clean empty scene
bpy.ops.wm.read_factory_settings(use_empty=True)

usdaFP = "finalScene.usda"
bpy.ops.wm.usd_import(filepath=usdaFP)

#Calculate the bounding box of visible scene geometry, ignoring objects far from origin
def getSceneBounds():
    minCoord = mathutils.Vector((float("inf"), float("inf"), float("inf")))
    maxCoord = mathutils.Vector((float("-inf"), float("-inf"), float("-inf")))
    
    # Set a threshold to ignore objects that are too far away
    distanceThreshold = 500  # Ignore objects further than this from origin
    
    for obj in bpy.context.scene.objects:
        #Skip objects that are too far from origin
        if obj.type == "MESH" and obj.matrix_world.translation.length < distanceThreshold:
            #Get world coordinates of the objects bounding box
            world_verts = [obj.matrix_world @ mathutils.Vector(v) for v in obj.bound_box]
            for vert in world_verts:
                if vert.length < distanceThreshold:
                    minCoord.x = min(minCoord.x, vert.x)
                    minCoord.y = min(minCoord.y, vert.y)
                    minCoord.z = min(minCoord.z, vert.z)
                    maxCoord.x = max(maxCoord.x, vert.x)
                    maxCoord.y = max(maxCoord.y, vert.y)
                    maxCoord.z = max(maxCoord.z, vert.z)
    
    #If no valid objects found, we just take some random coords
    if minCoord.x == float("inf"):
        minCoord = mathutils.Vector((-1, -1, -1))
        maxCoord = mathutils.Vector((1, 1, 1))
    
    return minCoord, maxCoord

#Get scene bounds
minCoord, maxCoord = getSceneBounds()
sceneCenter = (minCoord + maxCoord) / 2
sceneSize = maxCoord - minCoord
maxDimension = max(sceneSize.x, sceneSize.y, sceneSize.z)

if maxDimension <= 0.0:
    maxDimension = 1.0

#for light in bpy.data.lights:
#    bpy.data.lights.remove(light)



#lightData = bpy.data.lights.new(name="uniqueLightSun", type="POINT")
#lightData.energy = 1000


#lightObject = bpy.data.objects.new(name="sunLigh", object_data=lightData)

#lightObject.location = (0, 0, 10.0)
#bpy.context.collection.objects.link(lightObject)



lightData = bpy.data.lights.new(name="Sun", type="SUN")
lightData.energy = 1.0
lightData.angle = 0.526
lightObject = bpy.data.objects.new(name="Sun", object_data=lightData)
#Link light to the scene
bpy.context.collection.objects.link(lightObject)

lightObject.rotation_euler = (math.radians(45), math.radians(0), math.radians(30))

#lightObject.location = (50, 4, 107)



cameraDistance = maxDimension * 2
cameraHeight = maxDimension * 0.8

bpy.ops.object.camera_add(
    location=(
        sceneCenter.x + cameraDistance * math.cos(math.pi/4),
        sceneCenter.y - cameraDistance * math.sin(math.pi/4),
        sceneCenter.z + cameraHeight
    )
)

camera = bpy.context.active_object
bpy.context.scene.camera = camera

direction = sceneCenter - camera.location
rotQuat = direction.to_track_quat("-Z", "Y")
camera.rotation_euler = rotQuat.to_euler()

camera.data.clip_end = cameraDistance * 3

scene = bpy.context.scene
scene.render.fps = 100
#scene.render.engine = "BLENDER_EEVEE_NEXT"
scene.render.engine = "BLENDER_EEVEE_NEXT"
#scene.eevee.taa_samples = 128

eevee = scene.eevee

try:
    eevee.taa_samples = 64
    eevee.taa_render_samples = 1024
except Exception:
    pass


try:
    eevee.use_soft_shadows = True
    eevee.shadow_cube_size = 4096
    eevee.shadow_cascade_size = 4096
except Exception:
    pass

try:
    eevee.use_gtao = True
    eevee.gtao_distance = maxDimension * 0.25
    eevee.gtao_factor = 1.0
except Exception:
    pass

try:
    eevee.use_ssr = True
    eevee.use_ssr_refraction = True
    try:
        eevee.ssr_quality = 0.9
        eevee.ssr_max_roughness = 0.8
    except Exception:
        pass
except Exception:
    pass


try:
    eevee.gi_cubemap_resolution = 512
    eevee.gi_irradiance_smoothing = 0.5
except Exception:
    pass

scene.render.resolution_x = 1920
scene.render.resolution_y = 1080

print(f"Scene dimensions: {sceneSize}")
print(f"Maximum dimension: {maxDimension}")

scene.render.image_settings.file_format = "PNG"
scene.render.image_settings.color_mode = "RGBA"

framesDirRel = "//frames"
framesDir = bpy.path.abspath(framesDirRel)
os.makedirs(framesDir, exist_ok=True)

startFrame = scene.frame_start
endFrame = scene.frame_end

startFrame = int(startFrame)
endFrame = int(endFrame)
if endFrame < startFrame:
    endFrame = startFrame

print(f"Rendering frames {startFrame} .. {endFrame} to {framesDirRel} ({framesDir})")

for frame in range(startFrame, endFrame + 1):
    filename = os.path.join(framesDir, f"frame_{frame:04d}.png")
    if os.path.exists(filename):
        print(f"[{frame}] Skipping (already exists): {filename}")
        continue

    scene.frame_set(frame)
    scene.render.filepath = filename
    bpy.ops.render.render(write_still=True)
    print(f"[{frame}] Rendered -> {filename}")

print("Image sequence render complete (or stopped).")
print("When all frames are present combine into MP4 with ffmpeg, e.g.:")
#print(f"ffmpeg -framerate {scene.render.fps} -i {os.path.join(framesDir, "frame_%04d.png")} -c:v libx264 -pix_fmt yuv420p output_animation.mp4")
print(f"ffmpeg -framerate 100 -i frame_%04d.png -s 1920x1080 -c:v libx264 -preset fast -crf 18 -pix_fmt yuv420p output_video.mp4")
