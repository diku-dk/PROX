import xml.etree.ElementTree as ET
import numpy as np
import trimesh
from scipy.spatial import ConvexHull
from collections import Counter
import math
import sys
import os
from pxr import Usd, UsdGeom, Gf, Sdf

#From https://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToMatrix/index.htm
def quaternionToMatrix(qw, qx, qy, qz):
    w, x, y, z = qw, qx, qy, qz
    return np.array([[1-2*(y*y+z*z), 2*(x*y - z*w), 2*(z*x + y*w)],
                     [2*(x*y + z*w), 1-2*(x*x+z*z), 2*(y*z - x*w)],
                     [2*(x*z - y*w), 2*(y*z + x*w), 1-2*(x*x+y*y)]])

def parseTransform(elem):
    #If we have non-tetrahedral, return identity matrix, identity translation
    if elem is None:
        return np.eye(3), np.zeros(3)
    translation = np.array([float(elem.get("x",0)), float(elem.get("y",0)), float(elem.get("z",0))])
    rotationMatrix = quaternionToMatrix(float(elem.get("qs",1)), float(elem.get("qx",0)),
                                        float(elem.get("qy",0)), float(elem.get("qz",0)))
    return rotationMatrix, translation

def compose(rot_a, trans_a, rot_b, trans_b):
    rot = rot_a.dot(rot_b)
    trans = rot_a.dot(trans_b) + trans_a
    return rot, trans

#Convert tetrahedrals to displayable triangles (fast, duplicates allowed)
def tetsToTriangles(tets):
    faces = []
    for (a, b, c, d) in tets:
        faces.append((a, b, c))
        faces.append((a, b, d))
        faces.append((a, c, d))
        faces.append((b, c, d))
    return np.array(faces, dtype=int)

def sceneToObjects(xmlPath, outDir):
    tree = ET.parse(xmlPath)
    root = tree.getroot()

    #collect geometries
    geoms = {}
    for geom in root.findall(".//geometry"):
        name = geom.get("name")
        #Which transform do we use?
        geomRotation, geomTranslation = parseTransform(geom.find("transform"))
        conv = geom.find("convex")
        convexPoints = None
        convexRotation, convexTranslation = np.eye(3), np.zeros(3)
        if conv is not None:
            #We parse and collect all points in our convex.
            convexPoints = [(float(p.get("x")), float(p.get("y")), float(p.get("z"))) for p in conv.findall("point")]
            #Save transform
            cr, ct = parseTransform(conv.find("transform"))
            convexRotation, convexTranslation = cr, ct

        # --- sphere primitive (optional) ---
        sph = geom.find("sphere")
        sphereRadius = None
        sphereRotation, sphereTranslation = np.eye(3), np.zeros(3)
        if sph is not None:
            sphereRadius = float(sph.get("radius", 1.0))
            sr, st = parseTransform(sph.find("transform"))
            sphereRotation, sphereTranslation = sr, st

        # tetramesh (optional) — do NOT raise if missing (may be sphere/convex)
        tet = geom.find("tetramesh")
        verts, tets = None, None
        if tet is not None:
            vmap = {}
            for v in tet.findall("vertex"):
                idx = int(v.get("idx")); vmap[idx] = (float(v.get("x")),float(v.get("y")),float(v.get("z")))
            # build ordered list (assumes indices 0..N-1 present)
            verts = np.array([vmap[i] for i in range(len(vmap))], dtype=float)
            tets = np.array([[int(t.get(k)) for k in ("i","j","k","m")] for t in tet.findall("tetrahedron")], dtype=int)

        geoms[name] = dict(geomRotation=geomRotation, geomTranslation=geomTranslation,
                           convexPoints=np.array(convexPoints) if convexPoints else None,
                           convexRotation=convexRotation, convexTranslation=convexTranslation,
                           # include sphere data
                           sphereRadius=sphereRadius, sphereRotation=sphereRotation, sphereTranslation=sphereTranslation,
                           verts=verts, tets=tets)

    os.makedirs(outDir, exist_ok=True)

    #Finally we can write out objects!
    for obj in root.findall(".//object"):
        name = obj.get("name")
        shape = obj.find("shape"); ifshape = shape is not None
        ref = shape.get("ref") if ifshape else None
        state = obj.find("state")
        objRotation, objTranslation = parseTransform(state.find("transform") if state is not None else None)
        if ref not in geoms:
            print("(ERROR?) skipped: ", name); continue
        geometries = geoms[ref]

        #Convex transform (if any) with our body transform.
        tmpRotation, tmpTranslation = compose(geometries["geomRotation"], geometries["geomTranslation"],
                                              geometries["convexRotation"], geometries["convexTranslation"])
        totalRotation, totalTranslation = compose(objRotation, objTranslation, tmpRotation, tmpTranslation)
        totalRotation = np.eye(3)
        totalTranslation = np.zeros(3)
        # --- Decide which representation to use ---
        mesh = None

        # 1) tetramesh if present
        if geometries["tets"] is not None and geometries["verts"] is not None:
            faces = tetsToTriangles(geometries["tets"])
            vertsWorld = (totalRotation.dot(geometries["verts"].T).T + totalTranslation)
            mesh = trimesh.Trimesh(vertices=vertsWorld, faces=faces, process=False)

        # 2) sphere primitive when no tets are present
        elif geometries.get("sphereRadius", None) is not None:
            # compose geometry-level + sphere-local transforms
            tmpR_s, tmpT_s = compose(geometries["geomRotation"], geometries["geomTranslation"],
                                      geometries["sphereRotation"], geometries["sphereTranslation"])
            totalR_s, totalT_s = compose(objRotation, objTranslation, tmpR_s, tmpT_s)
            # create icosphere and apply transform
            mesh = trimesh.creation.icosphere(subdivisions=3, radius=geometries["sphereRadius"])
            H = np.eye(4)
            H[:3,:3] = totalR_s
            H[:3,3] = totalT_s
            mesh.apply_transform(H)

        # 3) convex points -> convex hull
        elif geometries["convexPoints"] is not None:
            ptsWorld = (totalRotation.dot(geometries["convexPoints"].T).T + totalTranslation)
            hull = ConvexHull(ptsWorld)
            faces = hull.simplices
            mesh = trimesh.Trimesh(vertices=ptsWorld, faces=faces, process=False)

        else:
            print(f"Warning: geometry '{ref}' for object '{name}' has no supported shape -> skipping.")
            continue

        #Unfortunately Blender uses Z as the up-axis and OpenGL uses Y. Thus we
        # have to rotate 90 deg in the x-axis to get a Y-to-Z conversion.
        rotationMatrix = np.array([
            [1, 0, 0, 0],
            [0, 0, 1, 0],
            [0, -1, 0, 0],
            [0, 0, 0, 1]
        ])

        mesh.apply_transform(rotationMatrix)
        mesh.export(os.path.join(outDir, f"{name}.obj"))

        print("wrote", name + ".obj")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convertXMLToMeshes.py inputScene.xml outDir")
        sys.exit(1)
    input_xml = sys.argv[1]
    outDir = sys.argv[2]
    sceneToObjects(input_xml, outDir)
