import glob
import os
import runpy
import re
import sys
import math
from pxr import Usd, UsdGeom, UsdShade, Gf, Sdf, Vt

outUSDA = "finalScene.usda"
sceneNames = "sliding_point_rigidBodiesData_*.py"
meshUsdDir = "out"
gizmoMeshPathA = "assets/gizmoBlue.usdc"
gizmoMeshPathB = "assets/gizmoRed.usdc"
FPS = 100

#Set the below to true if you want to apply materials 
gApplyTextures = True
#Set the below to true if you wish to apply no materials but set a custom colour
# below overrides above
gUseSingleCoolTextures = True 

textureDir = "textures"
textureLookup = {
    "Ground": "groundTexture.png",
    "Stone":  "stoneTexture.png",
    "Cannonball": "cannonballTexture.png"
}
g_textureColLookup = {
    "Ground": Gf.Vec3f(1.0, 1.0, 1.0),
    "Cannonball": Gf.Vec3f(0.1, 0.1, 0.1),
    "Stone":  Gf.Vec3f(0.662745098, 0.662745098, 0.662745098),
}

fallbackColour = Gf.Vec3f(0.6, 0.0, 0.8)

#If false, B-gizmos are not created
gRenderGizmosB = False

#Distance epsilon for deduplication
gGizmoDedupEps = 0.2
#Normal epsilon for deduplication
gGizmoNormalEps = 0.2

def _squared_dist(a, b):
    return (a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2

#Roughly o(n log n) -- I use sorting, should be ok even in python.
def dedupePositionsAndNormals(positions, normals, epsPos, epsNorm):
    if not positions:
        return [], [] if normals else None

    epsPos2 = epsPos * epsPos
    epsNorm2 = epsNorm * epsNorm if normals else None

    #pack (pos, norm, idx) so we can sort by position
    items = []
    for i, p in enumerate(positions):
        n = normals[i] if normals else None
        items.append((p, n, i))
    items.sort(key=lambda x: (round(x[0][0]/epsPos), 
                              round(x[0][1]/epsPos), 
                              round(x[0][2]/epsPos)))

    uniqPos, uniqNorm = [], [] if normals else None

    def close(a, b, thresh2):
        dx = a[0] - b[0]; dy = a[1] - b[1]; dz = a[2] - b[2]
        return dx*dx + dy*dy + dz*dz <= thresh2

    for p, n, _ in items:
        if uniqPos:
            up = uniqPos[-1]
            if close(p, up, epsPos2):
                if normals:
                    if close(n, uniqNorm[-1], epsNorm2):
                        #dupe
                        continue
                else:
                    #dupe
                    continue
        uniqPos.append((float(p[0]), float(p[1]), float(p[2])))
        if normals:
            uniqNorm.append((float(n[0]), float(n[1]), float(n[2])))

    return uniqPos, uniqNorm

def loadFrames(sceneGlob):
    files = sorted(glob.glob(sceneGlob))
    if not files:
        raise RuntimeError("No [scene]_rigidBodiesData_*.py files found.")
    #files = files[0:5]
    frames = []
    for fname in files:
        m = re.search(r'_(\d+)\.py$', fname)
        if not m:
            raise RuntimeError(f"Filename {fname} doesn't match pattern scene_####.py")
        idx = int(m.group(1))
        ns = runpy.run_path(fname)

        posKey = f"positions_{idx}"
        oriKey = f"orientations_{idx}"
        namesKey = f"rigidNames_{idx}"
        if posKey not in ns or oriKey not in ns or namesKey not in ns:
            raise RuntimeError(f"{fname} missing {posKey}, {oriKey}, or {namesKey}")

        cxKey = f"CX_{idx}"
        cyKey = f"CY_{idx}"
        czKey = f"CZ_{idx}"
        nxKey = f"NX_{idx}"
        nyKey = f"NY_{idx}"
        nzKey = f"NZ_{idx}"

        contactPoints = None
        normals = None
        if all(key in ns for key in [cxKey, cyKey, czKey]):
            contactPoints = list(zip(ns[cxKey], ns[cyKey], ns[czKey]))
        if all(key in ns for key in [nxKey, nyKey, nzKey]):
            normals = list(zip(ns[nxKey], ns[nyKey], ns[nzKey]))


        bxKey = f"BX_{idx}"
        byKey = f"BY_{idx}"
        bzKey = f"BZ_{idx}"
        axKey = f"AX_{idx}"
        ayKey = f"AY_{idx}"
        azKey = f"AZ_{idx}"

        print(f"\rProcessing {idx} frames!", end="")
        
        bPoints = None
        bNormals = None
        if all(key in ns for key in [bxKey, byKey, bzKey]):
            bPoints = list(zip(ns[bxKey], ns[byKey], ns[bzKey]))
        if all(key in ns for key in [axKey, ayKey, azKey]):
            bNormals = list(zip(ns[axKey], ns[ayKey], ns[azKey]))

        contactPointsDeduped, normals_deduped = dedupePositionsAndNormals(
            contactPoints, normals, gGizmoDedupEps, gGizmoNormalEps
        )
        """bPointsDeduped, bNormalsDeduped = dedupePositionsAndNormals(
            bPoints, bNormals, gGizmoDedupEps, gGizmoNormalEps
        )"""
        

        frames.append({
            "index": idx,
            "positions": ns[posKey],
            "orientations": ns[oriKey],
            "rigidNames": ns[namesKey],
            "contactPoints": contactPointsDeduped,
            "normals": normals_deduped,
            "bPoints": bPoints,
            "bNormals": bNormals,
        })
    return frames

def getRigidNames(firstSceneFile="sliding_point_rigidBodiesData_0000.py"):
    ns = runpy.run_path(firstSceneFile)
    for key in ns:
        if key.startswith("rigidNames_"):
            return ns[key]
    raise RuntimeError("Couldn't find rigidNames_* in first frame file.")

def rotationQuatFromNormal(normalTuple):
    forward = Gf.Vec3d(0.0, 0.0, 1.0)
    nx, ny, nz = float(normalTuple[0]), float(normalTuple[1]), float(normalTuple[2])
    normalVec = Gf.Vec3d(nx, ny, nz)
    if normalVec.GetLength() == 0.0:
        return Gf.Quatf(1.0, 0.0, 0.0, 0.0)
    normalVec = normalVec.GetNormalized()
    rot = Gf.Rotation(forward, normalVec)
    q = rot.GetQuat()
    return Gf.Quatf(q.GetReal(), q.GetImaginary()[0], q.GetImaginary()[1], q.GetImaginary()[2])

def makeMaterialsForScene(stage, materialsList):
    materialCache = {}
    looksRootPath = "/Looks"
    UsdGeom.Xform.Define(stage, looksRootPath)

    uniqueNames = []
    for m in materialsList:
        if m not in uniqueNames:
            uniqueNames.append(m)

    print("[material] Preparing materials (gApplyTextures=%s), textureDir=%s" % (gApplyTextures, textureDir))
    for matName in uniqueNames:
        matPath = f"{looksRootPath}/{matName}"
        print(f"[material] Creating material '{matName}' at {matPath}")
        material = UsdShade.Material.Define(stage, matPath)

        pbrPath = f"{matPath}/PreviewSurface"
        pbr = UsdShade.Shader.Define(stage, pbrPath)
        pbr.CreateIdAttr("UsdPreviewSurface")
        pbr.CreateOutput("surface", Sdf.ValueTypeNames.Token)
        diffInput = pbr.CreateInput("diffuseColor", Sdf.ValueTypeNames.Color3f)

        texFile = None
        if gApplyTextures:
            tex_fname = textureLookup.get(matName)
            if tex_fname:
                candidate = os.path.join(textureDir, tex_fname)
                if os.path.exists(candidate):
                    texFile = candidate
                    print(f"[material] Found texture for '{matName}': {candidate}")
                else:
                    print(f"[material] Texture for '{matName}' expected at {candidate} but not found -> will use fallback color")

        if texFile and not gUseSingleCoolTextures:
            texPath = f"{matPath}/Texture"
            texShader = UsdShade.Shader.Define(stage, texPath)
            texShader.CreateIdAttr("UsdUVTexture")
            texShader.CreateInput("file", Sdf.ValueTypeNames.Asset).Set(Sdf.AssetPath(texFile))
            texShader.CreateOutput("rgb", Sdf.ValueTypeNames.Color3f)

            primvarPath = f"{matPath}/PrimvarReader_st"
            primvarReader = UsdShade.Shader.Define(stage, primvarPath)
            primvarReader.CreateIdAttr("UsdPrimvarReader_float2")
            primvarReader.CreateInput("varname", Sdf.ValueTypeNames.Token).Set("st")
            primvarReader.CreateOutput("result", Sdf.ValueTypeNames.Float2)

            texInputSt = texShader.CreateInput("st", Sdf.ValueTypeNames.Float2)
            texInputSt.ConnectToSource(primvarReader.GetOutput("result"))

            diffInput.ConnectToSource(texShader.GetOutput("rgb"))
        else:
            if gUseSingleCoolTextures:
                someMatF = g_textureColLookup.get(matName);
                if (someMatF):
                    diffInput.Set(someMatF)
                else:
                    diffInput.Set(fallbackColour)
                    print(f"[material] Using fallback color for '{matName}'")
            else:
                diffInput.Set(fallbackColour)
                if gApplyTextures:
                    print(f"[material] Using fallback color for '{matName}'")

        surfaceOut = material.CreateSurfaceOutput()
        surfaceOut.ConnectToSource(pbr.GetOutput("surface"))

        materialCache[matName] = material

    return materialCache

def buildUsd(frames, rigidNames, outUsd, meshUsdDir, gizmoMeshPath, gizmoMeshPathB, fps=100):
    stage = Usd.Stage.CreateNew(outUsd)
    stage.SetMetadata("metersPerUnit", 1.0)
    stage.SetMetadata("upAxis", "Y")
    start = 0.0
    end = float(len(frames) - 1)
    stage.SetStartTimeCode(start)
    stage.SetEndTimeCode(end)
    stage.SetTimeCodesPerSecond(float(fps))
    stage.SetFramesPerSecond(float(fps))

    world = stage.DefinePrim("/World", "Xform")
    stage.SetDefaultPrim(world)

    gizmosParent = stage.DefinePrim("/World/Gizmos", "Xform")
    gizmosBParent = None
    if gRenderGizmosB:
        gizmosBParent = stage.DefinePrim("/World/GizmosB", "Xform")

    firstSceneFiles = sorted(glob.glob(sceneNames))
    materials = None
    if firstSceneFiles:
        ns0 = runpy.run_path(firstSceneFiles[0])
        for key, val in ns0.items():
            if key.startswith("materials_"):
                materials = val
                print(f"[buildUsd] Found materials variable '{key}' in {firstSceneFiles[0]}")
                break
    if materials is None:
        materials = ["Default"] * len(rigidNames)
        print("[buildUsd] No materials_* found; using defaults")

    materialsMap = makeMaterialsForScene(stage, materials)

    for i, name in enumerate(rigidNames):
        primPath = f"/World/{name}"
        prim = stage.DefinePrim(primPath, "Xform")

        meshUsdFname = os.path.join(meshUsdDir, f"{name}.usd")
        if not os.path.exists(meshUsdFname):
            meshUsdFname2 = os.path.join(meshUsdDir, f"{name}.usda")
            if os.path.exists(meshUsdFname2):
                meshUsdFname = meshUsdFname2
            else:
                raise RuntimeError(f"Missing converted mesh USD for {name}: tried {meshUsdFname}")

        prim.GetReferences().AddReference(meshUsdFname)

        xformable = UsdGeom.Xformable(prim)
        translateOp = xformable.AddTranslateOp()
        orientOp = xformable.AddOrientOp()

        if gApplyTextures:
            matName = materials[i] if i < len(materials) else "Default"
            mat = materialsMap.get(matName)
            if mat:
                print(f"[buildUsd] Binding material '{matName}' -> prim {primPath}")
                UsdShade.MaterialBindingAPI(prim).Bind(mat)

        for frameNumber, frame in enumerate(frames):
            pos = frame["positions"][i]
            ori = frame["orientations"][i]
            #x, y, z, w = ori[0], ori[1], ori[2], ori[3]
            w, x, y, z = ori[0], ori[1], ori[2], ori[3]
            usdQuat = Gf.Quatf(w, x, y, z)
            tcode = Usd.TimeCode(frameNumber)
            translateOp.Set(Gf.Vec3d(pos[0], pos[1], pos[2]), tcode)
            orientOp.Set(usdQuat, tcode)

    # --- Contact-point gizmos with PointInstancer ---
    gizmoSomeRoot = "/World/Gizmos/Prototypes"
    UsdGeom.Xform.Define(stage, gizmoSomeRoot)

    protoApath = f"{gizmoSomeRoot}/GizmoProtoA"
    protoA = stage.DefinePrim(protoApath, "Xform")
    protoA.GetReferences().AddReference(gizmoMeshPath)
    protoA.SetInstanceable(True)

    protoB_path = f"{gizmoSomeRoot}/GizmoProtoB"
    if gRenderGizmosB:
        protoB = stage.DefinePrim(protoB_path, "Xform")
        protoB.GetReferences().AddReference(gizmoMeshPathB)
        protoB.SetInstanceable(True)

    #Create a PointInstancer for contact gizmos
    instancer_path = "/World/Gizmos/Instancer"
    inst_prim = stage.DefinePrim(instancer_path, "PointInstancer")
    inst = UsdGeom.PointInstancer(inst_prim)

    protosRel = inst.CreatePrototypesRel()
    protosRel.AddTarget(Sdf.Path(protoApath))
    if gRenderGizmosB:
        protosRel.AddTarget(Sdf.Path(protoB_path))


    inst.CreatePositionsAttr()
    inst.CreateOrientationsAttr()
    inst.CreateScalesAttr()
    inst.CreateProtoIndicesAttr()
    inst.CreateIdsAttr()

    for frameNumber, frame in enumerate(frames):
        print(f"\rProcessing {frameNumber} frame (out of {len(frames)})!", end="")
        tcode = Usd.TimeCode(frameNumber)
        pts = frame.get("contactPoints") or []
        norms = frame.get("normals") or []

        if pts:
            posArray = Vt.Vec3fArray([Gf.Vec3f(*p) for p in pts])
        else:
            posArray = Vt.Vec3fArray()

        orientsFloat = []
        for i, p in enumerate(pts):
            qf = Gf.Quatf(1.0, 0.0, 0.0, 0.0)
            if norms and i < len(norms):
                n = norms[i]
                if n != (0, 0, 0):
                    qf = rotationQuatFromNormal(n)
            orientsFloat.append(qf)
        #convert to Gf.Quath and Vt.QuathArray (USD expects Quath for PointInstancer.orientations)
        if orientsFloat:
            orientsQuath = [Gf.Quath(q.GetReal(),
                                      q.GetImaginary()[0],
                                      q.GetImaginary()[1],
                                      q.GetImaginary()[2]) for q in orientsFloat]
            orientArray = Vt.QuathArray(orientsQuath)
        else:
            orientArray = Vt.QuathArray()

        #We scale down, o.w. it is too large
        scalesArray = Vt.Vec3fArray([Gf.Vec3f(0.1, 0.1, 0.1) for _ in pts]) if pts else Vt.Vec3fArray()

        #protoIndices: choose prototype A (0) for all contact gizmos here
        protoIndices = Vt.IntArray([0] * len(pts)) if pts else Vt.IntArray()

        ids = Vt.Int64Array([frameNumber * 100000 + i for i in range(len(pts))]) if pts else Vt.Int64Array()

        inst.GetPositionsAttr().Set(posArray, tcode)
        inst.GetOrientationsAttr().Set(orientArray, tcode)
        inst.GetScalesAttr().Set(scalesArray, tcode)
        inst.GetProtoIndicesAttr().Set(protoIndices, tcode)
        inst.GetIdsAttr().Set(ids, tcode)


    if gRenderGizmosB:
        protoRootB = "/World/GizmosB/Prototypes"
        UsdGeom.Xform.Define(stage, protoRootB)
        #reference prototype B under B prototypes (we already defined protoB above under /World/Gizmos/Prototypes
        #well reuse protoBPath)
        instancerB_path = "/World/GizmosB/Instancer"
        instBPrim = stage.DefinePrim(instancerB_path, "PointInstancer")
        instB = UsdGeom.PointInstancer(instBPrim)
        protos_rel_B = instB.CreatePrototypesRel()
        protos_rel_B.AddTarget(Sdf.Path(protoB_path))
        instB.CreatePositionsAttr()
        instB.CreateOrientationsAttr()
        instB.CreateScalesAttr()
        instB.CreateProtoIndicesAttr()
        instB.CreateIdsAttr()

        for frameNumber, frame in enumerate(frames):
            tcode = Usd.TimeCode(frameNumber)
            bpts = frame.get("bPoints") or []
            bn = frame.get("bNormals") or []

            posArray = Vt.Vec3fArray([Gf.Vec3f(*p) for p in bpts]) if bpts else Vt.Vec3fArray()

            orientsFloat = []
            for i, p in enumerate(bpts):
                qf = Gf.Quatf(1.0, 0.0, 0.0, 0.0)
                if bn and i < len(bn):
                    bnormal = bn[i]
                    if bnormal != (0,0,0):
                        qf = rotationQuatFromNormal(bnormal)
                orientsFloat.append(qf)
            if orientsFloat:
                orientsQuath = [Gf.Quath(q.GetReal(),
                                          q.GetImaginary()[0],
                                          q.GetImaginary()[1],
                                          q.GetImaginary()[2]) for q in orientsFloat]
                orientArray = Vt.QuathArray(orientsQuath)
            else:
                orientArray = Vt.QuathArray()

            scalesArray = Vt.Vec3fArray([Gf.Vec3f(0.1, 0.1, 0.1) for _ in bpts]) if bpts else Vt.Vec3fArray()
            protoIndices = Vt.IntArray([0] * len(bpts)) if bpts else Vt.IntArray()
            ids = Vt.IntArray([frameNumber * 100000 + i for i in range(len(bpts))]) if bpts else Vt.IntArray()

            instB.GetPositionsAttr().Set(posArray, tcode)
            instB.GetOrientationsAttr().Set(orientArray, tcode)
            instB.GetScalesAttr().Set(scalesArray, tcode)
            instB.GetProtoIndicesAttr().Set(protoIndices, tcode)
            instB.GetIdsAttr().Set(ids, tcode)

    stage.GetRootLayer().Save()
    print("FINISHED")

# -------------------------
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python buildUSDAnimation.py [scene]rigidBodiesData_*.py")
        sys.exit(1)
    sceneNames = sys.argv[1]
    frames = loadFrames(sceneNames)
    firstScene = sorted(glob.glob(sceneNames))[0]
    rigidNames = getRigidNames(firstScene)
    buildUsd(frames, rigidNames, outUSDA, meshUsdDir, gizmoMeshPathA, gizmoMeshPathB, FPS)
