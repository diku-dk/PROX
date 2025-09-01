import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  #noqa: F401

csvPath = sys.argv[1] if len(sys.argv) > 1 else "sdf.csv"
df = pd.read_csv(csvPath)
assert set(["x", "y", "z", "sdf"]).issubset(df.columns), "CSV must have columns x,y,z,sdf"

xs = df["x"].to_numpy()
ys = df["y"].to_numpy()
zs = df["z"].to_numpy()
sdfVals = df["sdf"].to_numpy()

#infer grid resolution from unique coordinates
uniqX = np.unique(xs)
uniqY = np.unique(ys)
uniqZ = np.unique(zs)
nx, ny, nz = len(uniqX), len(uniqY), len(uniqZ)
total = nx * ny * nz
if total != len(df):
    print(f"Warning: unique grid sizes {nx}x{ny}x{nz} != {len(df)} rows. Proceeding but results may be wrong.")
print(f"Detected grid: nx={nx}, ny={ny}, nz={nz} (total {total})")

#reshape SDF into (nz, ny, nx) to match (k,j,i) loop order
try:
    vol = sdfVals.reshape((nz, ny, nx))
except Exception:
    print("Could not reshape directly; attempting to reorder by (z,y,x) sorted order.")
    dfSorted = df.sort_values(by=["z", "y", "x"]).reset_index(drop=True)
    vol = dfSorted["sdf"].to_numpy().reshape((nz, ny, nx))

#scatter plot of SDF points
fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(121, projection="3d")
maxPoints = 280000
idxs = np.arange(len(xs))
if len(xs) > maxPoints:
    rng = np.random.default_rng(1234)
    idxs = rng.choice(idxs, size=maxPoints, replace=False)
sc = ax.scatter(xs[idxs], ys[idxs], zs[idxs], c=sdfVals[idxs], cmap="coolwarm", s=6, marker="o")
ax.set_title("Node point cloud (colored by SDF)")
ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
cbar = plt.colorbar(sc, ax=ax, shrink=0.6)
cbar.set_label("sdf")

#second view: only negative or zero SDF points
ax2 = fig.add_subplot(122, projection="3d")
mask = sdfVals <= 0.0
negXs, negYs, negZs, negSdf = xs[mask], ys[mask], zs[mask], sdfVals[mask]
sc2 = ax2.scatter(negXs, negYs, negZs, c=negSdf, cmap="coolwarm", s=6, marker="o")
ax2.set_title("Negative / Zero SDF points")
ax2.set_xlabel("x"); ax2.set_ylabel("y"); ax2.set_zlabel("z")
cbar2 = plt.colorbar(sc2, ax=ax2, shrink=0.6)
cbar2.set_label("sdf")

plt.tight_layout()
plt.show()
