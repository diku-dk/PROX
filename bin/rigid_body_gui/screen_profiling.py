




import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['lines.linewidth'] = 0.75
mpl.rcParams['lines.markersize'] = 6

f_type = 'Times';
f_size = 20;


import numpy as np
import matplotlib.pyplot as plt
import warnings

def toNp(x):
    return np.asarray(x, dtype=float)


Ekin = toNp(Ekin)
Epot = toNp(Epot)
Emech = Epot + Ekin
frames = np.arange(1, Ekin.size + 1)


plt.figure(100)
e1, = plt.plot(frames, Ekin, 'r-', label='Kinetic')
e2, = plt.plot(frames, Epot, 'g-', label='Potential')
e3, = plt.plot(frames, Emech, 'b-', label='Mechanical')
plt.title('Kinetic and potential energy of total system', fontsize=f_size, fontname=f_type)
plt.legend(handles=[e1, e2, e3], fontsize=f_size)
plt.ylabel('Energy (Joules)', fontsize=f_size, fontname=f_type)
plt.xlabel('Frame', fontsize=f_size, fontname=f_type)
plt.grid(True, linestyle=':', linewidth=0.5)
plt.tight_layout()
plt.savefig('energy_py.eps', format='eps', bbox_inches='tight')
plt.savefig('energy_py.png', dpi=300, bbox_inches='tight')


plt.figure(200)
avg_pen = toNp(average_penetration)
plt.plot(np.abs(avg_pen) * 100)
plt.title('Penetrations', fontsize=f_size, fontname=f_type)
plt.xlabel('Frame', fontsize=f_size, fontname=f_type)
plt.ylabel('Penetration depth (cm)', fontsize=f_size, fontname=f_type)
plt.grid(True, linestyle=':', linewidth=0.5)
plt.tight_layout()
plt.savefig('penetration_py.eps', format='eps', bbox_inches='tight')
plt.savefig('penetration_py.png', dpi=300, bbox_inches='tight')


plt.figure(300)
plt.plot(toNp(contacts))
plt.title('Contacts', fontsize=f_size, fontname=f_type)
plt.xlabel('Frame', fontsize=f_size, fontname=f_type)
plt.ylabel('Number of contacts', fontsize=f_size, fontname=f_type)
plt.grid(True, linestyle=':', linewidth=0.5)

plt.tight_layout()
plt.savefig('contacts_py.eps', format='eps', bbox_inches='tight')
plt.savefig('contacts_py.png', dpi=300, bbox_inches='tight')



plt.figure(400)
contactsNp = toNp(contacts)
solverNp = toNp(solver)

mask = contactsNp > 0
xVals = contactsNp[mask]

if solverNp.size == xVals.size:
    yVals = solverNp
elif solverNp.size == contactsNp.size:
    yVals = solverNp[mask]
else:
    minLen = min(solverNp.size, xVals.size)
    xVals = xVals[:minLen]
    yVals = solverNp[:minLen]

plt.plot(xVals, yVals, '.', markersize=6)
plt.title('Time used by solver vs number of contacts', fontsize=f_size, fontname=f_type)
plt.xlabel('Number of contacts', fontsize=f_size, fontname=f_type)
plt.ylabel('Solver (ms)', fontsize=f_size, fontname=f_type)
plt.grid(True, linestyle=':', linewidth=0.5)
plt.tight_layout()

plt.savefig('solver_time_vs_contacts_py.eps', format='eps', bbox_inches='tight')
plt.savefig('solver_time_vs_contacts_py.png', dpi=300, bbox_inches='tight')



plt.figure(500)
if solverNp.size == frames.size:
    xVal = frames[mask]
    yVal = solverNp[mask]
elif solverNp.size == mask.sum():
    xVal = frames[mask]
    yVal = solverNp
else:
    minLen = min(frames[mask].size, solverNp.size)
    print(f'warning: solver length ({solverNp.size}) not matching frames(mask) length ({frames[mask].size}). Truncating to {minLen}.')
    xVal = frames[mask][:minLen]
    yVal = solverNp[:minLen]

plt.plot(xVal, yVal)
plt.title('Time used by solver per frame', fontsize=f_size, fontname=f_type)
plt.xlabel('Frame', fontsize=f_size, fontname=f_type)
plt.ylabel('Solver (ms)', fontsize=f_size, fontname=f_type)
plt.grid(True, linestyle=':', linewidth=0.5)
plt.tight_layout()
plt.savefig('solver_time_per_frame_py.eps', format='eps', bbox_inches='tight')
plt.savefig('solver_time_per_frame_py.png', dpi=300, bbox_inches='tight')



plt.figure(600)
for i, ylist in enumerate(convergence, start=1):
    y = np.asarray(ylist, dtype=float)


    bad = ~np.isfinite(y) | (y <= 0)
    if np.any(bad):
        yClean = y.copy()
        yClean[bad] = np.nan
    else:
        yClean = y

    if yClean.size == 0:
        continue

    x = np.arange(0, yClean.size + 0)
    plt.semilogy(x, yClean, alpha=0.8)

plt.title('Convergence', fontsize=f_size, fontname=f_type)
plt.xlabel('Solver iteration', fontsize=f_size, fontname=f_type)
plt.ylabel('Natural merit function', fontsize=f_size, fontname=f_type)
plt.grid(True, which='both', linestyle=':', linewidth=0.5)

plt.savefig('convergence_py.eps', format='eps', bbox_inches='tight')
plt.savefig('convergence_py.png', dpi=300, bbox_inches='tight')

plt.show()


plt.figure(700)
for i, arr in enumerate(rfactor, start=1):
    y = toNp(arr)
    if y.size == 0:
        continue
    x = np.arange(0, y.size)
    plt.plot(x, y, alpha=0.8)

plt.title('R-factor development', fontsize=f_size, fontname=f_type)
plt.xlabel('Solver iteration', fontsize=f_size, fontname=f_type)
plt.grid(True, linestyle=':', linewidth=0.5)
plt.tight_layout()
plt.savefig('rfactors_py.png', dpi=300, bbox_inches='tight')
plt.show()

maxLen = max((len(c) for c in convergence), default=0)
if maxLen == 0:
    raise ValueError('all convergence entries are empty')

data = np.full((len(convergence), maxLen), np.nan, dtype=float)
for i, ylist in enumerate(convergence):
    y = np.asarray(ylist, dtype=float)
    if y.size == 0:
        continue
    data[i, :y.size] = y
    bad = ~np.isfinite(data[i, :y.size]) | (data[i, :y.size] <= 0)
    if np.any(bad):
        data[i, :y.size][bad] = np.nan


avgCurve = np.nanmean(data, axis=0)

x = np.arange(0, avgCurve.size)
plt.figure(800)
plt.semilogy(x, avgCurve, linewidth=0.75)
plt.title('Average Convergence', fontsize=f_size, fontname=f_type)
plt.xlabel('Solver iteration', fontsize=f_size, fontname=f_type)
plt.ylabel('Natural merit function', fontsize=f_size, fontname=f_type)
plt.grid(True, which='both', linestyle=':', linewidth=0.5)

plt.savefig('convergence_avg.eps', format='eps', bbox_inches='tight')
plt.savefig('convergence_avg.png', dpi=300, bbox_inches='tight')
plt.show()


convergenceData = np.array(convergence)
iterationsData = convergenceData.T

q1 = np.percentile(iterationsData, 25, axis=1)
median = np.percentile(iterationsData, 50, axis=1)
q3 = np.percentile(iterationsData, 75, axis=1)

plt.figure(figsize=(10, 6))
plt.plot(q1, label='Q1 (25th Percentile)', color='blue', linestyle='--')
plt.plot(median, label='Median (50th Percentile)', color='green', linestyle='-')
plt.plot(q3, label='Q3 (75th Percentile)', color='red', linestyle='--')

plt.title('Quartile Plot for Convergence Across Multiple Frames', fontsize=f_size, fontname=f_type)
plt.xlabel('Solver iteration', fontsize=f_size, fontname=f_type)
plt.ylabel('Natural merit function', fontsize=f_size, fontname=f_type)
plt.yscale('log')
plt.legend()

plt.savefig('convergence_quartile.eps', format='eps', bbox_inches='tight')
plt.savefig('convergence_quartile.png', dpi=300, bbox_inches='tight')
plt.show()
