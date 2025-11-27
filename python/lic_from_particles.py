import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import sys

dimensions = 2
points = []
velocities = []
with open(sys.argv[1], "r") as f:
    n = int(f.readline())
    for i in range(n):
        line = f.readline().split()
        dimensions = (len(line) - 1) // 2
        if i == 0:
            for j in range(dimensions):
                velocities.append([])
        p = []
        v = []
        for j in range(dimensions):
            p.append(float(line[1 + j]))
            velocities[j].append(float(line[1 + dimensions + j]))
        points.append(p)

grid_x, grid_y = np.mgrid[0:1:100j, 0:1:200j]
points = np.array(points)

values = []
for i in range(dimensions):
    values.append(np.array(velocities[i]))

grid_z0 = griddata(points, values[0], (grid_x, grid_y), method='nearest')
grid_z1 = griddata(points, values[0], (grid_x, grid_y), method='linear')
grid_z2 = griddata(points, values[0], (grid_x, grid_y), method='cubic')

plt.plot(points[:,0], points[:,1], 'k.', ms=1)
plt.title('Original')
plt.show()

