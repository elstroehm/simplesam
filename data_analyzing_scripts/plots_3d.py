from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import math

xfile = open("../data/x_coords.txt", "r")

x = []
y = []
z = []

for line in xfile:
    x.append(float(line.split()[0]))
    y.append(float(line.split()[1]))
    z.append(float(line.split()[2]))

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

ax.set_xlim([-1.0,1.0])
ax.set_ylim([-1.0,1.0])
ax.set_zlim([-1.0,1.0])

ax.set_xlabel("x_1")
ax.set_ylabel("x_2")
ax.set_zlabel("x_3")


ax.scatter(x, y, z, s=10, alpha=0.3)

#plt.scatter([0.5],[0.5], c="red")
#plt.scatter([0.0],[0.0], c="red")
#plt.scatter([-0.5],[0.5], c="red")
#plt.scatter([-0.5],[-0.5], c="red")
#plt.scatter([0.5],[-0.5], c="red")
#plt.scatter([0.8],[0.8], c="red")

plt.show()
