from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import math

wfile = open("../data/w_coords.txt", "r")
xfile = open("../data/x_coords.txt", "r")

wx = []
wy = []
xx = []
xy = []
w = []
x = []

for line in wfile:
    w.append(line)

for line in xfile:
    x.append(line)

for i in range(0, len(w), 1):
    wx.append(float(w[i].split()[0][0:8]))
    wy.append(float(w[i].split()[1][0:8]))

for i in range(0, len(wx), 1):
    if (math.fabs(wx[i]) > 1):
        wx[i] /= 100
    if (math.fabs(wy[i]) > 1):
        wy[i] /= 100

for i in range(0, len(x), 1):
    xx.append(float(x[i].split()[0][0:8]))
    xy.append(float(x[i].split()[1][0:8]))

for i in range(0, len(xx), 1):
    if (math.fabs(xx[i]) > 1):
        xx[i] /= 100
    if (math.fabs(xy[i]) > 1):
        xy[i] /= 100
        
axes = plt.gca()
axes.set_xlim([-1.0,1.0])
axes.set_ylim([-1.0,1.0])

plt.xlabel("$x_1$")
plt.ylabel("$x_2$")
plt.scatter(xx, xy, s=10, alpha=0.3)
#plt.scatter(xx, xy, s=10, alpha=1.0)
#plt.scatter(wx, wy, s=10, alpha=0.3)

#plt.scatter([0.5],[0.5], c="red")
#plt.scatter([0.0],[0.0], c="red")
#plt.scatter([-0.5],[0.5], c="red")
#plt.scatter([-0.5],[-0.5], c="red")
#plt.scatter([0.5],[-0.5], c="red")
#plt.scatter([0.8],[0.8], c="red")

plt.show()
