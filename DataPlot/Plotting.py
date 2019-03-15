import os
import sys
import subprocess
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
#-----------------------------#
xSize = float(sys.argv[1])
ySize = float(sys.argv[2])

xMin = float(sys.argv[3])
xMax = float(sys.argv[4])
xStep = (xMax - xMin) / xSize

yMin = float(sys.argv[5])
yMax = float(sys.argv[6])
yStep = (yMax - yMin) / ySize

X = np.arange(xMin, xMax, xStep)
Y = np.arange(yMin, yMax, yStep)
X, Y = np.meshgrid(X, Y)
#-----------------------------#
DataPath = "./Data/"

for FileName in os.listdir(DataPath):
    Z = []

    for Line in open(DataPath + FileName, "r"):
        Z.append([float(Number) for Number in Line.split()])

    Z = np.asarray(Z)

    for i in range(0, int(len(Z) / ySize)):
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.set_zlim(np.amin(Z), np.amax(Z))
        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

        surf = ax.plot_surface(X, Y, Z[int(i * ySize):int((i + 1) * ySize)], cmap=cm.coolwarm, linewidth=0, antialiased=False, vmax=9.999, vmin=10.001)
        fig.colorbar(surf, shrink=0.5, aspect=5)

        print(i)

        plt.savefig("Images/" + "Image" + "_%03d.png" % i)
        surf.remove()
        plt.close(fig)

    subprocess.call(["./MakeVideo", FileName[:-4]])

    for file_name in glob.glob("Images/*.png"):
        os.remove(file_name)

