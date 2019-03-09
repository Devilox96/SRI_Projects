from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import subprocess
import glob
import os
#-----------------------------#
X = np.arange(-1, 1, 0.01)
Y = np.arange(-1, 1, 0.01)
X, Y = np.meshgrid(X, Y)
#-----------------------------#
for i in range(0, 300):
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    ax.set_zlim(9.9, 10.1)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    Height = np.zeros((200, 200))

    Index = 0

    for line in open("Output_" + str(i) + ".dat", 'r'):
        values = [float(s) for s in line.split()]

        Height[Index] = values
        Index += 1

    surf = ax.plot_surface(X, Y, Height, cmap=cm.coolwarm, linewidth=0, antialiased=False, vmax=9.999, vmin=10.001)
    fig.colorbar(surf, shrink=0.5, aspect=5)

    print(i)

    # plt.savefig("Images/Image_%03d.png" % i, dpi=300)
    plt.savefig("Images/Image_%03d.png" % i)
    surf.remove()
    plt.close(fig)
#-----------------------------#
subprocess.call(['./MakeVideo'])

for file_name in glob.glob("Images/*.png"):
        os.remove(file_name)
