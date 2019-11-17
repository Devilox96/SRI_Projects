import numpy as np
#-----------------------------#
DataPath = "./Data/"

xVel = np.fromfile(DataPath + "xVelocity.dat", dtype=float, sep="\t")
yVel = np.fromfile(DataPath + "yVelocity.dat", dtype=float, sep="\t")

ResVel = np.add(xVel / 2.0, yVel / 2.0)

xField = np.fromfile(DataPath + "xField.dat", dtype=float, sep="\t")
yField = np.fromfile(DataPath + "yField.dat", dtype=float, sep="\t")

ResField = np.add(xField / 2.0, yField / 2.0)
#-----------------------------#
#np.save(DataPath + "ResVel.dat", ResVel.reshape((-1, 4)), delimiter='\t')
#np.save(DataPath + "ResField.dat", ResField.reshape((-1, 4)), delimiter='\t')
np.save(DataPath + "ResVel.dat", ResVel.reshape((-1, 4)))
np.save(DataPath + "ResField.dat", ResField.reshape((-1, 4)))
