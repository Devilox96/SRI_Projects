from matplotlib import pyplot as plt
import numpy as np

FrequiencyArray = []
AbsorbtionArray = []

Temperature = 182.5
LeftLimit = 5050
RightLimit = 5160
AbsorbtionLimit = 4.0e-3

for line in open(str(Temperature) + '.txt', 'r'):
    values = [float(s) for s in line.split()]

    if (values[0] > LeftLimit and values[0] < RightLimit and values[1] < AbsorbtionLimit):
        FrequiencyArray.append(values[0])
        AbsorbtionArray.append(values[1])

FrequiencyArray = np.asarray(FrequiencyArray)
AbsorbtionArray = np.asarray(AbsorbtionArray)

plt.plot(FrequiencyArray, AbsorbtionArray)
plt.show()