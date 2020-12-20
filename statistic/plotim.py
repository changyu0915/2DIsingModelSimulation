#coding by @Changyu 2020.12.20
#plot

import numpy as np
import matplotlib.pyplot as plt 

lat_05 = np.loadtxt("/Users/imac/Desktop/github/2DIsingModelSimulation/Statistic/lattice_05.text")
plt.subplot(1,5,1)
plt.imshow(lat_05,cmap=plt.cm.gray)

lat_16 = np.loadtxt("/Users/imac/Desktop/github/2DIsingModelSimulation/Statistic/lattice_16.text")
plt.subplot(1,5,2)
plt.imshow(lat_16,cmap=plt.cm.gray_r)

lat_27 = np.loadtxt("/Users/imac/Desktop/github/2DIsingModelSimulation/Statistic/lattice_27.text")
plt.subplot(1,5,3)
plt.imshow(lat_27,cmap=plt.cm.gray_r)

lat_38 = np.loadtxt("/Users/imac/Desktop/github/2DIsingModelSimulation/Statistic/lattice_38.text")
plt.subplot(1,5,4)
plt.imshow(lat_38,cmap=plt.cm.gray_r)

lat_50 = np.loadtxt("/Users/imac/Desktop/github/2DIsingModelSimulation/Statistic/lattice_50.text")
plt.subplot(1,5,5)
plt.imshow(lat_50,cmap=plt.cm.gray_r)

plt.show()

plt.savefig("/Users/imac/Desktop/github/2DIsingModelSimulation/Statistic/Change_with_Tmp")

