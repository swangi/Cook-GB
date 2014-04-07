import time
import numpy as np
import matplotlib.pyplot as plt
# import GBFe100
f = open('all-GBs.txt', 'r')
data = np.genfromtxt('all-GBs.txt',names = True)
fig=plt.figure()
plt.ylabel('$\mathrm{Grain\ boundary\ energy\ (J/m^2)}$')
plt.xlabel('$\mathrm{Angle\ (degrees)}$')
plt.axis([0,90,0,1500])
plt.grid()
plt.plot(data['Angle'],data['GB_energy'],'-ro')
fig.savefig('allGB.png', transparent=True, bbox_inches='tight', \
                        pad_inches=0, dpi=300)
plt.show()

