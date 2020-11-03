import numpy as np
import matplotlib.pyplot as plt

a = np.genfromtxt('potential_2_columns.txt', delimiter='\t')

plt.title("Potential of the Low Atmosphere")
plt.xlabel(r'$\varphi$  [V]')
plt.ylabel(r'$z$  [km]')
plt.plot(a[:, 1], a[:,0],'r-', label = 'with a source') 
plt.xlabel(r'$\varphi$  [V]')
plt.ylabel(r'$z$  [km]')
plt.plot(a[:, 2], a[:,0],'b-', label = 'without any sources')
plt.legend()
plt.grid()
plt.savefig('potential_2_columns.png', bbox_inches = 'tight', dpi = 128) 