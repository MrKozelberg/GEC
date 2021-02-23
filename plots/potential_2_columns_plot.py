import numpy as np
import matplotlib.pyplot as plt

a = np.genfromtxt('potential_2_columns.txt', delimiter='\t')

plt.title("Potential of the Lower Atmosphere")
plt.xlabel(r'$\varphi$  [V]')
plt.ylabel(r'$z$  [km]')
plt.plot(a[:, 1], a[:,0],'r-', label = '*column number*') 
plt.xlabel(r'$\varphi$  [V]')
plt.ylabel(r'$z$  [km]')
plt.legend()
plt.grid()

plt.show()