import matplotlib.pyplot as plt
import numpy as np

a = np.genfromtxt('cond_plot.txt', delimiter='\t')

plt.title("Conductivity Comparison")
plt.xlabel(r'$\sigma$  [S/cm]')
plt.ylabel(r'$z$  [km]')
plt.plot(a[:, 1], a[:,0],'r-', label = 'exp_conductivity') 
plt.plot(a[:, 2], a[:,0],'b-', label = 'conductivity')
plt.legend()
plt.grid()
plt.tight_layout()

plt.show()