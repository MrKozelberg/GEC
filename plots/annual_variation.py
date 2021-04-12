import matplotlib.pyplot as plt
import numpy as np

a = np.genfromtxt('annual_variation.txt', delimiter='\t')

plt.title("Annual Variation")
plt.xlabel(r'$t$  [Years]')
plt.ylabel(r'$\varphi$  [V/m]')
#plt.plot(a[:, 1], a[:,0],'r-', label = 'exp_conductivity') 
plt.plot(a[:, 0], a[:,1],'b-', label = 'potential')
plt.legend()
plt.grid()
plt.tight_layout()

plt.show()