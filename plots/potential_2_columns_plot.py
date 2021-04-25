import numpy as np
import matplotlib.pyplot as plt

phi = np.genfromtxt('potential_2_columns.txt', delimiter='\t')
phi_s = np.genfromtxt('potential_2_columns_sourse.txt', delimiter='\t')


E = np.zeros(np.size(phi[:,1]))
E_s = np.zeros(np.size(phi_s[:,1]))
for i in range(0,np.size(phi[:,1])-1):
    E[i] = - (phi[i+1,1] - phi[i-1,1]) / 2 
    E_s[i] = - (phi_s[i+1,1] - phi_s[i-1,1]) / 2
    if (i == 0):
        E[i] = - (phi[i+1,1] - phi[i,1])
        E_s[i] = - (phi_s[i+1,1] - phi_s[i,1])
    if (i == np.size(phi[:,1])-1):
        E[i] = - (phi[i,1] - phi[i-2,1]) / 2
        E_s[i] = - (phi_s[i,1] - phi_s[i-2,1]) / 2

plt.title("Potential of the Lower Atmosphere")
plt.xlabel(r'$\varphi$  [kV]')
plt.ylabel(r'$z$  [km]')
plt.plot(phi[:, 1], phi[:,0],'b-', label = 'column without sourses') 
plt.plot(phi_s[:, 1], phi_s[:,0],'r-', label = 'column with sourses')
plt.legend()
plt.grid()
plt.show()

plt.title("Electric Intensity of the Lower Atmosphere")
plt.xlabel(r'$E$ $[{V}{m^{-1}}]$')
plt.ylabel(r'$z$  $[km]$')
plt.plot(E, phi[:,0],'b-', label = 'column without sourses') 
plt.plot(E_s, phi[:,0],'r-', label = 'column with sourses')
plt.legend()
plt.grid()
plt.show()

