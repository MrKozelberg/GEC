import numpy as np
import matplotlib.pyplot as plt

a = np.genfromtxt('potential_2_columns.txt', delimiter='\t')
c = np.genfromtxt('potential_2_columns_sourse.txt', delimiter='\t')
b = np.zeros(np.size(a[:,1]))
d = np.zeros(np.size(a[:,1]))

b[0] = a[0,1] - 2*a[1,1] + a[2,1]
b[np.size(a[:,1])-1] = a[np.size(a[:,1])-1,1] - 2*a[np.size(a[:,1])-2,1] + a[np.size(a[:,1])-2,1]
for i in range(1,np.size(a[:,1])-2):
    b[i] = a[i-1,1] - 2*a[i,1] + a[i+1,1]

d[0] = c[0,1] - 2*c[1,1] + c[2,1]
d[np.size(c[:,1])-1] = c[np.size(a[:,1])-1,1] - 2*c[np.size(a[:,1])-2,1] + c[np.size(a[:,1])-2,1]
for i in range(1,np.size(a[:,1])-2):
    d[i] = c[i-1,1] - 2*c[i,1] + c[i+1,1]

plt.title("Potential of the Lower Atmosphere")
plt.xlabel(r'$\varphi$  [kV]')
plt.ylabel(r'$z$  [km]')
plt.plot(a[:, 1], a[:,0],'b-', label = 'column without sourses') 
plt.plot(c[:, 1], c[:,0],'r-', label = 'column with sourses')
plt.legend()
plt.grid()
plt.show()

plt.title("Charge Density of the Lower Atmosphere")
plt.xlabel(r'$\rho$  [C]')
plt.ylabel(r'$z$  [km]')
plt.plot(-b[:], a[:,0],'b-', label = 'column without sourses') 
plt.plot(-d[:], c[:,0],'r-', label = 'column with sourses')
plt.legend()
plt.grid()
plt.show()

Q1 = 0.0
for k in range(0,12):
    Q1 = Q1 - b[k]
Q1 = (180 - 20)*360*Q1

Q2 = 0.0
for k in range(0,12):
    Q2 = Q2 - d[k]
Q2 = 20*360*Q2

Q = Q1 + Q2
print("Charge of the lower atmosphere is about {0:6.2f} [C]".format(Q))
