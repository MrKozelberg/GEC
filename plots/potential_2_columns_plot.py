import numpy as np
import matplotlib.pyplot as plt

a = np.genfromtxt('potential_2_columns.txt', delimiter='\t')
c = np.genfromtxt('potential_2_columns_sourse.txt', delimiter='\t')
b = np.zeros(np.size(a[:,1]))
d = np.zeros(np.size(a[:,1]))
E = np.zeros(np.size(a[:,1]))
E_s = np.zeros(np.size(a[:,1]))

#dd
b[1] = a[1,1] - 2*a[2,1] + a[3,1]
b[np.size(a[:,1])-1] = a[np.size(a[:,1])-1,1] - 2*a[np.size(a[:,1])-2,1] + a[np.size(a[:,1])-2,1]
for i in range(2,np.size(a[:,1])-2):
    b[i] = a[i-1,1] - 2*a[i,1] + a[i+1,1]

d[1] = c[1,1] - 2*c[2,1] + c[3,1]
d[np.size(c[:,1])-1] = c[np.size(a[:,1])-1,1] - 2*c[np.size(a[:,1])-2,1] + c[np.size(a[:,1])-2,1]
for i in range(2,np.size(c[:,1])-2):
    d[i] = c[i-1,1] - 2*c[i,1] + c[i+1,1]

#d
E[1] = (a[3,1] - a[1,1])/2
E[np.size(a[:,1])-1] = (a[np.size(a[:,1])-1,1] - a[np.size(a[:,1])-2,1])/2
for i in range(2,np.size(a[:,1])-2):
    E[i] = (a[i+1,1] - a[i-1,1])/2

E_s[1] = (c[3,1] - c[1,1])/2
E_s[np.size(c[:,1])-1] = (c[np.size(c[:,1])-1,1] - c[np.size(c[:,1])-2,1])/2
for i in range(2,np.size(c[:,1])-2):
    E_s[i] = (c[i+1,1] - c[i-1,1])/2

plt.title("Potential of the Lower Atmosphere")
plt.xlabel(r'$\varphi$  [kV]')
plt.ylabel(r'$z$  [km]')
plt.plot(a[:, 1], a[:,0],'b-', label = 'column without sourses') 
plt.plot(c[:, 1], c[:,0],'r-', label = 'column with sourses')
plt.legend()
plt.grid()
plt.show()

plt.title("Elecric Field of the Lower Atmosphere")
plt.xlabel(r'$E_z$  [C]')
plt.ylabel(r'$z$  [km]')
plt.plot(-E[:], a[:,0],'b-', label = 'column without sourses') 
plt.plot(-E_s[:], c[:,0],'r-', label = 'column with sourses')
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

#charge
# Q_e = integrate(div(E),(z,0,12)km)*S_e
Q1 = 0.0
for k in range(0,20):
    Q1 = Q1 - b[k]

Q2 = 0.0
for k in range(0,20):
    Q2 = Q2 - d[k]

print("Charge per unit area in the column without sourses is about {0:6.2f} E-6 [C m^-2]".format(Q1))
print("Charge per unit area in the column with sourses is about {0:6.2f} E-6 [C m^-2]".format(Q2))
