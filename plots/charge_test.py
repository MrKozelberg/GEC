a = np.genfromtxt('cond.txt', delimiter='\t')

z = a[:,0]
sigma = a[:,1]
mu = a[:,2]

E = np.zeros(np.size(z))
E[0] = 100 # [V/m]
for i in range(1,np.size(E)):
    E[i] = E[i-1]*sigma[i]/sigma[i-1]

u = mu*E

b = np.ones(np.size(z))
for i in range(1,np.size(z)-1):
    if (i % 2 == 1):
        b[i] = 4
    else:
        b[i] = 2

print('Время, за которое положительные ионы преодолевают тропосферу:', sum(b/u)*1000)
print('тогда заряд атмосферы будет оцениваться в диапазоне от ', 0.016*sum(b/u)*1000, ' до ', 0.08*sum(b/u)*1000, ' [Кл]')