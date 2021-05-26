import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

a = np.genfromtxt('phi_wos.txt')
#b = np.genfromtxt('phi_180*90+1.txt')

fig, ax = plt.subplots()

ax.set_xlim([0,250])
ax.set_ylim([0,70])
ax.set_ylabel(r'$z$ [km]')
ax.set_xlabel(r'$\varphi$ [V]')
ax.plot(a[:,1], a[:,0], 'b-', linewidth = 2, label = 'without sources')
#ax.plot(b[:,1], b[:,0], 'r-', linewidth = 2, label = 'with sources')

ax.legend()

ax.xaxis.set_major_locator(ticker.MultipleLocator(50))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(10))
ax.yaxis.set_major_locator(ticker.MultipleLocator(10))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(5))

ax.grid(which='major',
        color = 'k')


ax.minorticks_on()
ax.grid(which='minor',
        color = 'gray',
        linestyle = ':')

fig.set_figwidth(8)
fig.set_figheight(6)

plt.show()