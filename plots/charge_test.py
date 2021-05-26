import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

alt = np.load("alt.npy")
phi = np.load("phi.npy")

fig, ax = plt.subplots()

M = 360
N = 180

n = 80
m = 180

num = 2 * n * M + 2 * m + 1

ax.set_title("(2015-12-31) Potential of {} cell".format(num))
ax.set_ylabel(r'$z$ [km]')
ax.set_xlabel(r'$\phi$ [V]')

ax.plot(phi[num,:], alt, '-r', linewidth = 2)
ax.plot(phi[num-1,:], alt, '-b', linewidth = 2)