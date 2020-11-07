#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

a = np.loadtxt('cond.txt')

plt.figure(figsize=(10, 5))

plt.title("Conductivity Parametrization")

plt.xlabel(r'conductivity [S/cm]')
plt.ylabel(r'altitude [km]')
plt.plot(a[:,1], a[:,0], 'r-', label = '$xi = 0$')
plt.plot(a[:,2], a[:,0], 'b-', label = '$xi = 0.5$')
plt.plot(a[:,3], a[:,0], 'g-', label = '$xi = 1$')
plt.legend()
plt.grid()

plt.savefig('cond.png', bbox_inches = 'tight', dpi = 200)

plt.show()