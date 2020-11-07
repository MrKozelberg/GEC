#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

a = np.loadtxt('cond.txt')

plt.figure(figsize=(10, 5))

plt.title("Conductivity Parametrization")

plt.xlabel(r'conductivity [S/cm]')
plt.ylabel(r'altitude [km]')
plt.plot(a[:,1], a[:,0], 'r-', label = 'my realization')
plt.legend()
plt.grid()

plt.savefig('cond.png', bbox_inches = 'tight', dpi = 200)

plt.show()