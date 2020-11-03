#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

a = np.loadtxt('atm76.txt')
b = np.loadtxt('my_atm.txt')


# Линейная зависимость
x = np.linspace(0, 10, 50)
y1 = x

# Квадратичная зависимость
y2 = [i**2 for i in x]

# Построение графиков
plt.figure(figsize=(10, 5))

plt.suptitle("Standard Atmosphere")

plt.subplot(1, 2, 1)
plt.title("Pressure Comparison")
plt.xlabel(r'pressure [Pa]')
plt.ylabel(r'altitude [km]')
plt.plot( 101325*a[:, 2], a[:,0], 'b-', label = 'USSA 1976')
plt.plot(b[:,1], b[:,0], 'r-', label = 'my realization')
plt.legend()
plt.grid()
                
plt.subplot(1, 2, 2)
plt.title("Temperature Comparison")
plt.xlabel(r'temperature [K]')
plt.ylabel(r'altitude [km]')
plt.plot( 288.15*a[:, 3], a[:,0], 'b-', label = 'USSA 1976')
plt.plot(b[:,2], b[:,0], 'r-', label = 'my realization')
plt.legend()
plt.grid()

plt.savefig('my_atm.png', bbox_inches = 'tight', dpi = 200)