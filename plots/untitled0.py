import numpy as np
from sympy import *

z, sigma_0, j_0, H_0, H, sigma, j, S, n =  symbols('z, sigma_0, j_0, H_0, H, sigma, j, S, n')

sigma_0 = 6*10**(-14)
H_0 = 6
j_0 = 6.4*10**(-9)
H = 70

S = 1

def sigma(z):
    return sigma_0 *  exp(z / H_0)

def f1(z):
    return j_0 / sigma(z)

def f2(z):
    return 1 / sigma(z)

sum_up = summation(S *  integrate(f1(z), (z, 5, 10)) / integrate(f2(z), (z, 0, H)), (n,1,10*360))
sum_down = summation(S / integrate(f2(z), (z, 0, H)), (n,1,180*360))

ans = sum_up / sum_down

print(simplify(ans.evalf()))