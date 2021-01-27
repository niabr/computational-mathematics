#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
from math import sin,cos,log,exp,pi

Nt = 2*10**4 - 416
Nh = 100
L = 1
T = 1
dt = T/(Nt-1)
dh = L/(Nh-1)

u = np.zeros((Nt,Nh))
x = np.array([i*dh for i in range(Nh)])

def f(nt,nh):
    """Right hand side"""
    return pi**2*sin(pi*(nh*dh))

def Lambda(nt,nx):
    """Second difference"""
    ans = (u[nt][nx-1]-2*u[nt][nx]+u[nt][nx+1])/(dh**2)
    return ans

def Sol(nt,nx):
    """Exact solution"""
    ans = (1 - exp(-pi**2*nt*dt))*sin(pi*nx*dh)
    return ans

for i in range(Nt-1):
    for j in range(1,Nh-1):
        u[i+1][j] = u[i][j] + (Lambda(i,j) + f(i,j))*dt

eps = 0
for i in range(Nh-1):
    if Sol(Nt-1,i) == 0:
        continue
    else:
        eps = max(eps, abs(u[Nt-1][i]/Sol(Nt-1,i)-1))

y = [Sol(Nt-1,n) for n in range(Nh)]

print(eps)
print(Sol(Nt-1,Nh//2))
print(u[Nt-1][Nh//2])
"""print(u[Nt-1][Nh-5:Nh])
print((Nh-1)*dh)
print(sin(pi*(Nh-1)*dh))
print(u[Nt-1][Nh-1])"""

plt.plot(x,u[Nt-1],x,y)
plt.grid(which='both',axis='both')
plt.legend(['Последний слой, ε = %.8f' % eps, 'Разностная проекция точного решения'], loc='best',fontsize=18)
plt.show()

    
