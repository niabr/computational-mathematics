#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
from math import sin,cos,log,exp,pi

k = 0
a = 0.5 + k/10
Nti = 101
Nhi = 100
L = 1
T = 1


def Sol(nt,nx):
    """Exact solution"""
    ans = 1/(nt*dt+nx*dx+a)
    return ans

def Lambda(nt,nx):
    """Second difference"""
    ans = (u[nt,nx-1]-2*u[nt,nx]+u[nt,nx+1])/dx**2
    return ans

for z in range(3):
    if z==0:
        Nt = Nti
        Nh = Nhi
    elif z==1:
        Nt = Nti*2 - 1
        Nh = Nhi*2 - 1
    elif z==2:
        Nt = Nti*4 - 3
        Nh = Nhi*4 - 3
    # elif z==3:
    #     Nh = Nhi*2 - 1
    #     Nt = Nti*2 - 1
    # elif z==4:
    #     Nh = Nhi*4 - 3
    #     Nt = Nti*4 - 3
    dt = T/(Nt-1)
    dx = L/(Nh-1)

    if z!=0:
        del u,x
    u = np.zeros((Nt,Nh))
    x = np.array([i*dx for i in range(Nh)])


    """Initial value"""
    for i in range(Nh):
        u[0][i] = Sol(0,i)

    """First step"""
    for i in range(1,Nh-1):
        u[1][i] = u[0][i] - 1/(i*dx+a)**2*dt + Lambda(0,i)*dt**2/2

    for i in range(1,Nt-1):
        """Boundary condition"""
        u[i][0] = Sol(i,0)
        u[i][Nh-1] = Sol(i,Nh-1)
        for j in range(1,Nh-1):
            u[i+1][j] = 2*u[i][j] - u[i-1,j] + Lambda(i,j)*dt**2

    u[Nt-1,0] = Sol(Nt-1,0)
    u[Nt-1,Nh-1] = Sol(Nt-1,Nh-1)
    eps = 0
    for i in range(Nh-1):
        if Sol(Nt-1,i) == 0:
            continue
        else:
            eps = max(eps, abs(u[Nt-1][i]/Sol(Nt-1,i)-1))

    y = [Sol(Nt-1,n) for n in range(Nh)]

    if z==0:
        u0 = u[Nt-1][Nhi//2]
        #print(Nhi//2*dh)
    elif z==1:
        u1 = u[Nt-1][(Nhi//2)*2]
        #print((Nhi//2)*2*dh)
    elif z==2:
        u2 = u[Nt-1,(Nhi//2)*4]
        #print((Nhi//2)*4*dh)
    # elif z==3:
    #     u2 = u[Nt-1][Nh//2]
    # else:
    #     u22 = u[Nt-1,Nh//2]

p = log((u1-u0)/(u2-u1))/log(2)
# ph = (u2-u0)/(u22-u2)
# print(u0,u1,u12,u2,u22)
#print(p)
#print(u0,u1,u2)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('Экспериментальное значение порядка аппроксимации p = %.2f' % p, fontsize=20)

ax.plot(x,u[Nt-1],x,y)
ax.grid(which='both',axis='both')
ax.legend(['Последний слой, ε = %.8f' % eps, 'Разностная проекция точного решения'], loc='best',fontsize=18)
plt.show()

    
