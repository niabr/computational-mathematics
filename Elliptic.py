#!/usr/bin/python
import numpy as np
from math import log,exp
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter


eps = 1.0e-6
X = 1
Y = 1
def f(nx,ny):
    """Right hand side"""
    ans = -32*(nx*dx*(1-nx*dx) + ny*dy*(1-ny*dy))
    return ans

def vnorm(u):
    max_abs=0.0
    for x in range(u.size):
        max_abs+=abs(u[x])
    return max_abs


def Seidel(method):
    u=np.zeros((N**2))
    us=np.zeros((N**2))
    us1=np.zeros((N**2))
    count=0
    q=1.0

    while True:
        if method == 'regular':
            for m in range(1,N-1):
                for l in range(1,N-1):
                    us1[m*N+l] = ((us1[(m-1)*N+l]+us[(m+1)*N+l])/dx**2 + (us1[m*N+l-1]+us[m*N+l+1])/dy**2 - f(m,l))/(2/dx**2+2/dy**2)
        else:
            for m in range(1,N-1):
                for l in range(1+m%2,N-1,2):
                    us1[m*N+l] = ((us[(m-1)*N+l]+us[(m+1)*N+l])/dx**2 + (us[m*N+l-1]+us[m*N+l+1])/dy**2 - f(m,l))/(2/dx**2+2/dy**2)
            for m in range(1,N-1):
                for l in range(1+(m+1)%2,N-1,2):
                    us1[m*N+l] = ((us1[(m-1)*N+l]+us1[(m+1)*N+l])/dx**2 + (us1[m*N+l-1]+us1[m*N+l+1])/dy**2 - f(m,l))/(2/dx**2+2/dy**2)
            
        if vnorm(us-us1)<=eps:
            break
        # if vnorm(us-us1)==0:
        #     print('Error!')
        #     print(vnorm(us-us1),q)
        #     break
        #print(q)
        #q=vnorm(u-us)/vnorm(us-us1)
        u=us.copy()
        us=us1.copy()
        count+=1
    return (count,us1)
  
    
if __name__=="__main__":
    
    N_min = 5                   # Is a variable
    N_max = 41
    N_step = 5
    x = []
    y = []
    
    for N in range(N_min,N_max,N_step):
        x.append(N)
        dx = X/N
        dy = Y/N

        u=np.zeros((N**2))
        #count,ansSeidel=Seidel(method='regular')
        #print("Размер — %d, число итераций обычного метода Зейделя: " % N,count)
        count,ansSeidel = Seidel(method='chess')
        print("Размер — %d, число итераций шахматного метода Зейделя: " % N,count)
        y.append(count)

    X = np.arange(0, X, dx)
    Y = np.arange(0, Y, dy)
    X, Y = np.meshgrid(X, Y)
    Z = np.empty((N,N))
    for i in range(N**2):
        Z[i//N,i%N] = ansSeidel[i]
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                           linewidth=0, antialiased=True)
    # Customize the z axis.
    #ax.set_zlim(-1.01, 1.01)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)

    #print(x,y)
    xx = np.array(x)
    yy = np.array(y)
    logv = np.vectorize(log, otypes=[np.float])
    xx = logv(xx)
    yy = logv(yy)
    
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    A = np.vstack([xx, np.ones(len(x))]).T
    m, c = np.linalg.lstsq(A, yy, rcond=None)[0]
    ax2.set_title('Скорость сходимости, логарифмическая шкала. Степень %.2f, коэффициент %.2f' % (m, exp(c)), fontsize=20)

    plt.plot(xx, yy, 'o', label='Экспериментальные данные', markersize=10)
    plt.plot(xx, m*xx + c, 'r', label='Интерполяция')
    plt.legend()
    plt.grid(which='both',axis='both')
    #ax2.legend(['Последний слой, ε = %.8f' % eps, 'Разностная проекция точного решения'], loc='best',fontsize=18)
    plt.show()
    

