#!/usr/bin/python3
from pylab import *
from scipy.spatial import ConvexHull

def xlnx(x):
    with np.errstate(divide='ignore', invalid='ignore'):
        return where(x, x*log(x), 0)

T0 = 273.15
R = 8.314
def Gmix(x, T):
    return T*R*( xlnx(x) + xlnx(1-x) )

N = 150
x = linspace(0, 1, N)

plot(x, Gmix(x, T0+20)  )
plot(x, Gmix(x, T0+20)  + 1000*x*(1-x) )
plot(x, Gmix(x, T0+20)  + 3000*x*(1-x)*(1-x-x) )
show()

