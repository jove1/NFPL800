#!/usr/bin/env python3

from pylab import *
from scipy.optimize import root

gca().add_artist(Circle((5,6), 5, fill=False))

x = np.array([0, 4])
plot(x, 3*x+2)

def func(p):
    print(".", end="", flush=True)
    x, y = p
    return (3*x+2-y, (x-5)**2 + (y-6)**2 - 25)

sol = root(func, x0=(1,1))
xs, ys = sol.x
plot([xs], [ys], "o")
print(xs, ys, func(sol.x))

gca().set_aspect("equal")
show()
