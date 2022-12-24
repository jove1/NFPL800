#!/usr/bin/python3

from pylab import *

from scipy.integrate import solve_ivp

g, vx0, vy0 = -10, 10, 10

def func(t, y):
    print(".", end="", flush=True)
    x, y, vx, vy = y
    return (vx, vy, 0, g)

sol = solve_ivp(func, (0,3), y0=(0, 0, vx0, vy0), max_step=0.05 )
x, y, _, _ = sol.y
plot(x, y, ".")

x = linspace(0,30)
t = x/vx0
plot(x, t*vy0 + 0.5*g*t*t, "-")

gca().set_aspect("equal")
show()
