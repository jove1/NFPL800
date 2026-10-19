#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.rc("axes", grid=True)

from common import *


func = FT(0, {T:1, T2:1e-3, TLNT:2}, 1000)

t = np.linspace(0, 1000)

plt.figure()
plt.plot(t, func(t))

eps = 1e-4

plt.figure()
plt.plot(t, (func(t+eps)-func(t-eps))/2/eps, ".")
plt.plot(t, func(t, d=1))

plt.figure()
plt.plot(t, NotImplementedError, ".")
plt.plot(t, func(t, d=2))



plt.show()
