#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.rc("axes", grid=True)

from common import wt2at, at2wt

print( at2wt(["H", "O"], [2, 1]) )
print("should be:", [0.11190053, 0.88809947])

print( wt2at(["H", "O"], [0.11190053, 0.88809947]) )
print("should be:", [0.66666666, 0.33333334])

x = np.linspace(0, 1)
plt.plot(x, wt2at(["Cu", "Al"], [x, 1-x]) )
plt.show()
