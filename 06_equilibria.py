#!/usr/bin/python3

from pylab import *

#
# plot phase diagram from previous calculations as background
#

data = 'H4sIAO4cp2MC/52YC1BUVRjHQVLToSzMSO0BjTH2kMHEInWEMspxk+FhoCmjyzvksSwLtEoIgjyWFbbh2ohPBsGsxBRUTEHQLAbMBAMMdRADQQQWzEfSjNU9507Nufd+Mx9z7jDszJnvd8695/y/10l/bFvTczb0CRZmm6fFJcXqjK6h8fpw19ikGMOnWr1eaxTMExJDtTFavVAkmMdTE/E3zGDUhZORcREegsZGY5sj+AuzNXZm20WCj4/Psn/Eh/7T2BiEkCWPS6uMeArZol2UXdSk/4ZsbFRD7l50KEc2iKFeFWtf9R3w5Ee9edAbmVMWP+zS8K/qy4NatlfPm9UbwP/CgTzo3PjIpT9fXcWDZnxS51/SvoYHbfF2DDa0rOXepggtD/pOQbdHQXcI9w6vC+NB84wLl4R0hnMfzupIHrQjrGj5Wx1R3OcaGC1Hf/rzw8kuiTo1KptN1H7Q5Nb1PCiVhF8sD0o3+Jc4HpSqyUXHg9IN7uFCqRD3JfCgT1A56XlQqmGnRB70WFq5V8V1EE1HUPq+eww8qBgiRFEkydGndANV068Ps3bAbNRzZibzoOOjHy3LuIKjX99dP+X7g7meSqcrTuFBad4I+owHpf76rHGM32phhz4mTxsXSl1d2MCDivsrvvJGFD1xv8M8XCioosTTqTxoOcnNzalj3OHtqgBT8Lkcdepqj5x3GkPFGCGKIg1Fd9x59MXmxl3skDMRsf0myHNKWDtgtvuviypu2oSuCsxGw9r76Sha98By3pRRyg5dJNG0TIHO9T07c/lu4HDKVBHRLoMH3ZfvevJWCI4evveyz4o9+1XB9CyOSud6gB1KyT5Xs9t5sxwt/eNo8wunvmHt3q2vuBiROqyKw+kKVDr/gwiqIY7TpUCLRj4I6Lt8CEFpCPfMhNR0mLXzcyvOTF8tQ19KHTkzca8CfS9lMHrlhiMISqP/3wpUkk4lgt5NFrc4OEuOLkye/0ZtbhVrF1KS9faOBcOqxFGrQCXpHEXQhoTnxb8tclSSznHWTk8jx7Aq52xUoHMMDYeyy6sRdGdMZaPLVQUqbsjeCw+taiGeUKWrBdmQEE+ydsBstPQvzobUJLN7Ub9KdIpTqkw3qkAdE8i6tawdMBup1Swrc6BVT7N2ZevenGhbbVUlyWoFSr50xpdWdVirZ4ccqSZy5WhPJVHAFdYOmI3m12QFmjt121dbrL+xdkRLdQ/qUVRSUztrJ6WhOnaIFNJ/tSnQH70mxbsnyd7u2pFvtf3e13jQNvomnaqCYH4eejiSM7WiKLBq32Wipy5ZLhH7yO+EPChKXFInyRoULZeUyNoBHkErmHsKNDToh0v5HjdYO+AbQBRwTuk9Lsg8glQEASY5StOr5+/IqiAaQ8XejaC05Ko0QfJvUifJ4yja773rwKI5PciqJL5Md8hH43Cce1LakHMDio4ey1kTaH+OtQPSEK0RE/LRKEFTydSbKArkV8DVaQ/ajK9q90pNS+voTRQFgi5QXtCi1s0sRwuf+ajXqasXCWsgCtRXgCTIbVNmoRnKdH1ICAfRpl/JMd5CUFqFW81oBQMUtSAKLEG80N+tnx0iV0YePluhEL4f6ThAFMiI06gm+lVtQ8VW6HAEpB0EURLU/Iy3kYKAXEAO2ReM8YZAdmkEojqq2AFkVdrnxBRA8h9A6iZS+peeV6BAVSclyUFVi/RaoRw1ZTh0Rt0eRGpE2qqbFOjSM2FPjosbQoraRhex1OtXoBNsqxfvvDOEFPAzHMSKeJYFza9A/0IbsxUWtLwE2rws6rAWtGsGvpX2dFU4Ctz8WO3FxqH7fzTc9V8qdrWQ1RgAAA=='

import base64, gzip, pickle
for x1, x2, T in pickle.loads(gzip.decompress(base64.b64decode(data))):
    plot([x1, x2], [T, T], "k-", lw=0.2)

#
# now start with the equilibria calculations
#

from dataCuAg import *
from dataCuAg_deriv import *
from scipy.optimize import root

#
# Melting points
#

def func(T, x):
    print(".", end="", flush=True)
    return G_AgCu_fcc(T0+T, x) - G_AgCu_liq(T0+T, x)

sol = root(func, 1000, args=(0,))
print(sol.message)
TAg = sol.x[0]
print("Ag melting point", TAg)
plot([0], [TAg], "o")
        
sol = root(func, 1000, args=(1,))
print(sol.message)
TCu = sol.x[0]
print("Cu melting point", TCu)
plot([1], [TCu], "o")


#
# Solidus and liquidus temperature at given composition
#

def func(x1, x2, T, a, b):
    print(".", end="", flush=True)
    return [
            a+b*x1 - G_AgCu_fcc(T0+T, x1), 
            a+b*x2 - G_AgCu_liq(T0+T, x2), 
            b - dx_G_AgCu_fcc(T0+T, x1), 
            b - dx_G_AgCu_liq(T0+T, x2), 
    ]

xx, T = 0.09, 800
axvline(xx, ls="--")

b0 = dx_G_AgCu_liq(T0+T, xx)
a0 = G_AgCu_liq(T0+T, xx) - xx*b0

sol = root( lambda p: func(p[0], xx, p[1], p[2], p[3]), x0=[0.05, T, a0, b0])
print(sol.message)
x1, x2, T = sol.x[0], xx, sol.x[1]
print("Liquidus", x1, x2, T)
plot([x1, x2], [T, T], "o-")

sol = root( lambda p: func(xx, p[0], p[1], p[2], p[3]), x0=[0.15, T, a0, b0])
print(sol.message)
x1, x2, T = xx, sol.x[0], sol.x[1]
print("Solidus", x1, x2, T)
plot([x1, x2], [T, T], "o-")


#
# Eutectic point
#

def func(x):
    print(".", end="", flush=True)
    x1, x2, x3, T, a, b = x
    return [
        a+b*x1 - G_AgCu_fcc(T0+T, x1), 
        a+b*x2 - G_AgCu_fcc(T0+T, x2), 
        a+b*x3 - G_AgCu_liq(T0+T, x3), 
        b - dx_G_AgCu_fcc(T0+T, x1),  
        b - dx_G_AgCu_fcc(T0+T, x2), 
        b - dx_G_AgCu_liq(T0+T, x3),
    ]

x3, T = 0.5, 800
b0 = dx_G_AgCu_liq(T0+T, x3)
a0 = G_AgCu_liq(T0+T, x3) - x3*b0

if 1:
    # this is very dependent on starting conditions,
    # often jumps into x<0 or x>1, where entropy term is not defined
    x0 = (0.1, 0.9, x3, T, a0, b0)
    sol = root( func, x0=x0)
else:
    # (ab)use constrained optimization, to set bounds on x variables
    x0 = (0.25, 0.8, x3, T, a0, b0)
    from scipy.optimize import minimize, Bounds
    sol = minimize(lambda p: 0, x0=x0,
        constraints=[dict(type="eq", fun=func)], 
        bounds=Bounds([0,0,0,-np.inf,-np.inf,-np.inf],[1,1,1,np.inf,np.inf,np.inf]))

print(sol.message)
x1, x2, x3, T = sol.x[:4] # don't care about a, b
print("Eutectic point", x1, x2, x3, T)
plot([x1, x2, x3], [T, T, T], "o-")

xlim(0, 1)
grid(True)
show()
