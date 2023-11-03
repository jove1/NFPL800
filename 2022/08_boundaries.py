#!/usr/bin/env python3

from pylab import *

#
# Plot phase diagram from previous calculations as background
#

data = 'H4sIAO4cp2MC/52YC1BUVRjHQVLToSzMSO0BjTH2kMHEInWEMspxk+FhoCmjyzvksSwLtEoIgjyWFbbh2ohPBsGsxBRUTEHQLAbMBAMMdRADQQQWzEfSjNU9507Nufd+Mx9z7jDszJnvd8695/y/10l/bFvTczb0CRZmm6fFJcXqjK6h8fpw19ikGMOnWr1eaxTMExJDtTFavVAkmMdTE/E3zGDUhZORcREegsZGY5sj+AuzNXZm20WCj4/Psn/Eh/7T2BiEkCWPS6uMeArZol2UXdSk/4ZsbFRD7l50KEc2iKFeFWtf9R3w5Ee9edAbmVMWP+zS8K/qy4NatlfPm9UbwP/CgTzo3PjIpT9fXcWDZnxS51/SvoYHbfF2DDa0rOXepggtD/pOQbdHQXcI9w6vC+NB84wLl4R0hnMfzupIHrQjrGj5Wx1R3OcaGC1Hf/rzw8kuiTo1KptN1H7Q5Nb1PCiVhF8sD0o3+Jc4HpSqyUXHg9IN7uFCqRD3JfCgT1A56XlQqmGnRB70WFq5V8V1EE1HUPq+eww8qBgiRFEkydGndANV068Ps3bAbNRzZibzoOOjHy3LuIKjX99dP+X7g7meSqcrTuFBad4I+owHpf76rHGM32phhz4mTxsXSl1d2MCDivsrvvJGFD1xv8M8XCioosTTqTxoOcnNzalj3OHtqgBT8Lkcdepqj5x3GkPFGCGKIg1Fd9x59MXmxl3skDMRsf0myHNKWDtgtvuviypu2oSuCsxGw9r76Sha98By3pRRyg5dJNG0TIHO9T07c/lu4HDKVBHRLoMH3ZfvevJWCI4evveyz4o9+1XB9CyOSud6gB1KyT5Xs9t5sxwt/eNo8wunvmHt3q2vuBiROqyKw+kKVDr/gwiqIY7TpUCLRj4I6Lt8CEFpCPfMhNR0mLXzcyvOTF8tQ19KHTkzca8CfS9lMHrlhiMISqP/3wpUkk4lgt5NFrc4OEuOLkye/0ZtbhVrF1KS9faOBcOqxFGrQCXpHEXQhoTnxb8tclSSznHWTk8jx7Aq52xUoHMMDYeyy6sRdGdMZaPLVQUqbsjeCw+taiGeUKWrBdmQEE+ydsBstPQvzobUJLN7Ub9KdIpTqkw3qkAdE8i6tawdMBup1Swrc6BVT7N2ZevenGhbbVUlyWoFSr50xpdWdVirZ4ccqSZy5WhPJVHAFdYOmI3m12QFmjt121dbrL+xdkRLdQ/qUVRSUztrJ6WhOnaIFNJ/tSnQH70mxbsnyd7u2pFvtf3e13jQNvomnaqCYH4eejiSM7WiKLBq32Wipy5ZLhH7yO+EPChKXFInyRoULZeUyNoBHkErmHsKNDToh0v5HjdYO+AbQBRwTuk9Lsg8glQEASY5StOr5+/IqiAaQ8XejaC05Ko0QfJvUifJ4yja773rwKI5PciqJL5Md8hH43Cce1LakHMDio4ey1kTaH+OtQPSEK0RE/LRKEFTydSbKArkV8DVaQ/ajK9q90pNS+voTRQFgi5QXtCi1s0sRwuf+ajXqasXCWsgCtRXgCTIbVNmoRnKdH1ICAfRpl/JMd5CUFqFW81oBQMUtSAKLEG80N+tnx0iV0YePluhEL4f6ThAFMiI06gm+lVtQ8VW6HAEpB0EURLU/Iy3kYKAXEAO2ReM8YZAdmkEojqq2AFkVdrnxBRA8h9A6iZS+peeV6BAVSclyUFVi/RaoRw1ZTh0Rt0eRGpE2qqbFOjSM2FPjosbQoraRhex1OtXoBNsqxfvvDOEFPAzHMSKeJYFza9A/0IbsxUWtLwE2rws6rAWtGsGvpX2dFU4Ctz8WO3FxqH7fzTc9V8qdrWQ1RgAAA=='

import base64, gzip, pickle
for x1, x2, T in pickle.loads(gzip.decompress(base64.b64decode(data))):
    plot([x1, x2], [T, T], "k-", lw=0.2)

TAg = 961.7795532622599
plot([0], [TAg], "o")

TCu = 1084.6199854338533
plot([1], [TCu], "o")

x1, x2, T = 0.04245288466265213, 0.09, 912.4762124401511
plot([x1, x2], [T, T], "o-")

x1, x2, T = 0.09, 0.21369526454692406, 851.1059357249133
plot([x1, x2], [T, T], "o-")

x1, x2, x3, T = 0.13272305627978498, 0.9546058274141793, 0.4025244602393358, 779.140887553754
plot([x1, x2, x3], [T, T, T], "o-")

#
# Now start with the phase boundary tracing
#

from dataCuAg import *
from dataCuAg_deriv import *
from scipy.integrate import solve_ivp

def func(T, x, *args):
    print(".", end="", flush=True)
    x1, x2 = x
    (dT1, dxdT1, dxdx1), (dT2, dxdT2, dxdx2) = args
    db = (dT1(T, x1) - dT2(T, x2))/(x1-x2)
    ret = (
        (db - dxdT1(T, x1))/dxdx1(T, x1), 
        (db - dxdT2(T, x2))/dxdx2(T, x2),
    )
    return ret

der_fcc = dT_G_AgCu_fcc, dxdT_G_AgCu_fcc, dxdx_G_AgCu_fcc
der_liq = dT_G_AgCu_liq, dxdT_G_AgCu_liq, dxdx_G_AgCu_liq

print("fcc1+fcc2")
sol = solve_ivp(func, (T0+T, T0), [x1, x2], args=(der_fcc, der_fcc), max_step=5)
print(sol.message)

plot(sol.y[0], sol.t-T0, "-")
plot(sol.y[1], sol.t-T0, "-")

print("fcc1+liq")
sol = solve_ivp(func, (T0+T, T0+TCu-1), [x2, x3], args=(der_fcc, der_liq), max_step=5)
print(sol.message)

plot(sol.y[0], sol.t-T0, "-")
plot(sol.y[1], sol.t-T0, "-")

print("fcc2+liq")
sol = solve_ivp(func, (T0+T, T0+TAg-1), [x1, x3], args=(der_fcc, der_liq), max_step=5)
print(sol.message)

plot(sol.y[0], sol.t-T0, "-")
plot(sol.y[1], sol.t-T0, "-")

xlim(0, 1)
grid(True)
show()
