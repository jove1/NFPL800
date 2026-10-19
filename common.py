
import numpy as np

class FT:
    def __init__(self, *args):
        self.limits = args[0::2]
        self.parts = args[1::2]

    def __call__(self, T, d=0):
        return np.select(
            [ T<lim for lim in self.limits[1:] ], 
            [ sum( f(T, d=d) * v for f, v in part.items() ) for part in self.parts ], 
            np.nan)

def C(T, d=0): 
    if d == 0: return 1 
    raise NotImplementedError
 
def T(T, d=0): 
    if d == 0: return T 
    raise NotImplementedError
 
def Tk(k): 
    def f(T, d=0): 
        if d == 0: return T**k 
        raise NotImplementedError
    return f 
 
def TLNT(T, d=0): 
    if d == 0: return T*np.log(T) 
    raise NotImplementedError
 
T2 = Tk(2) 
T3 = Tk(3) 
T4 = Tk(4) 
T7 = Tk(7) 
Tm1 = Tk(-1)
Tm9 = Tk(-9)

