#!/usr/bin/python3

#
# Gibs energies for CuAg system
# 
# data from NIST Solder database
# https://www.metallurgy.nist.gov/phase/solder/solder.html
#

import numpy as np
R = 8.314
T0 = 273.15
LN = np.log

def xlnx(x):
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(x, x*np.log(x), 0)


def G_Ag_ref(T):
    return np.where(T<1234.93, 
        -7209.512+118.202013*T-23.8463314*T*LN(T)-.001790585*T**2-3.98587E-07*T**3-12011*T**(-1),
        -15095.252+190.266404*T-33.472*T*LN(T)+1.412E+29*T**(-9)
    )

def G_Cu_ref(T):
    return np.where(T<1357.77,
        7209.512+118.202013*T-23.8463314*T*LN(T)-.001790585*T**2-3.98587E-07*T**3-12011*T**(-1),
        15095.252+190.266404*T-33.472*T*LN(T)+1.412E+29*T**(-9)
    )

#
# fcc phase
#

def L0_AgCu_fcc(T):
    return +36061.88-10.44288*T

def L1_AgCu_fcc(T):
    return -4310.12

def G_AgCu_fcc(T, x):
    xCu, xAg = x, 1-x
    return (
        xAg * G_Ag_ref(T)  + xCu * G_Cu_ref(T) # reference surface
        + R*T*(xlnx(xAg) + xlnx(xCu)) # mixing entropy
        + xAg*xCu*( L0_AgCu_fcc(T) + L1_AgCu_fcc(T)*(xAg-xCu) ) # excess terms
    ) 

#
# Liquid phase
#

def G_Ag_liq(T):
    return np.where(T<1235.08,
        +11025.076-8.891021*T-1.034E-20*T**7,
        +11508.141-9.301747*T-1.412E+29*T**(-9) 
    ) + G_Ag_ref(T)

def G_Cu_liq(T):
    return np.where(T<1358.02,
        +12964.736-9.511904*T-5.849E-21*T**7,
        +13495.481-9.922344*T-3.642E+29*T**(-9)
    ) + G_Cu_ref(T)

def L0_AgCu_liq(T):
    return +17323.40-4.46819*T

def L1_AgCu_liq(T):
    return +1654.38-2.35285*T

def G_AgCu_liq(T, x):
    xCu, xAg = x, 1-x
    return (
        xAg * G_Ag_liq(T)  + xCu * G_Cu_liq(T) # reference surface
        + R*T*(xlnx(xAg) + xlnx(xCu)) # mixing entropy
        + xAg*xCu*( L0_AgCu_liq(T) + L1_AgCu_liq(T)*(xAg-xCu) ) # excess terms
    ) 
