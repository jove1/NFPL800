#!/usr/bin/python3

#
# derivatives of Gibs energies for CuAg system
# 

from dataCuAg import *

def dG_Ag_ref(T):
    return np.where(T<1234.93, 
        +118.202013-23.8463314*(LN(T)+1)-.001790585*2*T-3.98587E-07*3*T**2-12011*-1*T**(-2),
        +190.266404-33.472*(LN(T)+1)+1.412E+29*-9*T**(-10)
    )

def dG_Cu_ref(T):
    return np.where(T<1357.77,
        +118.202013-23.8463314*(LN(T)+1)-.001790585*2*T-3.98587E-07*3*T**2-12011*-1*T**(-2),
        +190.266404-33.472*(LN(T)+1)+1.412E+29*-9*T**(-10)
    )

#
# fcc phase
#

def dL0_AgCu_fcc(T):
    return -10.44288

def dL1_AgCu_fcc(T):
    return 0

def dT_G_AgCu_fcc(T, x):
    xCu, xAg = x, 1-x
    return (
        xAg * dG_Ag_ref(T)  + xCu * dG_Cu_ref(T) # reference surface
        + R*(xlnx(xAg) + xlnx(xCu)) # mixing entropy
        + xAg*xCu*( dL0_AgCu_fcc(T) + dL1_AgCu_fcc(T)*(xAg-xCu) ) # excess terms
    ) 

def dx_G_AgCu_fcc(T, x):
    return (
        -G_Ag_ref(T) + G_Cu_ref(T)
        + R*T*(np.log(x) - np.log(1-x))
        + (1-2*x)*L0_AgCu_fcc(T)
        + (6*x*x - 6*x + 1)*L1_AgCu_fcc(T)
    )

def dxdT_G_AgCu_fcc(T, x):
    return (
        -dG_Ag_ref(T) + dG_Cu_ref(T)
        + R*(np.log(x) - np.log(1-x))
        + (1-2*x)*dL0_AgCu_fcc(T)
        + (6*x*x - 6*x + 1)*dL1_AgCu_fcc(T)
    )

def dxdx_G_AgCu_fcc(T, x):
    return (
        + R*T*(1/x + 1/(1-x))
        + -2*L0_AgCu_fcc(T)
        + (12*x - 6)*L1_AgCu_fcc(T)
    )


#
# Liquid phase
#

def dG_Ag_liq(T):
    return np.where(T<1235.08,
        -8.891021-1.034E-20*7*T**6,
        -9.301747-1.412E+29*-9*T**(-10) 
    ) + dG_Ag_ref(T)

def dG_Cu_liq(T):
    return np.where(T<1358.02,
        -9.511904-5.849E-21*7*T**6,
        -9.922344-3.642E+29*-9*T**(-10)
    ) + dG_Cu_ref(T)

def dL0_AgCu_liq(T):
    return -4.46819

def dL1_AgCu_liq(T):
    return -2.35285

def dT_G_AgCu_liq(T, x):
    xCu, xAg = x, 1-x
    return (
        xAg * dG_Ag_liq(T)  + xCu * dG_Cu_liq(T) # reference surface
        + R*(xlnx(xAg) + xlnx(xCu)) # mixing entropy
        + xAg*xCu*( dL0_AgCu_liq(T) + dL1_AgCu_liq(T)*(xAg-xCu) ) # excess terms
    ) 

def dx_G_AgCu_liq(T, x):
    return (
        -G_Ag_liq(T) + G_Cu_liq(T)
        + R*T*(np.log(x) - np.log(1-x))
        + (1-2*x)*L0_AgCu_liq(T)
        + (6*x*x - 6*x + 1)*L1_AgCu_liq(T)
    )

def dxdT_G_AgCu_liq(T, x):
    return (
        -dG_Ag_liq(T) + dG_Cu_liq(T)
        + R*(np.log(x) - np.log(1-x))
        + (1-2*x)*dL0_AgCu_liq(T)
        + (6*x*x - 6*x + 1)*dL1_AgCu_liq(T)
    )

def dxdx_G_AgCu_liq(T, x):
    return (
        + R*T*(1/x + 1/(1-x))
        + -2*L0_AgCu_liq(T)
        + (12*x - 6)*L1_AgCu_liq(T)
    )


