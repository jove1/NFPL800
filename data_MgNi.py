
try:
    import autograd.numpy as np
except ImportError:
    import numpy as np

R = 8.314
def xlnx(x):
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(x, x*np.log(x), 0)

class Phase:
    pass

#
# FUN GHSERMG 298.15 -8367.34+143.675547*T-26.1849782*T*LN(T)
#     +4.858E-04*T**2-1.393669E-06*T**3+78950*T**(-1); 923 Y
#      -14130.185+204.716215*T-34.3088*T*LN(T)
#
def GHSERMG(T):
    return np.where(T < 923.0, -8367.34 + 143.675547*T + -26.1849782*T*np.log(T) + 0.0004858*T**2 + -1.393669e-06*T**3 + 78950.0*T**-1, 
                             -14130.185 + 204.716215*T + -34.3088*T*np.log(T) + 1.038192e+28*T**-9)
#
# FUN GLIQMG 298.15 +8202.243-8.83693*T-8.01759E-20*T**7
#     +GHSERMG#; 923 Y
#      -5439.869+195.324057*T-34.3088*T*LN(T); 3000 N !
#
def GLIQMG(T):
    return np.where(T < 923.0, 8202.243 + -8.83693*T + -8.01759e-20*T**7 + GHSERMG(T), 
                              -5439.869 + 195.324057*T + -34.3088*T*np.log(T))
#
# FUN GHSERNI 298.15 -5179.159+117.854*T-22.096*T*LN(T)
#     -.0048407*T**2; 1728 Y
#      -27840.655+279.135*T-43.1*T*LN(T)+1.12754E+31*T**(-9); 3000 N !
#
def GHSERNI(T):
    return np.where(T < 1728.0, -5179.159 + 117.854*T + -22.096*T*np.log(T) + -0.0048407*T**2, 
                               -27840.655 + 279.135*T + -43.1*T*np.log(T) + 1.12754e+31*T**-9)
#
# FUN GLIQNI 298.15 +16414.686-9.397*T-3.82318E-21*T**7
#     +GHSERNI#; 1728 Y
#      +18290.88-10.537*T-1.12754E+31*T**(-9)+GHSERNI#; 6000 N !
#
def GLIQNI(T):
    return np.where(T < 1728.0, 16414.686 + -9.397*T + -3.82318e-21*T**7 + GHSERNI(T), 
                                18290.88 + -10.537*T + -1.12754e+31*T**-9 + GHSERNI(T))
#
# FUN GFCCMG 298.15 +2600-.9*T+GHSERMG#; 6000 N !
#
def GFCCMG(T):
    return (2600.0 + -0.9*T + GHSERMG(T))

#
# PHASE LIQUID:L %  1  1.0
#  > Metallic Liquid Solution: Redlich-Kister_Muggianu Model. !
# CONST LIQUID:L :AG,AL,....,MG,....,NI,...,ZN,ZR : !
#
# (Mg,Ni)1
#
MgNi_liq = Phase()
MgNi_liq.X = lambda y: y[1]
MgNi_liq.G = lambda T, y: (
#
# (mechanical mixture)
#
# PARAM G(LIQUID,MG;0) 298.15 +GLIQMG#; 3000 N 91Din !
# PARAM G(LIQUID,NI;0) 298.15 +GLIQNI#; 6000 N 91Din !
# 
        y[0] * (GLIQMG(T)) + 
        y[1] * (GLIQNI(T)) + 
#
# (subregular solution interaction)
#
# PARAM G(LIQUID,MG,NI;0) 298.15 -42304.49+7.45704*T; 6000 N 98Jac2 !
# PARAM G(LIQUID,MG,NI;1) 298.15 -15611.66+9.11885*T; 6000 N 98Jac2 !
#
        y[0]*y[1] * (-42304.49 + 7.45704*T) + 
        y[0]*y[1]*(y[0]-y[1]) * (-15611.66 + 9.11885*T) + 
#
# (ideal solution entropy)
#
        1.0*R*T*xlnx(y[0]) + 1.0*R*T*xlnx(y[1])
    )

#
# PHASE FCC_A1  %F  2 1   1 !
# CONST FCC_A1 : AG%,AL%,...,MG,...,NI,...,ZN,ZR : B,C,N,O,VA : !
#
# (Mg,Ni)1 (Va)1
#
MgNi_fcc = Phase()
MgNi_fcc.X = lambda y: y[1]
MgNi_fcc.G = lambda T, y: (
#
# (mechanical mixture)
#
# PARAM G(FCC_A1,MG:VA;0) 298.15 +GFCCMG#; 3000 N 91Din !
# PARAM G(FCC_A1,NI:VA;0) 298.15 +GHSERNI#; 6000 N 91Din !
#
        y[0] * (GFCCMG(T)) + 
        y[1] * (GHSERNI(T)) + 
#
# (regular solution interaction)
#
# PARAM G(FCC_A1,MG,NI:VA;0) 298.15 +80*T; 6000 N 98Jac2 !
#
        y[0]*y[1] * (80.0*T) + 
#
# (ideal solution entropy)
#
        1.0*R*T*xlnx(y[0]) + 1.0*R*T*xlnx(y[1])
    )
#
# (ignored magnetic contribution)
#
# PARAM TC(FCC_A1,NI:VA;0) 298.15 633; 3000 N 01Din !
# PARAM BMAGN(FCC_A1,NI:VA;0) 298.15 .52; 3000 N 01Din  !
#



#
# PHASE HCP_A3  %F  2 1   .5 !
# CONST HCP_A3  : AG,AL,...,MG%,...,NI,...,ZN,ZR% : B,C,N,O,VA : !
#
# (Mg,Ni)1 (Va)0.5
#
MgNi_hcp = Phase()
MgNi_hcp.X = lambda y: y[1]
MgNi_hcp.G = lambda T, y: (
#
# (mechanical mixture)
#
# PARAM G(HCP_A3,MG:VA;0) 298.15 +GHSERMG#; 3000 N 01Din !
# PARAM G(HCP_A3,NI:VA;0) 298.15 +1046+1.255*T+GHSERNI#; 3000 N 01Din !
#
        y[0] * (GHSERMG(T)) + 
        y[1] * (1046.0 + 1.255*T + GHSERNI(T)) + 
#
# (regular solution interaction)
#
#  PARAM G(HCP_A3,MG,NI:VA;0) 298.15 +80*T; 6000 N 98Jac2 !
#
        y[0]*y[1] * (80.0*T) + 
#
# (ideal solution entropy)
#
        1.0*R*T*xlnx(y[0]) + 1.0*R*T*xlnx(y[1])
    )
#
# (ignored magnetic contribution)
#
# PARAM TC(HCP_A3,NI:VA;0) 298.15 633; 6000 N 88Gui2 !
# PARAM BMAGN(HCP_A3,NI:VA;0) 298.15 .52; 6000 N 88Gui2 !
#



#
# PHASE C36_LAVES  %  2 2   1 !
# CONST C36_LAVES : CO%,CR%,MG, NI%,TA,TI, ZN,ZR : CO, CR, MG%,NI, TA,TI%,ZN,ZR% : !
#
# (Mg,Ni)2 (Mg,Ni)1
#
MgNi2 = Phase()
MgNi2.X = lambda y: (2*y[1] + y[3])/3
MgNi2.G = lambda T, y: (
#
# (mechanical mixture)
#
# PARAM G(C36_LAVES,MG:MG;0) 298.15 +15000+3*GHSERMG#; 6000 N 98Jac2 !
# PARAM G(C36_LAVES,NI:MG;0) 298.15 -74136+293.9216*T-54.35385*T*LN(T)
#    -.0333*T**2-99*T**(-1)+5.14203E-06*T**3; 6000 N 98Jac2 !
# PARAM G(C36_LAVES,MG:NI;0) 298.15 +104136-293.9216*T+54.35385*T*LN(T)
#    +.0333*T**2+99*T**(-1)-5.14203E-06*T**3; 6000 N 98Jac2 !
#
        y[0] * y[2] * (15000.0 + 3.0*GHSERMG(T)) + 
        y[1] * y[2] * (-74136.0 + 293.9216*T + -54.35385*T*np.log(T) + -0.0333*T**2 + -99.0*T**-1 + 5.14203e-06*T**3) + 
        y[0] * y[3] * (104136.0 + -293.9216*T + 54.35385*T*np.log(T) + 0.0333*T**2 + 99.0*T**-1 + -5.14203e-06*T**3) + 
#
# (reguar solution interation on each sublattice)
#
# PARAM G(C36_LAVES,MG,NI:MG;0) 298.15 50000; 6000 N 98Jac2 !
# PARAM G(C36_LAVES,MG,NI:NI;0) 298.15 50000; 6000 N 98Jac2 !
# PARAM G(C36_LAVES,MG:MG,NI;0) 298.15 50000; 6000 N 98Jac2 !
# PARAM G(C36_LAVES,NI:MG,NI;0) 298.15 50000; 6000 N 98Jac2 !
#
        y[0]*y[1] * y[2] * (50000.0) + 
        y[0]*y[1] * y[3] * (50000.0) + 
        y[0] * y[2]*y[3] * (50000.0) + 
        y[1] * y[2]*y[3] * (50000.0) + 
#
# (ideal solution entropy)
#
        2.0*R*T*xlnx(y[0]) + 2.0*R*T*xlnx(y[1]) + 1.0*R*T*xlnx(y[2]) + 1.0*R*T*xlnx(y[3])
    ) / 3



#
# PHASE MG2NI  %  2 2   1 !
# CONST MG2NI  : MG : NI : !
#
# (Mg)2 (Ni)1
#
Mg2Ni = Phase()
Mg2Ni.X = lambda y: np.array([1/3])
Mg2Ni.G = lambda T, y: (
#
# PARAM G(MG2NI,MG:NI;0) 298.15 -82211.2+571.0183*T-95.992*T*LN(T);
#    6000 N 98Jac2 !
#
        -82211.2 + 571.0183*T + -95.992*T*np.log(T)
    )/3
