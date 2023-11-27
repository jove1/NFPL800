
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
#FUN GHSERAG  298.15   -7209.512+118.202013*T-23.8463314*T*LN(T)
#                        -.001790585*T**2-3.98587E-07*T**3-12011*T**(-1);
#            1234.93 Y -15095.252+190.266404*T-33.472*T*LN(T)
#                        +1.412E+29*T**(-9);
#            3000.00 N 91Din !
#
def GHSERAG(T):
    return np.where(T < 1234.93, -7209.512 + 118.202013*T + -23.8463314*T*np.log(T) + -0.001790585*T**2 + -3.98587e-07*T**3 + -12011.0*T**-1, 
                                -15095.252 + 190.266404*T + -33.472*T*np.log(T) + 1.412e+29*T**-9)
#
#FUN GHSERCU  298.15   -7770.458+130.485235*T-24.112392*T*LN(T)
#                        -.00265684*T**2+1.29223E-07*T**3+52478*T**(-1);
#            1357.77 Y -13542.026+183.803828*T-31.38*T*LN(T)
#                        +3.642E+29*T**(-9);
#            3200.00  N 91Din !
#
def GHSERCU(T):
    return np.where(T < 1357.77, -7770.458 + 130.485235*T + -24.112392*T*np.log(T) + -0.00265684*T**2 + 1.29223e-07*T**3 + 52478.0*T**-1, 
                                -13542.026 + 183.803828*T + -31.38*T*np.log(T) + 3.642e+29*T**-9)
#
#PHASE LIQUID:L %  1  1.0 !
#   CONSTITUENT LIQUID:L  :AG,BI,CU,PB,SB,SN : !
#
CuAg_liq = Phase()
CuAg_liq.X = lambda y: y
CuAg_liq.G = lambda T, y: (
#
# (mechanical mixture)
#
#PARA G(LIQUID,AG;0)  298.15   +11025.076-8.891021*T-1.034E-20*T**7
#                                +GHSERAG#;
#                    1235.08 Y +11508.141-9.301747*T-1.412E+29*T**(-9)
#                                +GHSERAG#;
#                    3000.00  N 91Din !
        y[0] * np.where(T < 1235.08, 11025.076 + -8.891021*T + -1.034e-20*T**7 + GHSERAG(T), 
                                     11508.141 + -9.301747*T + -1.412e+29*T**-9 + GHSERAG(T)) + 
#
#PARA G(LIQUID,CU;0)  298.15   +12964.736-9.511904*T-5.849E-21*T**7
#                                 +GHSERCU#;
#                    1358.02 Y +13495.481-9.922344*T-3.642E+29*T**(-9)
#                                 +GHSERCU#;  3200.00  N 91Din !
#
        y[1] * np.where(T < 1358.02, 12964.736 + -9.511904*T + -5.849e-21*T**7 + GHSERCU(T), 
                                     13495.481 + -9.922344*T + -3.642e+29*T**-9 + GHSERCU(T)) + 
#
# (subregular solution interaction)
#
#PARA L(LIQUID,AG,CU;0)  298.15  +17323.40-4.46819*T;,,N 86Hay !
#PARA L(LIQUID,AG,CU;1)  298.15  +1654.38-2.35285*T;,,N 86Hay !
#
        y[0]*y[1] * (17323.4 + -4.46819*T) + 
        y[0]*y[1]*(y[0]-y[1]) * (1654.38 + -2.35285*T) + 
#
# (ideal solution entropy)
# 
        1.0*R*T*xlnx(y[0]) + 1.0*R*T*xlnx(y[1])
    )

#
#PHASE FCC_A1 %A  2  1.0 1.0 !
#   CONSTITUENT FCC_A1  :AG%,BI,CU%,PB%,SB,SN :VA : !
#
CuAg_fcc = Phase()
CuAg_fcc.X = lambda y: y
CuAg_fcc.G = lambda T, y: (
#
# (mechanical mixture)
#
#PARA G(FCC_A1,AG:VA;0)  298.15  +GHSERAG#;,,N 0 !
#PARA G(FCC_A1,CU:VA;0)  298.15  +GHSERCU#;,,N 0 !
#
        y[0] * (GHSERAG(T)) + 
        y[1] * (GHSERCU(T)) + 
#
# (subregular solution interaction)
# 
#PARA L(FCC_A1,AG,CU:VA;0)  298.15  +36061.88-10.44288*T;,,N 86Hay !
#PARA L(FCC_A1,AG,CU:VA;1)  298.15  -4310.12;,,N 86Hay !
#
        y[0]*y[1] * (36061.88 + -10.44288*T) + 
        y[0]*y[1]*(y[0]-y[1]) * (-4310.12) + 
#
# (ideal solution entropy)
#
        1.0*R*T*xlnx(y[0]) + 1.0*R*T*xlnx(y[1])
    )
