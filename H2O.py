#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.rc("axes", grid=True)

#$
#$ PSUB - TC Public Substances Database (Version 1.1)
#$ ****   *******************************************
#
#$------------------------------------------------------------------------
#$ELEM NAME STABLE_ELEMENT_REF         ATOMIC MASS  H298-H0     S298 !
#$------------------------------------------------------------------------
# ELEMENT H    1/2_MOLE_H2(G)          1.0079E+00  0.0000E+00  1.5603E+01!
# ELEMENT O    1/2_MOLE_O2(G)          1.5999E+01  0.0000E+00  2.4502E+01!
#
# SPECIES H2O1                        H2O1!
#
#$------------------------------------------------------------------------
#$PHASE  NAME:TYPE  MARKCODE  #SUBL  SITES_IN_EACH_SUBL. !
#$------------------------------------------------------------------------
#$PHASE GAS:G %  1  1.0  !
# PHASE GAS:G %  1  1.0  
#  > Gaseous Mixture, using the ideal gas model !
# CONST GAS:G 
#  : H,H2, O,O2,O3, H1O1,H1O2,H2O1,H2O2, ...
#  : !
# PARAM G(GAS,H2O1;0)  298.15 +F10378T#+R#*T*LN(1E-05*P);   
#    6000.00 N REF1!
#
#'F10378T    2.98140E+02  -250423.434+4.45470381*T-28.40916*T*LN(T)             
#2     -.00623741*T**2-6.01526167E-08*T**3-64163.45*T**(-1);  1.10000E+03  Y    
#3      -256145.88+30.1894688*T-31.43044*T*LN(T)-.007055445*T**2                
#4     +3.05535833E-07*T**3+1246309.5*T**(-1);  2.80000E+03  Y                  
#5      -268423.418+116.690198*T-42.96842*T*LN(T)-.003069987*T**2               
#6     +6.97594167E-08*T**3+2458230.5*T**(-1);  8.40000E+03  Y                  
#7      -489068.882+553.259882*T-92.4077*T*LN(T)+.0016703495*T**2               
#8     -1.32333233E-08*T**3+1.765625E+08*T**(-1);  1.80000E+04  Y               
#9      -165728.771+239.645644*T-59.77872*T*LN(T)+2.213599E-04*T**2             
#:     -1.2921095E-09*T**3-4.1931655E+08*T**(-1);  2.00000E+04  N !             
#
#$PHASE H2O1_L  %  1  1.0  !
# PHASE H2O_L  %  1  1.0  !
# CONST H2O_L  :H2O1 :  !
# PARAM G(H2O_L,H2O1;0)  298.15 +F10373T#;     6000.00 N REF1!
#
#*F10373T    2.98140E+02  -332319.672+1078.59563*T-186.8669*T*LN(T)             
#2     +.2320948*T**2-9.14296167E-05*T**3+978019*T**(-1);  5.00000E+02  Y       
#3      -62418.8793-3288.18729*T+495.1304*T*LN(T)-.504926*T**2                  
#4     +4.917665E-05*T**3-18523425*T**(-1);  5.40000E+02  Y                     
#5      -8528143.9+142414.45*T-22596.19*T*LN(T)+27.48508*T**2                   
#6     -.00631160667*T**3+5.63356E+08*T**(-1);  6.00000E+02  Y                  
#7      -331037.282+741.178606*T-117.41*T*LN(T);  6.01000E+02  N !              
#

from common import *

G_gas = FT(298.14, { C: -250423.434, T: 4.45470381, TLNT: -28.40916, T2: -.00623741, T3: -6.01526167E-08, Tm1: -64163.45}, 1100)
G_liq = FT(298.14, { C: -332319.672, T: 1078.59563, TLNT: -186.8669, T2: +.2320948, T3: -9.14296167E-05, Tm1: +978019}, 500)

# CRC Handbook
H_liq_CRC_298, G_liq_CRC_298, S_liq_CRC_298, cp_liq_CRC = -285.8e3, -237.1e3, 70.0, 75.3
H_gas_CRC_298, G_gas_CRC_298, S_gas_CRC_298, cp_gas_CRC = -241.8e3, -228.6e3, 188.8, 33.6

T0 = 273.15
T298 = 298.15
T = T0 + np.linspace(20, 200)

S_liq_CRC = lambda T: NotImplementedError # ENTROPY FROM CRC Handbook data
S_gas_CRC = lambda T: NotImplementedError # ENTROPY FROM CRC Handbook data

plt.figure()
plt.title("S(T)")
plt.plot(T-T0, NotImplementedError) # ENTROPY FROM G_liq
plt.plot(T-T0, S_liq_CRC(T), "--" )
plt.plot(T-T0, NotImplementedError) # ENTROPY FROM G_gas
plt.plot(T-T0, S_gas_CRC(T), "--" )

H_liq_CRC = lambda T: NotImplementedError # ENTHALPY FROM CRC Handbook data
H_gas_CRC = lambda T: NotImplementedError # ENTHALPY FROM CRC Handbook data

plt.figure()
plt.title("H(T)")
plt.plot(T-T0, NotImplementedError) # ENTHALPY FROM G_liq
plt.plot(T-T0, H_liq_CRC(T), "--" )
plt.plot(T-T0, NotImplementedError) # ENTHALPY FROM G_gas
plt.plot(T-T0, H_gas_CRC(T), "--" )


G_liq_CRC = lambda T: NotImplementedError # GIBBS ENERGY FROM CRC Handbook data
G_gas_CRC = lambda T: NotImplementedError # GIBBS ENERGY FROM CRC Handbook data

plt.figure()
plt.title("G(T)")
plt.plot(T-T0, G_liq(T) )
plt.plot(T-T0, G_liq_CRC(T), "--")
plt.plot(T-T0, G_gas(T) )
plt.plot(T-T0, G_gas_CRC(T), "--")


print("S298_liq", S_liq_CRC_298, NotImplementedError) # FILL IN EXPRESSION FROM G_liq
print("S298_gas", S_gas_CRC_298, NotImplementedError) # FILL IN EXPRESSION FROM G_gas

print("H298_liq", H_liq_CRC_298, NotImplementedError) # FILL IN EXPRESSION FROM G_liq
print("H298_gas", H_gas_CRC_298, NotImplementedError) # FILL IN EXPRESSION FROM G_gas

print("cp_liq", cp_liq_CRC, NotImplementedError) # FILL IN EXPRESSION FROM G_liq
print("cp_gas", cp_gas_CRC, NotImplementedError) # FILL IN EXPRESSION FROM G_gas
print()

print("G298_liq", G_liq_CRC_298, NotImplementedError) # CALCULATE USING OTHER CRC Handbook data
print("G298_gas", G_gas_CRC_298, NotImplementedError) # CALCULATE USING OTHER CRC Handbook data
print()

# CALCULATE BOILING POINT AND LATENT HEAT USING Thermocalc GIBBS FUNCTIONS G_liq and G_gas
from scipy.optimize import root_scalar
Tboil = root_scalar( lambda T: NotImplementedError, # WRITE YOUR OWN FUNCTION
        bracket=(T0+20,T0+200)).root

print("Tboil = 100 C", Tboil - T0 )
print("Hvap = 40.65 kJ/mol", NotImplementedError ) # CALCULATE FROM G_liq AND G_gas
print()

# CALCULATE BOILING POINT AND LATENT HEAT USING CRC Handbook GIBBS FUNCTIONS G_liq_CRC and G_gas_CRC
Tboil = root_scalar( lambda T: NotImplementedError, # WRITE YOUR OWN FUNCTION
        bracket=(T0+20,T0+200)).root


print("Tboil = 100 C", Tboil - T0 )
print("Hvap = 40.65 kJ/mol", NotImplementedError ) # CALCULATE FROM G_liq2 AND G_gas2
print()

plt.show()
