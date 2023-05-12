import rebutils as reb
import numpy as np
import pandas as pd
from scipy import optimize

# Path to directory for saving files
filepath = './reb_12_2/python/'

Tinit = 65 + 273.15;

# Given and known constants in consistent units
V = 1.0E4 # cm^3
V_ex = 1.4E3 # cm^3
U = 138.0 / 929.0 # cal cm^-2 min^-1 K^-1
A = 1.2E3 # cm^2
TexIn = 40.0 + 273.15 # K
mfr_e = 100 # g min^-1
rho = 1.0 # g cm^-3
rho_ex = rho
Cp = 1.0 # cal g^-1 K^-1
Cpex = Cp
CA0 = 5.0/1000.0 # mol cm^-3
CB0 = 7.0/1000.0 # mol cm^-3
dH1 = -1.67E4 # cal mol^-1
dH2 = -1.43E4 # cal mol^-1
k01 = 9.74E12 # cm^3 mol^-1 min^-1
E1 = 2.01E4 # cal mol^-1
k02 = 2.38E13 # min^-1
E2 = 2.53E4 # cal mol^-1
fAf = 0.5
t_rxn = 30.0 # min
R = 1.987 # cal mol^-1 K^-1

# Define the response function
def response(x):
    T0 = x[0]
    # Define the reactor model
    def reactor_model(current_ind, current_dep):
        # Evaluate the derivatives at the current values of the variables
        nA = current_dep[0]
        nB = current_dep[1]
        nX = current_dep[2]
        nY = current_dep[3]
        nZ = current_dep[4]
        T = current_dep[5]
        Tex = current_dep[6]
        CA = nA/V
        CB = nB/V
        r1 = k01*np.exp(-E1/R/T)*CA*CB
        r2 = k02*np.exp(-E2/R/T)*CA
        Q = U*A*(Tex - T)
        dnAdt = (-r1 - r2)*V
        dnBdt = -r1*V
        dnXdt = r1*V
        dnYdt = r1*V
        dnZdt = r2*V
        dTdt = (Q-V*(r1*dH1 + r2*dH2))/Cp/rho/V
        dTexdt = (-Q - mfr_e*Cpex*(Tex-TexIn))/rho_ex/V_ex/Cpex
        ddt = np.array([dnAdt, dnBdt, dnXdt, dnYdt, dnZdt, dTdt, dTexdt])
        # Return the derivatives/residuals
        return ddt
    
    # Solve the reactor model equations
    t0 = 0.0
    dep0 = np.array([CA0*V, CB0*V, 0.0, 0.0, 0.0, T0, TexIn])
    f_var = 0
    f_val = t_rxn
    soln = reb.solveIVODEs(t0, dep0, f_var, f_val, reactor_model, \
            odes_are_stiff=False)
    if not(soln.success):
        print("A solution was NOT obtained:")
        print(soln.message)

    # Extract the results
    nAf = soln.y[0,-1]
    nXf = soln.y[2,-1]
    nZf = soln.y[4,-1]
    Tf = soln.y[5,-1] - 273.15
    Texf = soln.y[6,-1] - 273.15

    # Calculate the response
    fAf = (CA0*V - nAf)/CA0/V*100.0
    sel = nXf/nZf

    print('\nTf: %.1f, Texf: %.1f, fA: %.1f, sel: %.3g'\
          %(Tf, Texf, fAf, sel))

    # Return the response
    return Tf, fAf, sel

Tf, fA, sel = response([Tinit])