import rebutils as reb
import numpy as np
import pandas as pd
from scipy import optimize

# Path to directory for saving files
filepath = './reb_12_2/python/'

# Given and known constants in consistent units
V = 1.0E4 # cm^3
V_ex = 1.4E3 # cm^3
U = 138.0 / 929.0 # cal cm^-2 min^-1 K^-1
A = 1.2E3 # cm^2
TexIn = 40.0 + 273.15 # K
mfr_e = 100.0 # g min^-1
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
fA = 0.45
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
        dnbdt = -r1*V
        dnXdt = r1*V
        dnYdt = r1*V
        dnZdt = r2*V
        dTdt = (Q-V*(r1*dH1 + r2*dH2))/Cp/rho/V
        dTexdt = (-Q - mfr_e*Cpex*(Tex-TexIn))/rho_ex/V_ex/Cpex
        ddt = np.array([dnAdt, dnbdt, dnXdt, dnYdt, dnZdt, dTdt, dTexdt])
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
    Tf = soln.y[5,-1]
    Texf = soln.y[6,-1]

    # Calculate the response
    fAf = (CA0*V - nAf)/CA0/V
    sel = nXf/nZf

    # Return the response
    return Tf, Texf, fAf, sel

# Calculate missing input to response function
def residuals(x):
    _, _, fAf, _ = response(x)
    return np.array([fA - fAf])

guess = np.array([TexIn])
soln = optimize.root(residuals, guess)
if not(soln.success):
    print("A solution was NOT obtained:")
    print(soln.message)
T0 = soln.x[0]

# Calculate the response
Tf, Texf, fAf, sel = response(np.array([T0]))

# Show and save the results
data = [['T0', f'{(T0 - 273.15):.3g}', '°C'],
    ['Tf', f'{(Tf - 273.15):.3g}', '°C'],
    ['Texf', f'{(Texf - 273.15):.3g}', '°C'],
    ['fA', f'{100.0*fAf:.3g}', '%'],
    ['S_X_per_Z', f'{sel:.3g}', ' mol X per mol Z']]
result = pd.DataFrame(data, columns=['item','value','units'])
print(result)
result.to_csv(filepath + "reb_12_2_results.csv", index=False)