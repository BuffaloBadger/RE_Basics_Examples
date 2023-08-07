import numpy as np
import pandas as pd
import scipy as sp
import rebutils as reb
import math
import random

# Constant inputs
D = 1.0 # cm
L = 10.0 # cm
P = 1.0 # atm
T = 1500.0 # K
R = 8.314E-3 # kJ/mol/K
mol_per_cc_stp = 1.0/1000.0/22.4 # mol/cc
Rpv = 82.06 # cc atm/mol/K

# kinetics parameters
k = 1.5E-3 # mol cm-3 min-1 atm-1

# Define rate expression
def rate(P_A):
    return k*P_A

# Adjusted inputs
VFR_Ain = np.array([30, 60, 90, 120, 150, 180, 210, 240, 270, 300])
VFR_Yin = np.array([0, 50, 100, 150])
VFR_Zin = np.array([0, 50, 100, 150])

# Create empty dataframe for the results
df = pd.DataFrame(columns=["A_in", "Y_in", "Z_in", "fA"])

# Calculate the responses
for VdotAin in VFR_Ain:
    nAin = VdotAin * mol_per_cc_stp
    for VdotYin in VFR_Yin:
        nYin = VdotYin * mol_per_cc_stp
        for VdotZin in VFR_Zin:
            nZin = VdotZin * mol_per_cc_stp
            
            # Define the mole balances
            def mole_balances(z,n):
                n_total = np.sum(n)
                PA = n[0]/n_total*P
                PY = n[1]/n_total*P
                PZ = n[2]/n_total*P
                r = rate(PA)
                ddz = np.array([-r, r, r]) *math.pi*D**2/4
                return ddz
            
            # Solve the mole balances
            z0 = 0
            n0 = np.array([nAin, nYin, nZin])
            f_var = 0
            f_val = L
            soln = reb.solveIVODEs(z0, n0, f_var, f_val, mole_balances)

            # calculate the response
            nA = soln.y[0,-1]
            fA = (nAin - nA)/nAin*100.0

            # add +/- 0.01 random "error"
            random_error = (2*random.random() - 1.0)*2.0
            fA += random_error

            # round to 3 decimal places
            fA = round(fA,1)

            # append the result to the dataframe
            df.loc[len(df.index)] = [VdotAin, VdotYin, VdotZin, fA]

# display the results
print("\n")
print(df)

# save the results
filename = "reb_10_1/Data/reb_10_1_data.csv"
print("\nSaving results to " + filename + "\n")
df.to_csv(filename,index=False)

filename = "../RE_Basics/Data/reb_10_1_data.csv"
print("\nSaving results to " + filename + "\n")
df.to_csv(filename,index=False)
