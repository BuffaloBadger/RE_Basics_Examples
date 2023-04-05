import numpy as np
import pandas as pd
import scipy as sp
import rebutils as reb
import math
import random

# Constant inputs
R = 8.314E-3 # kJ/mol/K
Rpv = 82.06 # cc atm/mol/K
V = 100 # cc
PA0 = 4.0 # atm

# kinetics parameters
k = 9.62E-7 # mol/cc/atm/s
E = 87.6 # kJ/mol
T_avg = 275 + 273.15
k0 = k*math.e**(E/R/T_avg)
print(k0)

# Define rate expression
def rate(P_A, temp):
    reaction_rate = k0*math.e**(-E/R/temp)*P_A
    return reaction_rate

# Adjusted inputs
T_expt = np.array([250, 275, 300]) + 273.15
t_reaction = np.arange(1.0, 25, 1.0)

# Create empty dataframe for the results
df = pd.DataFrame(columns=["T", "t", "P"])

# Calculate the responses
for T in T_expt:
    for time in t_reaction:
        # Define the mole balances
        def mole_balances(t,n):
            PA = n[0]*Rpv*T/V
            r = rate(PA,T)
            ddt = np.array([-2*r*V, r*V])
            return ddt
        
        # Solve the mole balances
        t0 = 0
        n0 = [PA0*V/Rpv/T, 0]
        f_var = 0
        f_val = time
        soln = reb.solveIVODEs(t0, n0, f_var, f_val, mole_balances)

        # calculate the response
        nA = soln.y[0,-1]
        nA2 = soln.y[1,-1]
        P = (nA + nA2)*Rpv*T/V

        # add +/- 0.05 atm random "error"
        random_error = (2*random.random() - 1.0)*0.05
        P = P + random_error

        # round to 2 decimal places
        P = round(P,2)

        # append the result to the dataframe
        df.loc[len(df.index)] = [T - 273.15, time, P]

# display the results
print("\n")
print(df)

# save the results
filename = "reb_8_3/reb_8_3_data.csv"
print("\nSaving results to " + filename + "\n")
#df.to_csv(filename,index=False)
