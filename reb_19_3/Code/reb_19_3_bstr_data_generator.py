import numpy as np
import pandas as pd
import rebutils as reb
import random

# Constant inputs
Rpv = 82.06 # cc atm/mol/K
V = 100 # cc
P0 = 6.0 # atm

# kinetics parameters
k = 1.62E-6 # mol/cc/atm^1.5/min
T = 275.0 + 273.15

# Define rate expression
def rate(P_A, P_B, temp):
    if P_A <= 0.0:
        reaction_rate = 0.0
    elif P_B <= 0.0:
        reaction_rate = 0.0
    else:
        reaction_rate = k*P_A*np.sqrt(P_B)
    return reaction_rate

# Adjusted inputs
PAin = np.array([2.0, 3.0, 4.0])
t_reaction = np.arange(1.0, 25, 1.0)

# Create empty dataframe for the results
df = pd.DataFrame(columns=["PA0", "t", "P"])

# Calculate the responses
for PA0 in PAin:
    PB0 = P0 - PA0
    for time in t_reaction:
        # Define the mole balances
        def mole_balances(t,n):
            PA = n[0]*Rpv*T/V
            PB = n[1]*Rpv*T/V
            r = rate(PA,PB,T)
            ddt = np.array([-r*V, -r*V, r*V])
            return ddt
        
        # Solve the mole balances
        t0 = 0
        n0 = [PA0*V/Rpv/T, PB0*V/Rpv/T, 0.0]
        f_var = 0
        f_val = time
        soln = reb.solveIVODEs(t0, n0, f_var, f_val, mole_balances)

        # calculate the response
        nA = soln.y[0,-1]
        nB = soln.y[1,-1]
        nZ = soln.y[2,-1]
        P = (nA + nB + nZ)*Rpv*T/V

        # add +/- 0.05 atm random "error"
        random_error = (2*random.random() - 1.0)*0.05
        P = P + random_error

        # round to 2 decimal places
        P = round(P,2)

        # append the result to the dataframe
        df.loc[len(df.index)] = [PA0, time, P]

# display the results
print("\n")
print(df)

# save the results
filename = "./reb_19_3/Data/reb_19_3_data.csv"
print("\nSaving results to " + filename + "\n")
df.to_csv(filename,index=False)

filename = "../RE_Basics/Data/reb_19_3_data.csv"
print("\nSaving results to " + filename + "\n")
df.to_csv(filename,index=False)
