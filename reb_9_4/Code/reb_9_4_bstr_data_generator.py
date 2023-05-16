import numpy as np
import pandas as pd
import scipy as sp
import rebutils as reb
import math
import random

# Constant inputs
V = 50.0 # ml

# kinetics parameters
Vm = 1.15E-4 # mmol/ml/min
Km = 0.0021 # mmol/ml

# Define rate expression
def rate(CS):
    reaction_rate = Vm*CS/(Km + CS)
    return reaction_rate

# Adjusted inputs
CS0 = np.array([15.0, 10.0, 5.0]) * 1E-3 # mmol/mL
t_reaction = np.linspace(5.0, 120.0, 24) # min

# Create empty dataframe for the results
df = pd.DataFrame(columns=["CS0", "t", "CP"])

# Calculate the responses
for CS_init in CS0:
    for time in t_reaction:
        # Define the mole balances
        def mole_balances(t,n): #n[0] = nS, n[1] = nP, n[2] = nH2O
            C_S = n[0]/V
            r = rate(C_S)
            ddt = np.array([-r*V, r*V, r*V])
            return ddt
                
        # Solve the mole balances
        t0 = 0
        n0 = np.array([CS_init*V, 0.0, 0.0])
        f_var = 0
        f_val = time
        soln = reb.solveIVODEs(t0, n0, f_var, f_val, mole_balances)

        # calculate the response
        CP = soln.y[1,-1]/V

        # add +/- 0.01 random "error"
        random_error = (2*random.random() - 1.0)*0.0004
        CP = CP + random_error

        # convert to mmol/L
        CS = 1000.0*CS_init
        CP = 1000.0*CP

        # round CP to 2 decimal places
        CP = round(CP,2)

        # append the result to the dataframe
        df.loc[len(df.index)] = [CS, time, CP]

# display the results
print("\n")
print(df)

# save the results
filename = "reb_9_4/Data/reb_9_4_data.csv"
print("\nSaving results to " + filename + "\n")
df.to_csv(filename,index=False)

filename = "../RE_Basics/Data/reb_9_4_data.csv"
print("\nSaving results to " + filename + "\n")
df.to_csv(filename,index=False)
