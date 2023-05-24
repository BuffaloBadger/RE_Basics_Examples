import numpy as np
import pandas as pd
import scipy as sp
import rebutils as reb
import math
import random

# Filenames for data
reb_examples_filename = 'reb_9_1/Data/reb_9_1_data.csv'
reb_filename = '../RE_Basics/Data/reb_9_1_data.csv'

# Constant inputs
V = 1.0 # L
R = 8.314E-3 # kJ/mol/K

# Define rate expression
def rate(conc_A, temp):
    k0 = 4.23E8
    E = 68.0
    reaction_rate = k0*math.e**(-E/R/temp)*conc_A
    return reaction_rate

# Adjusted inputs
T_expt = np.array([65, 73, 82, 90]) + 273.15
CA0_expt = np.array([0.5, 1.0, 1.5])
t_reaction = np.array([5.0, 10.0, 15.0, 20.0, 25.0, 30.0])

# Create empty dataframe for the results
df = pd.DataFrame(columns=["Experiment","T", "CA0", "t", "CA"])

# Calculate the responses
expt_number = 0
for CA0 in CA0_expt:
    for T in T_expt:
        expt_number += 1
        for time in t_reaction:
            # Define the mole balances
            def mole_balances(t,C):
                ddt = [-rate(C,T)]
                return ddt
            
            # Solve the mole balances
            t0 = 0
            f_var = 0
            soln = reb.solveIVODEs(t0, [CA0], f_var, time, mole_balances)

            # calculate the response
            CA = soln.y[0,-1]

            # add +/- 0.01 M random "error"
            random_error = (2*random.random() - 1.0)*0.01
            CA = CA + random_error

            # round to 2 decimal places
            CA = round(CA,2)

            # append the result to the dataframe
            df.loc[len(df.index)] = [expt_number, T - 273.15, CA0, time, CA]

# display the results
print("\n")
print(df)

# save the results

print("\nSaving results to " + reb_examples_filename + "\n")
df.to_csv(reb_examples_filename,index=False)

print("\nSaving results to " + reb_filename + "\n")
df.to_csv(reb_filename,index=False)

