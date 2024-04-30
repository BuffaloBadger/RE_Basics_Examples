"""Calculations for Reaction Engineering Basics Example 12.7.1"""

import reb_12_7_1
import numpy as np
import pandas as pd

# read in the results from the assignment
results_df = pd.read_csv('reb_12_7_1/python/results.csv')

# set the initial guess
init_guess = np.array([0.5*reb_12_7_1.nA_in, 0.5*reb_12_7_1.nA_in
                       , 0.1*reb_12_7_1.nA_in, 0.1*reb_12_7_1.nA_in
            , reb_12_7_1.T_in + 10.0])

# solve the design equations for a lower inlet temperature
reb_12_7_1.T_in = 325 # K
solution = reb_12_7_1.unknowns(init_guess)
# extract individual unknowns
nA = solution[0]
nD = solution[2]
nU = solution[3]
T = solution[4]

# calculate the other quantities of interest
fA = 100*(reb_12_7_1.nA_in - nA)/reb_12_7_1.nA_in
S_D_U = nD/nU

# add the new results
n_rows = len(results_df.index)
results_df.loc[n_rows] = ['Conversion(325 K)',f'{fA}','%']
results_df.loc[n_rows+1] = ['Selectivity(325 K)',f'{S_D_U}','mol D per mol U']
results_df.loc[n_rows+2] = ['T out (325)',f'{T}', 'K']

# solve the design equations for a higher inlet temperature
reb_12_7_1.T_in = 375 # K
solution = reb_12_7_1.unknowns(init_guess)
# extract individual unknowns
nA = solution[0]
nD = solution[2]
nU = solution[3]
T = solution[4]

# calculate the other quantities of interest
fA = 100*(reb_12_7_1.nA_in - nA)/reb_12_7_1.nA_in
S_D_U = nD/nU

# add the new results
n_rows = len(results_df.index)
results_df.loc[n_rows] = ['Conversion(375 K)',f'{fA}','%']
results_df.loc[n_rows+1] = ['Selectivity(375 K)',f'{S_D_U}','mol D per mol U']
results_df.loc[n_rows+2] = ['T out (375)',f'{T}', 'K']


# display the results
print(results_df)

# save the results
results_df.to_csv('reb_12_7_1/python/results.csv', index=False)

