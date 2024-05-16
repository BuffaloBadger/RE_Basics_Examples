"""Calculations for Reaction Engineering Basics Example 13.7.6 Discussion"""

# import libraries
import pandas as pd
import reb_13_7_6

reb_13_7_6.eta = 1.0
reb_13_7_6.Pconv = 0.0


# solve the reactor design equations
z, nA, nZ, P = reb_13_7_6.profiles()

# calculate the other quantities of interest
L = z[-1]
P_out = P[-1]

# read in the results from the assignment
results_df = pd.read_csv('reb_13_7_6/python/results.csv')

# add the new results
n_rows = len(results_df.index)
results_df.loc[n_rows] = ['effectiveness factor',1.0,'']
results_df.loc[n_rows+1] = ['reactor length',L,'cm']

# display the results
print(' ')
print(results_df)

# save the results
results_df.to_csv('reb_13_7_6/python/results.csv', index=False)
