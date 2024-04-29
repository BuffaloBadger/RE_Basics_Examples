"""Calculations for Reaction Engineering Basics Example I.6 Discussion"""

# import libraries
import numpy as np
import pandas as pd
import reb_I_6

# read in the results from the assignment
results_df = pd.read_csv('reb_I_6/python/results.csv')
sol = results_df['value'].to_numpy()
sol[2] = sol[2] + 273.15

# calculate the residuals
residuals = reb_I_6.residuals(sol)
sol[2] = sol[2] - 273.15

# add the new results
n_rows = len(results_df.index)
results_df.loc[n_rows] = ['residual 1',f'{residuals[0]}',' ']
results_df.loc[n_rows+1] = ['residual 2',f'{residuals[1]}',' ']
results_df.loc[n_rows+2] = ['residual 3',f'{residuals[2]}',' ']

# display the results
print(results_df)

# save the results
results_df.to_csv('reb_I_6/python/results.csv', index=False)
