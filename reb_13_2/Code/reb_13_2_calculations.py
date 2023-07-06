import numpy as np
import pandas as pd
from scipy import optimize
import response_function as resp

# Path to directory for saving files
filepath = './reb_13_2/Results/'

# Given and known constants in consistent units
specified_fA = 0.45
TexIn = 40.0 + 273.15 # K

# Define the trivial residual equation
def residuals(x):
    _, _, fAf, _ = resp.response(x)
    return np.array([specified_fA - fAf])

guess = np.array([TexIn])
soln = optimize.root(residuals, guess)
if not(soln.success):
    print("A solution was NOT obtained:")
    print(soln.message)
T0 = soln.x[0]

# Calculate the response
Tf, Texf, fAf, sel = resp.response(np.array([T0]))

# Show and save the results
data = [['$T_0$', f'{(T0 - 273.15):.3g}', '°C'],
    ['$T_f$', f'{(Tf - 273.15):.3g}', '°C'],
    ['$T_{ex,out}$', f'{(Texf - 273.15):.3g}', '°C'],
    ['$f_A$', f'{100.0*fAf:.3g}', '%'],
    ['$S_{X/Z}$', f'{sel:.3g}', ' mol X per mol Z']]
result = pd.DataFrame(data, columns=['item','value','units'])
print(result)
result.to_csv(filepath + "reb_13_2_results.csv", index=False)