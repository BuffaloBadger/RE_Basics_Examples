"""Calculations for Reaction Engineering Basics Example 19.5.1"""

# import libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from score_utils import Arrhenius_parameters

# Ideal gas constant
R = 8.314E-3 # kJ/mol/K

# read rate coefficients from predicted responses function
df = pd.read_csv('reb_19_5_1/python/reb_19_5_1_resp_fcn_params.csv')
    # columns: 'T', 'k', 'k_ll', 'k_ul', 'R_sq'
T = df['T'].to_numpy() + 273.15
k = df['k'].to_numpy()

# fit the Arrhenius expression to the data
k0, k0_ci, E, E_ci, r_squared = Arrhenius_parameters(k,T,R)

# report the results
print(' ')
print(f'k0: {k0:.3g} L/mol/min, 95% CI [{k0_ci[0]:.3g}, {k0_ci[1]:.3g}]')
print(f'E: {E:.3g} kJ/mol, 95% CI [{E_ci[0]:.3g}, {E_ci[1]:.3g}]')
print(f'R-squared: {r_squared:.3g}')
print(' ')

# save the results to a .csv file
data = [['k0', f'{k0:.3g}', 'min^-1^'],
    ['k0_lower_limit', f'{k0_ci[0]:.3g}', 'min^-1^'],
    ['k0_upper_limit', f'{k0_ci[1]:.3g}', 'min^-1^'],
    ['E', f'{E:.3g}', 'kJ mol^-1^'],
    ['E_lower_limit', f'{E_ci[0]:.3g}', 'kJ mol^-1^'],
    ['E_upper_limit', f'{E_ci[1]:.3g}', 'kJ mol^-1^'],
    ['R_squared', f'{r_squared:.3g}', '']]
result = pd.DataFrame(data, columns=['item','value','units'])
result.to_csv("reb_19_5_1/python/reb_19_5_1_Arrhenius_resp_fcn.csv", 
              index=False)

# create, show, and save an Arrhenius plot
y_pred = k0*np.exp(-E/R/T)
plt.figure()
plt.semilogy(1/T,k,color='r',marker='o', ls='none')
plt.semilogy(1/T,y_pred,color='k')
plt.xlabel('T$^{-1}$ (K$^{-1}$)')
plt.ylabel('k (min$^{-1}$)')
plt.xticks(rotation=25)
plt.tight_layout()
plt.savefig('reb_19_5_1/python/reb_19_5_1_Arrhenius_resp_fcn.png')
plt.show()

# read rate coefficients from differential analysis
df = pd.read_csv('reb_19_5_1/python/reb_19_5_1_diff_params.csv')
    # columns: 'T', 'k', 'k_ll', 'k_ul', 'R_sq'
T = df['T'].to_numpy() + 273.15
k = df['k'].to_numpy()

# fit the Arrhenius expression to the data
k0, k0_ci, E, E_ci, r_squared = Arrhenius_parameters(k,T,R)

# report the results
print(' ')
print(f'k0: {k0:.3g} L/mol/min, 95% CI [{k0_ci[0]:.3g}, {k0_ci[1]:.3g}]')
print(f'E: {E:.3g} kJ/mol, 95% CI [{E_ci[0]:.3g}, {E_ci[1]:.3g}]')
print(f'R-squared: {r_squared:.3g}')
print(' ')

# save the results to a .csv file
data = [['k0', f'{k0:.3g}', 'min^-1^'],
    ['k0_lower_limit', f'{k0_ci[0]:.3g}', 'min^-1^'],
    ['k0_upper_limit', f'{k0_ci[1]:.3g}', 'min^-1^'],
    ['E', f'{E:.3g}', 'kJ mol^-1^'],
    ['E_lower_limit', f'{E_ci[0]:.3g}', 'kJ mol^-1^'],
    ['E_upper_limit', f'{E_ci[1]:.3g}', 'kJ mol^-1^'],
    ['R_squared', f'{r_squared:.3g}', '']]
result = pd.DataFrame(data, columns=['item','value','units'])
result.to_csv("reb_19_5_1/python/reb_19_5_1_Arrhenius_diff.csv", 
              index=False)

# create, show, and save an Arrhenius plot
y_pred = k0*np.exp(-E/R/T)
plt.figure()
plt.semilogy(1/T,k,color='r',marker='o', ls='none')
plt.semilogy(1/T,y_pred,color='k')
plt.xlabel('T$^{-1}$ (K$^{-1}$)')
plt.ylabel('k (min$^{-1}$)')
plt.xticks(rotation=25)
plt.tight_layout()
plt.savefig('reb_19_5_1/python/reb_19_5_1_Arrhenius_diff.png')
plt.show()