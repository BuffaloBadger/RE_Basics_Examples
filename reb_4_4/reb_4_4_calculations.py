import rebutils as reb
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math

# given or known
R = 8.314E-3 # kJ/mol/K

# read the data file
df = pd.read_csv('./reb_4_4/reb_4_4_data.csv')
df.columns = ['i','T','k']
T = df['T'].to_numpy()
k = df['k'].to_numpy()

# create adjusted input vector, eqn 4
x = T + 273.15

# create experimental response vector, eqn 5
y_expt = k

# response model function
def calc_y(T_K,log_k0,E):
    return 10**log_k0*math.e**(-E/R/T_K) # eqn 3

# guess the values of the parameters
guess = np.array([0.,20.])

# fit the function to the data
use_relative_errors = False
beta, beta_ci, r_squared = reb.fitNLSR(guess,x,y_expt,calc_y,
    use_relative_errors)

# extract the results
k0 = 10**beta[0] # eqn 6
k0_lower = 10**beta_ci[0,0] # eqn 6
k0_upper = 10**beta_ci[0,1] #eqn 6
E = beta[1]
E_lower = beta_ci[1,0]
E_upper = beta_ci[1,1]

# report the results
print(' ')
print(f'k0: {k0:.3g} [{k0_lower:.3g}, {k0_upper:.3g}] (95%CI) L/mol/min')
print(f'E: {E:.3g} [{E_lower:.3g}, {E_upper:.3g}] (95%CI) kJ/mol')
print(f'R-squared: {r_squared:.3g}')
print(' ')

# save the results to a .csv file
data = [['k0', f'{k0:.3g}', 'L/mol/min'],
    ['k0_lower_limit', f'{k0_lower:.3g}', 'L/mol/min'],
    ['k0_upper_limit', f'{k0_upper:.3g}', 'L/mol/min'],
    ['E', f'{E:.3g}', 'kJ/mol'],
    ['E_lower_limit', f'{E_lower:.3g}', 'kJ/mol'],
    ['E_upper_limit', f'{E_upper:.3g}', 'kJ/mol'],
    ['R_squared', f'{r_squared:.3g}', '']]
result = pd.DataFrame(data, columns=['item','value','units'])
if use_relative_errors:
    result.to_csv("./reb_4_4/reb_4_4_Python_results_rel.csv", index=False)
else:
    result.to_csv("./reb_4_4/reb_4_4_Python_results_abs.csv", index=False)

# generate the parity plot
y_pred = calc_y(x,beta[0],beta[1]) # eqn 7
parity_range = np.array([0.9*np.min(y_expt), 1.1*np.max(y_expt)])
plt.figure(1)
plt.loglog(y_expt,y_pred,color='r',marker='o', ls='none', label='Actual Model')
plt.loglog(parity_range,parity_range,color='k', label='Perfect Model')
plt.xlabel('Measured k (L mol$^{-1}$ min$^{-1}$)')
plt.ylabel('Predicted k (L mol$^{-1}$ min$^{-1}$)')
plt.legend()
# save and show the figure
if use_relative_errors:
    plt.savefig('reb_4_4/reb_4_4_Python_parity_plot_rel.png')
else:
    plt.savefig('reb_4_4/reb_4_4_Python_parity_plot_abs.png')
plt.show()

# generate the residuals plot
residual = y_expt - y_pred # eqn 8
plt.figure(2)
plt.plot(x,residual,color='r',marker='o',ls='none')
plt.axhline(y=0, color = 'k')
plt.xlabel('T (K)')
plt.ylabel('Residual (L mol$^{-1}$ min$^{-1}$)')
# save and show the figure
if use_relative_errors:
    plt.savefig('reb_4_4/reb_4_4_Python_residual_plot_rel.png')
else:
    plt.savefig('reb_4_4/reb_4_4_Python_residual_plot_abs.png')
plt.show()
