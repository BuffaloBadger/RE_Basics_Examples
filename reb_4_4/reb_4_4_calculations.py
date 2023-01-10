import rebutils as reb
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math

# given or known
R = 8.314E-3 # kJ/mol/K

# read the data file
df = pd.read_csv('./reb_4_4/reb_4_4_data.csv')
df.columns = ['T_C','k']
T_C = df['T_C'].to_numpy()
k = df['k'].to_numpy()

# convert temperatures to K
T_K = T_C + 273.15

# define the function to fit to the data
def calc_k(T,log_k0,E):
    return 10**log_k0*math.e**(-E/R/T)

# guess the values of the parameters
guess = np.array([0.,20.])

# fit the function to the data
beta, beta_ci, r_squared = reb.fitNLSR(guess,T_K,k,calc_k,True)

# extract the results
k0 = 10**beta[0]
k0_lower = 10**beta_ci[0,0]
k0_upper = 10**beta_ci[0,1]
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
result.to_csv("./reb_4_4/reb_4_4_Python_results.csv", index=False)

# plot the results
T = np.linspace(min(T_K),max(T_K),100)
k_pred = calc_k(T,beta[0],beta[1])
plt.semilogy(T_K,k,color='r',marker='o', ls='none', label='Measured')
plt.semilogy(T,k_pred,color='k', label='Predicted')
plt.xlabel('T (K)')
plt.ylabel('Rate Coefficient (L mol$^{-1}$ min$^{-1}$)')
plt.legend()
# save and show the figure
plt.savefig('reb_4_4/reb_4_4_Python_fig_1.png')
plt.show()
