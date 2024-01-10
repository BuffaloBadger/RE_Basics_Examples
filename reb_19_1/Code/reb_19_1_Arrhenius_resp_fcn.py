import pandas as pd
import numpy as np
import rebutils as reb
import matplotlib.pyplot as plt

# filepaths
filepath_to_data = './reb_19_1/Data/'
filepath_to_results = './reb_19_1/Results/'
filepath_to_figures = '../RE_Basics/Graphics/'

# Known and given constants
R = 8.314E-3 # kJ/mol/K

# read the rate coefficient datadata from the .csv file (not in filepath)
df = pd.read_csv(filepath_to_results + 'reb_19_1_resp_fcn_params.csv')
    # columns: 'T', 'k', 'k_ll', 'k_ul', 'R_sq'

# fit the Arrhenius expression to the estimated rate coefficients
T = df['T'].to_numpy() + 273.15
k = df['k'].to_numpy()
beta, beta_ci, r_squared = reb.fit_Arrhenius_model(k,T,R)

# extract the results
k0 = beta[0]
k0_lower = beta_ci[0,0]
k0_upper = beta_ci[0,1]
E = beta[1]
E_lower = beta_ci[1,0]
E_upper = beta_ci[1,1]

# report the results
print(' ')
print(f'k0: {k0:.3g} L/mol/min, 95% CI [{k0_lower:.3g}, {k0_upper:.3g}]')
print(f'E: {E:.3g} kJ/mol, 95% CI [{E_lower:.3g}, {E_upper:.3g}]')
print(f'R-squared: {r_squared:.3g}')
print(' ')

# save the results to a .csv file
data = [['k0', f'{k0:.3g}', 'min^-1^'],
    ['k0_lower_limit', f'{k0_lower:.3g}', 'min^-1^'],
    ['k0_upper_limit', f'{k0_upper:.3g}', 'min^-1^'],
    ['E', f'{E:.3g}', 'kJ mol^-1^'],
    ['E_lower_limit', f'{E_lower:.3g}', 'kJ mol^-1^'],
    ['E_upper_limit', f'{E_upper:.3g}', 'kJ mol^-1^'],
    ['R_squared', f'{r_squared:.3g}', '']]
result = pd.DataFrame(data, columns=['item','value','units'])
result.to_csv(filepath_to_results + "reb_19_1_Arrhenius_resp_fcn.csv", \
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
filename = 'reb_19_1_Arrhenius_resp_fcn.png'
plt.savefig(filepath_to_results + filename)
plt.savefig(filepath_to_figures + filename)
plt.show()
