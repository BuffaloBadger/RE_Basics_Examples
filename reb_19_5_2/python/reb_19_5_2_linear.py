"""Calculations for Reaction Engineering Basics Example 19.5.2"""

# import libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from score_utils import Arrhenius_parameters

# given or known
V = 500.0 # cc
Re = 1.987E-3 # kcal/mol/K
Rpv = 82.06 # cc atm/mol/K

# read the data
df = pd.read_csv('reb_19_5_2/reb_19_5_2_data.csv')

# get the temperatures of the data blocks
block_temperatures = df['T'].unique()

# create a dataframe for the fitting results
parameter_estimates = pd.DataFrame(columns=['T', 'k', 'k_ll', 'k_ul', 'R_sq'])

# process the data blocks
for T in block_temperatures:
    # create the block
    block = df[df['T'] == T]
    n_data = len(block.index)

    # extract the data as arrays
    T = block['T'].to_numpy()
    T_K = T + 273.15
    PA0 = block['PA0'].to_numpy()
    PB0 = block['PB0'].to_numpy()
    tf = block['tf'].to_numpy()
    fA = block['fA'].to_numpy()

    # allocate storage for x and y
    x = np.zeros(n_data)
    y = np.zeros(n_data)
    
    # calculate x and y 
    for i in range(0,n_data):
        nA0 = PA0[i]*V/Rpv/T_K[i]
        nB0 = PB0[i]*V/Rpv/T_K[i]
        nA = nA0*(1-fA[i])
        nB = nB0 - nA0 + nA
        x[i] = -tf[i]*(Rpv*T_K[i])**2/V
        if nA0 == nB0:
            y[i] = 1/nA0 - 1/nA
        else:
            y[i] = 1/(nA0-nB0)*np.log(nA0*nB/nB0/nA)

    # fit y = mx to the data
    model = sm.OLS(y,x)
    res = model.fit()
    m = res.params[0]
    m_ci = res.conf_int(alpha=0.05)

    # add the fitting results to the dataframe
    parameter_estimates.loc[len(parameter_estimates)] = [T[0]
        , m, m_ci[0,0], m_ci[0,1], res.rsquared]

    # calculate the model-predicted y
    y_model = m*x
    
    # create, show, and save a model plot
    plt.figure() 
    plt.plot(x, y_model, color = 'k', ls='-' , label = 'model')
    plt.plot(x, y, color = 'r', marker='o', ls='', label = 'experiment')
    plt.xlabel("x")
    plt.ylabel("y")
    T_as_text = format(T[0],'.0f')
    plt.legend(title='T = ' + T_as_text + '°C')
    plt.ticklabel_format(axis='x',style='plain')
    plt.xticks(rotation=25)
    plt.tight_layout()
    filename = 'reb_19_5_2_model_' + T_as_text + '.png'
    plt.savefig('reb_19_5_2/python/' + filename)
    plt.show()

# show and save the results
print('\nParameter Estimates:\n')
print(parameter_estimates)
parameter_estimates.to_csv('reb_19_5_2/python/reb_19_5_2_linear_results.csv',\
    index=False)

# fit the Arrhenius expression to the k vs T data
T_block = parameter_estimates['T'].to_numpy() + 273.15
k = parameter_estimates['k'].to_numpy()
# fit the Arrhenius expression to the data
k0, k0_ci, E, E_ci, r_squared = Arrhenius_parameters(k,T_block,Re)

# report the results
print(' ')
print(f'k0: {k0:.3g} /min, 95% CI [{k0_ci[0]:.3g}, {k0_ci[1]:.3g}]')
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
result.to_csv("reb_19_5_2/python/reb_19_5_2_Arrhenius_params.csv", 
            index=False)

# create, show, and save an Arrhenius plot
y_pred = k0*np.exp(-E/Re/T_block)
plt.figure()
plt.semilogy(1/T_block,k,color='r',marker='o', ls='none')
plt.semilogy(1/T_block,y_pred,color='k')
plt.xlabel('T$^{-1}$ (K$^{-1}$)')
plt.ylabel('k (min$^{-1}$)')
plt.xticks(rotation=25)
plt.tight_layout()
plt.savefig('reb_19_5_2/python/reb_19_5_2_Arrhenius_plot.png')
plt.show()
