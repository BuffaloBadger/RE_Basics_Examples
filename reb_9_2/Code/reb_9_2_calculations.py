import pandas as pd
import numpy as np
import rebutils as reb
import matplotlib.pyplot as plt

# path for saving files
filepath = './reb_9_2/Results/'

# Constant inputs
V = 500.0 # cc
R = 1.987E-3 # kcal/mol/K
Rpv = 82.06 # cc atm/mol/K

# read the data from the .csv file (not in filepath)
df = pd.read_csv('./reb_9_2/Data/reb_9_2_data.csv')
        # columns: T, PA0, PB0, t, fA

# create a dataframe for the fitting results
phase_1_results_df = pd.DataFrame(columns=['T', 'k', 'k_ll', 'k_ul', 'R_sq'])

# get the temperatures used in the experiments
block_temperatures = df['T'].unique() + 273.15

# process each block separately
for T_K in block_temperatures:
    # create the block
    block = df[df['T'] == T_K - 273.15]

    # extract the data as arrays
    PA0 = block['PA0'].to_numpy()
    PB0 = block['PB0'].to_numpy()
    t = block['t'].to_numpy()
    f_A = block['fA'].to_numpy()

    # allocate storage for x and y
    x = np.zeros(len(block.index))
    y = np.zeros(len(block.index))
    
    # calculate x and y 
    for i in range(0,len(x)):
        nA0 = PA0[i]*V/Rpv/T_K
        nB0 = PB0[i]*V/Rpv/T_K
        nA = nA0*(1-f_A[i])
        nB = nB0 - nA0 + nA
        x[i] = -t[i]*(Rpv*T_K)**2/V
        if nA0 == nB0:
            y[i] = 1/nA0 - 1/nA
        else:
            y[i] = 1/(nA0-nB0)*np.log(nA0*nB/nB0/nA)

    # fit y = mx to the data
    beta, beta_ci, r_squared = reb.fitLinSR(x,y,False)

    # add the fitting results to the dataframe
    phase_1_results_df.loc[len(phase_1_results_df)] = [T_K - 273.15, beta[0],
            beta_ci[0,0], beta_ci[0,1], r_squared]
    
    # calculate the model-predicted y
    y_model = beta[0]*x
    
    # create, show, and save a model plot
    plt.figure() 
    plt.plot(x, y_model, color = 'k', ls='-' , label = 'model')
    plt.plot(x, y, color = 'r', marker='o', ls='', label = 'experiment')
    plt.xlabel("x")
    plt.ylabel("y")
    T_as_text = format(T_K-273.15,'.0f')
    plt.legend(title='T = ' + T_as_text + '°C')
    plt.ticklabel_format(axis='x',style='plain')
    plt.xticks(rotation=25)
    plt.tight_layout()
    filename = filepath + 'reb_9_2_model_T_equals_' + T_as_text + '.png'
    plt.savefig(filename)
    plt.savefig('../RE_Basics/Graphics/reb_9_2_model_T_equals_'+ T_as_text + \
                '.png')
    plt.show()

# show and save the fitting results
phase_1_results_df.to_csv(filepath + 'reb_9_2_phase_1_results.csv',index=False)

# fit the Arrhenius expression to the estimated rate coefficients
T = phase_1_results_df['T'].to_numpy() + 273.15
k = phase_1_results_df['k'].to_numpy()
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
print(f'k0: {k0:.3g} cc/mol/min/atm^2, 95% CI [{k0_lower:.3g}, {k0_upper:.3g}]')
print(f'E: {E:.3g} kcal/mol, 95% CI [{E_lower:.3g}, {E_upper:.3g}]')
print(f'R-squared: {r_squared:.3g}')
print(' ')

# save the results to a .csv file
data = [['k0', f'{k0:.3g}', 'mol cm^-3^ min^-1^ atm^-2^'],
    ['k0_lower_limit', f'{k0_lower:.3g}', 'mol cm^-3^ min^-1^ atm^-2'],
    ['k0_upper_limit', f'{k0_upper:.3g}', 'mol cm^-3^ min^-1^ atm^-2'],
    ['E', f'{E:.3g}', 'kJ mol^-1^'],
    ['E_lower_limit', f'{E_lower:.3g}', 'kcal mol^-1^'],
    ['E_upper_limit', f'{E_upper:.3g}', 'kcal mol^-1^'],
    ['R_squared', f'{r_squared:.3g}', '']]
result = pd.DataFrame(data, columns=['item','value','units'])
result.to_csv(filepath + "reb_9_2_Arrhenius_results.csv", index=False)

# create, show, and save an Arrhenius plot
y_pred = k0*np.exp(-E/R/T)
plt.figure()
plt.semilogy(1/T,k,color='r',marker='o', ls='none')
plt.semilogy(1/T,y_pred,color='k')
plt.xlabel('T$^{-1}$ (K$^{-1}$)')
plt.ylabel('k (min$^{-1}$)')
plt.xticks(rotation=25)
plt.tight_layout()
plt.savefig(filepath + 'reb_9_2_Arrhenius_plot.png')
plt.savefig('../RE_Basics/Graphics/reb_9_2_Arrhenius_plot.png')
plt.show()
