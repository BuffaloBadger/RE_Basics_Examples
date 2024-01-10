import pandas as pd
import numpy as np
import response_function as rf
import rebutils as reb
import matplotlib.pyplot as plt

# set filepaths
filepath_to_data = './reb_19_3/Data/'
filepath_to_results = './reb_19_3/Results/'
filepath_to_figures = '../RE_Basics/Graphics/'

# given and known constants
V = 100.0 # cm^3
P0 = 6.0 # atm
T = 275 + 273.15 # K
R = 82.06 # cm^3 atm/mol/K

# Read the experimental data into a dataframe
df = pd.read_csv(filepath_to_data + 'reb_19_3_data.csv')
        # columns: PA0, t, P

# extract the data as arrays
PA0 = df['PA0'].to_numpy()
t = df['t'].to_numpy()
P = df['P'].to_numpy()

# combine the adjusted inputs in a matrix
adjusted_inputs = np.transpose(np.array([PA0, t]))

# make a guess for log_10 of k
par_guess = [-6.0]

# check the guess
check = rf.response(adjusted_inputs[0:1,:],par_guess[0])
print('\nPredicted pressure for first data point using guess: %.3g'%(check[0]))

# estimate log_10 of k
beta, beta_ci, r_squared = reb.fitNLSR(par_guess, adjusted_inputs, P, 
    rf.response, False)

# extract the results
k = 10.**beta[0]
k_ll = 10.**beta_ci[0,0]
k_ul = 10.**beta_ci[0,1]

# Display the results
print(' ')
print(f'k: {k:.3g} mol/cc/min/atm^1.5, 95% CI [{k_ll:.3g}, {k_ul:.3g}]')
print(f'R-squared: {r_squared:.3g}')
print(' ')

# Save the results
data = [['k', f'{k:.3g}', 'mol cm^-3^ min^-1^ atm^-1.5^'],
    ['k_lower_limit', f'{k_ll:.3g}', 'mol cm^-3^ min^-1^ atm^-1.5^'],
    ['k_upper_limit', f'{k_ul:.3g}', 'mol cm^-3^ min^-1^ atm^-1.5^'],
    ['R_squared', f'{r_squared:.3g}', '']]
result = pd.DataFrame(data, columns=['item','value','units'])
result.to_csv(filepath_to_results + '/reb_19_3_results.csv', index=False)

# calculate the model-predicted responses and the residuals
P_model = rf.response(adjusted_inputs, beta[0])
residuals = P - P_model

# create, display and save the parity plot
plt.figure(1) 
plt.plot(P, P_model, color = 'r', marker='o', ls='')
plt.plot([np.min(P), np.max(P)], [np.min(P), np.max(P)], color = 'k')
plt.xlabel("Measured Pressure (atm)")
plt.ylabel("Predicted Pressure (atm)")
plt.savefig(filepath_to_figures + 'reb_19_3_parity.png')
plt.savefig(filepath_to_results + 'reb_19_3_parity.png')
plt.show()

# create, display and save the residuals plots
# residuals vs. PA0
plt.figure(2) 
plt.plot(PA0, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Initial Partial Pressure of A (atm)")
plt.ylabel("Residual (atm)")
plt.savefig(filepath_to_figures + 'reb_19_3_PA0_residuals.png')
plt.savefig(filepath_to_results + 'reb_19_3_PA0_residuals.png')
plt.show()

# residuals vs. tRxn
plt.figure(2) 
plt.plot(t, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Reaction time (min)")
plt.ylabel("Residual (atm)")
plt.savefig(filepath_to_figures + 'reb_19_3_t_residuals.png')
plt.savefig(filepath_to_results + 'reb_19_3_t_residuals.png')
plt.show()