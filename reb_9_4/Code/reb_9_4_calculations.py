import pandas as pd
import numpy as np
import scipy as sp
import response_function as rf
import rebutils as reb
import matplotlib.pyplot as plt

# set filepaths
filepath_to_data = './reb_9_4/Data/'
filepath_to_results = './reb_9_4/Results/'
filepath_to_figures = '../RE_Basics/Graphics/'

# Known constant quantities
V = 50.0E-3 # L

# Read the experimental data into a dataframe
df = pd.read_csv("reb_9_4/Data/reb_9_4_data.csv")
        # columns CS0, t, CP

# extract the data as arrays
CS0 = df["CS0"].to_numpy()
t = df["t"].to_numpy()
CP = df["CP"].to_numpy()

# combine the adjusted inputs into a matrix
adjusted_inputs = np.transpose(np.array([CS0, t]))

# allocate storage for the responses
resp = np.zeros(len(t))

# Estimate the kinetics parameters
par_guess = [0.0, 0.0]
beta, beta_ci, r_squared = reb.fitNLSR(par_guess, adjusted_inputs, CP, 
        rf.response, False)

# extract the results
Vmax = 10**beta[0]
Vmax_lower = 10**beta_ci[0,0]
Vmax_upper = 10**beta_ci[0,1]
Km = 10**beta[1]
Km_lower = 10**beta_ci[1,0]
Km_upper = 10**beta_ci[1,1]

# report the results
print(' ')
print(f'Vmax: {Vmax:.3g} mmol/L/min, 95% CI ['\
      + f'{Vmax_lower:.3g}, {Vmax_upper:.3g}]')
print(f'Km: {Km:.3g} mmol/L, 95% CI [{Km_lower:.3g}, {Km_upper:.3g}]')
print(f'R-squared: {r_squared:.3g}')
print(' ')

# save the results to a .csv file
data = [['Vmax', f'{Vmax:.3g}', 'mmol L^-1^ min^-1^'],
    ['Vmax_lower_limit', f'{Vmax_lower:.3g}', 'mmol L^-1^ min^-1^'],
    ['Vmax_upper_limit', f'{Vmax_upper:.3g}', 'mmol L^-1^ min^-1^'],
    ['Km', f'{Km:.3g}', 'mmol L^-1^'],
    ['Km_lower_limit', f'{Km_lower:.3g}', 'mmol L^-1^'],
    ['Km_upper_limit', f'{Km_upper:.3g}', 'mmol L^-1^'],
    ['R_squared', f'{r_squared:.3g}', '']]
result = pd.DataFrame(data, columns=['item','value','units'])
result.to_csv(filepath_to_results + "reb_9_4_results.csv", index=False)

# calculate the model-predicted responses
CP_model = rf.response(adjusted_inputs,beta[0],beta[1])

# calculate the residuals
residuals = CP - CP_model

# create a parity plot
plt.figure(1) 
plt.plot(CP, CP_model, color = 'r', marker='o', ls='')
plt.plot([np.min(CP), np.max(CP)], [np.min(CP), np.max(CP)], color = 'k')
plt.xlabel("experimental response (mmol L$^{-1}$)")
plt.ylabel("model-predicted response (mmol L$^{-1}$)")

# save and show the parity plot
plt.savefig(filepath_to_results + 'reb_9_4_parity.png')
plt.savefig(filepath_to_figures + 'reb_9_4_parity.png')
plt.show()

# create a residuals plot for the reaction time
plt.figure(2) 
plt.plot(t, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Reaction time (s)")
plt.ylabel("Residual (mmol L$^{-1}$)")

# save and show the residuals plot for the reaction time
plt.savefig(filepath_to_results + 'reb_9_4_residuals_vs_t.png')
plt.savefig(filepath_to_figures + 'reb_9_4_residuals_vs_t.png')
plt.show()

# create a residuals plot for the initial substrate concentration
plt.figure(3) 
plt.plot(CS0, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Initial Substrate Concentration (mmol L$^{-1}$)")
plt.ylabel("Residual (mmol L$^{-1}$)")

# save and show the residuals plot for the initial substrate concentrationfilename = filepath + 'reb_9_4_residuals_vs_CS.png'
plt.savefig(filepath_to_results + 'reb_9_4_residuals_vs_CS0.png')
plt.savefig(filepath_to_figures + 'reb_9_4_residuals_vs_CS0.png')
plt.show()