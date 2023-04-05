import pandas as pd
import numpy as np
import scipy as sp
import rebutils as reb
import matplotlib.pyplot as plt

# Known constant quantities
V = 50.0E-3 # L

# path for saving files
filepath = './reb_8_4/python/'

# Read the experimental data into a dataframe
df = pd.read_csv("reb_8_4/reb_8_4_data.csv")
        # columns CS0, t, CP

# extract the data as arrays
CS0 = df["CS0"].to_numpy()
t = df["t"].to_numpy()
CP = df["CP"].to_numpy()

# combine the adjusted inputs into a matrix
adjusted_inputs = np.transpose(np.array([CS0, t]))

# allocate storage for the responses
resp = np.zeros(len(t))

# Define the response function
def response_function(inputs, log_V_max_guess, log_K_m_guess):
    # extract V_max_guess and K_m_guess
    V_max_guess = 10.0**log_V_max_guess
    K_m_guess = 10.0**log_K_m_guess

    # loop through the data points
    for i, input in enumerate(inputs):
        # define the reactor model
        def reactor_model(t,y):
            CS = y[0]/V
            r = V_max_guess*CS/(K_m_guess + CS)
            ddt = np.array([-r*V, r*V, r*V])
            return ddt

        # Solve the reactor model equations
        t0 = 0.0
        n0 = np.array([input[0]*V, 0.0, 0.0])
        f_var = 0
        f_val = input[1]
        soln = reb.solveIVODEs(t0, n0, f_var, f_val, reactor_model)
        if not(soln.success):
            print("A solution was NOT obtained:")
            print(soln.message)
        # Calculate the response
        resp[i] = soln.y[1,-1]/V
    return resp

# Estimate the kinetics parameters
par_guess = [0.0, 0.0]
beta, beta_ci, r_squared = reb.fitNLSR(par_guess, adjusted_inputs, CP, 
        response_function, False)

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
result.to_csv(filepath + "reb_8_4_results.csv", index=False)

# calculate the model-predicted responses
y_model = response_function(adjusted_inputs,beta[0],beta[1])

# calculate the residuals
residuals = CP - y_model

# create a parity plot
plt.figure(1) 
plt.plot(CP, y_model, color = 'r', marker='o', ls='')
plt.plot([np.min(CP), np.max(CP)], [np.min(CP), np.max(CP)], color = 'k')
plt.xlabel("experimental response (mmol L$^{-1}$)")
plt.ylabel("model-predicted response (mmol L$^{-1}$)")

# save and show the parity plot
filename = filepath + 'reb_8_4_parity.png'
plt.savefig(filename)
plt.show()

# create a residuals plot for the reaction time
plt.figure(2) 
plt.plot(t, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Reaction time (s)")
plt.ylabel("Residual (mmol L$^{-1}$)")

# save and show the residuals plot for the reaction time
filename = filepath + 'reb_8_4_residuals_vs_t.png'
plt.savefig(filename)
plt.show()

# create a residuals plot for the initial substrate concentration
plt.figure(3) 
plt.plot(CS0, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Initial Substrate Concentration (mmol L$^{-1}$)")
plt.ylabel("Residual (mmol L$^{-1}$)")

# save and show the residuals plot for the initial substrate concentration
filename = filepath + 'reb_8_4_residuals_vs_CS.png'
plt.savefig(filename)
plt.show()