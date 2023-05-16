import pandas as pd
import numpy as np
import scipy as sp
import rebutils as reb
import matplotlib.pyplot as plt

# Known constant quantities
R = 8.314E-3 # kJ/mol/K
Rpv = 82.06 # cc atm/mol/K
V = 100 # cc
PA0 = 4.0 # atm

# path for saving files
filepath = './reb_9_3/Results/'

# Read the experimental data into a dataframe
df = pd.read_csv("reb_9_3/Data/reb_9_3_data.csv")
        # columns T, t, P

# create a dataframe for the fitting results
phase_1_results_df = pd.DataFrame(columns=['T', 'k', 'k_ll', 'k_ul', 'R_sq'])

# get the temperatures used in the experiments
block_temperatures = df['T'].unique() + 273.15

# process each block separately
for T_K in block_temperatures:
    # create the block
    block = df[df['T'] == T_K - 273.15]

    # extract the data as arrays
    t = block["t"].to_numpy()
    P = block["P"].to_numpy()

    # allocate storage for the responses
    resp = np.zeros(len(t))

    # Define the response function
    def response_function(t, par1):
        # extract the rate coefficient
        k = 10.0**par1

        # loop through the data points
        for i in range(0,len(t)):
            # define the reactor model
            def reactor_model(t,y):
                PA = y[0]*Rpv*T_K/V
                r = k*PA
                ddt = np.array([-2*r*V, r*V])
                return ddt

            # Solve the reactor model equations
            t0 = 0.0
            n0 = np.array([PA0*V/Rpv/T_K, 0.0])
            f_var = 0
            f_val = t[i]
            soln = reb.solveIVODEs(t0, n0, f_var, f_val, reactor_model)
            if not(soln.success):
                print("A solution was NOT obtained:")
                print(soln.message)
            # Calculate the response
            nAf = soln.y[0,-1]
            nA2f = soln.y[1,-1]
            resp[i] = (nAf + nA2f)*Rpv*T_K/V
        return resp

    # Estimate the kinetics parameters
    if len(phase_1_results_df.index) > 0:
        par_guess = beta
    else:
        par_guess = -7.0
    beta, beta_ci, r_squared = reb.fitNLSR(par_guess, t, P, response_function, 
            False)

    # add the fitting results to the dataframe
    phase_1_results_df.loc[len(phase_1_results_df)] = [T_K - 273.15, 
            10.0**beta[0], 10.0**beta_ci[0,0], 10.0**beta_ci[0,1], r_squared]
    
    # calculate the model-predicted y
    y_model = response_function(t,beta[0])

    # calculate the residuals
    residuals = P - y_model

    # create, show, and save a parity plot and a residuals plot
    plt.figure(1) 
    plt.plot(P, y_model, color = 'r', marker='o', ls='')
    plt.plot([np.min(P), np.max(P)], [np.min(P), np.max(P)], color = 'k')
    plt.xlabel("experimental response (atm)")
    plt.ylabel("model-predicted response (atm)")

    # save and show the parity plot
    T_as_text = format(T_K-273.15,'.0f')
    filename = filepath + 'reb_9_3_parity_T_equals_' + T_as_text + '.png'
    plt.savefig(filename)
    plt.show()

    plt.figure(2) 
    plt.plot(t, residuals, color = 'r', marker='o', ls='')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("Reaction time (s)")
    plt.ylabel("Residual (atm)")

    # save and show the parity plot
    filename = filepath + 'reb_9_3_residuals_T_equals_' + T_as_text + '.png'
    plt.savefig(filename)
    plt.show()

# Save and show the phase 1 results
print(phase_1_results_df)
phase_1_results_df.to_csv(filepath + 'reb_9_3_phase_1_results.csv',index=False)

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
print(f'k0: {k0:.3g} mol/cc/s/atm, 95% CI [{k0_lower:.3g}, {k0_upper:.3g}]')
print(f'E: {E:.3g} kcal/mol, 95% CI [{E_lower:.3g}, {E_upper:.3g}]')
print(f'R-squared: {r_squared:.3g}')
print(' ')

# save the results to a .csv file
data = [['k0', f'{k0:.3g}', 'mol cm^-3^ s^-1^ atm^-1^'],
    ['k0_lower_limit', f'{k0_lower:.3g}', 'mol cm^-3^ s^-1^ atm^-1'],
    ['k0_upper_limit', f'{k0_upper:.3g}', 'mol cm^-3^ s^-1^ atm^-1'],
    ['E', f'{E:.3g}', 'kJ mol^-1^'],
    ['E_lower_limit', f'{E_lower:.3g}', 'kcal mol^-1^'],
    ['E_upper_limit', f'{E_upper:.3g}', 'kcal mol^-1^'],
    ['R_squared', f'{r_squared:.3g}', '']]
result = pd.DataFrame(data, columns=['item','value','units'])
result.to_csv(filepath + "reb_9_3_Arrhenius_results.csv", index=False)

# create, show, and save an Arrhenius plot
y_pred = k0*np.exp(-E/R/T)
plt.figure()
plt.semilogy(1/T,k,color='r',marker='o', ls='none')
plt.semilogy(1/T,y_pred,color='k')
plt.xlabel('T$^{-1}$ (K$^{-1}$)')
plt.ylabel('k (min$^{-1}$)')
plt.xticks(rotation=25)
plt.tight_layout()
plt.savefig(filepath + 'reb_9_3_Arrhenius_plot.png')
plt.show()
