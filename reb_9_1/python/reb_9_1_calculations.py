import pandas as pd
import numpy as np
import scipy as sp
import rebutils as reb
import matplotlib.pyplot as plt

# Given and known constants
V = 0.1 # L
T = 35 + 273.15 # K

# Read the data file
df = pd.read_csv('./reb_9_1/reb_9_1_data.csv')
        # Columns: VFR, CAin, CBin, CYin, CZin, CYout

# Extract the data into arrays and get the number of experiments
VFR = df['VFR'].to_numpy()
CAin = df['CAin'].to_numpy()
CBin = df['CBin'].to_numpy()
CYin = df['CYin'].to_numpy()
CZin = df['CZin'].to_numpy()
CYout = df['CYout'].to_numpy()
n_expts = len(VFR)

# create a matrix containing the adjusted inputs
adj_inputs = np.transpose(np.array([VFR, CAin, CBin, CYin, CZin]))

# allocate storage for the responses
resp = np.zeros(n_expts)

# define the response function
def response_function(inputs, log_k_guess):
    k_guess = 10**log_k_guess
    for i in range(0,n_expts):
        # extract the inputs
        Vdot = inputs[i,0]
        CA0 = inputs[i,1]
        CB0 = inputs[i,2]
        CY0 = inputs[i,3]
        CZ0 = inputs[i,4]
        # calculate inlet molar flow rates
        nAin = CA0*Vdot
        nBin = CB0*Vdot
        nYin = CY0*Vdot
        nZin = CZ0*Vdot

        # define the reactor model residuals
        def reactor_model(x):
            # extract the test values for the unknowns
            nA = x[0]
            nB = x[1]
            nY = x[2]
            nZ = x[3]

            # calculate the rate
            CA = nA/Vdot
            CB = nB/Vdot
            r = k_guess*CA*CB

            # evaluate the reactor model residuals
            r1 = nAin - nA - r*V
            r2 = nBin - nB - r*V
            r3 = nYin - nY + r*V
            r4 = nZin - nZ + r*V
            return [r1, r2, r3, r4]
        
        # guess the solution to the reactor model
        nY_guess = CYout[i]*Vdot
        extent = nY_guess - nYin
        nA_guess = nAin - extent
        nB_guess = nBin - extent
        nZ_guess = nZin + extent
        guess = [nA_guess, nB_guess, nY_guess, nZ_guess]

        # solve the reactor model equations
        soln = sp.optimize.root(reactor_model,guess)
        if not(soln.success):
            print("A solution was NOT obtained:")
            print(soln.message)
        
        # calculate the response, CYout
        resp[i] = soln.x[2]/Vdot
    return resp

# Estimate the parameters (here the base 10 log of k)
par_guess = [0.0]
beta, beta_ci, r_squared = reb.fitNLSR(par_guess, adj_inputs, CYout, 
        response_function, False)

# Extract the results
k = 10.0**beta[0]
k_ll = 10.0**beta_ci[0,0]
k_ul = 10.0**beta_ci[0,1]

# Display the results
print(' ')
print(f'k: {k:.3g} L/mol/min, 95% CI [{k_ll:.3g}, {k_ul:.3g}]')
print(f'R-squared: {r_squared:.3g}')
print(' ')

# Save the results to a .csv file
data = [['k', f'{k:.3g}', 'L mol^-1^ min^-1^'],
    ['k_lower_limit', f'{k_ll:.3g}', 'L mol^-1^ min^-1^'],
    ['k_upper_limit', f'{k_ul:.3g}', 'L mol^-1^ min^-1^'],
    ['R_squared', f'{r_squared:.3g}', '']]
result = pd.DataFrame(data, columns=['item','value','units'])
result.to_csv('./reb_9_1/python/reb_9_1_results.csv', index=False)

# calculate the model-predicted responses and the residuals
y_model = response_function(adj_inputs, beta[0])
residuals = CYout - y_model

# create, display and save the parity plot
plt.figure(1) 
plt.plot(CYout, y_model, color = 'r', marker='o', ls='')
plt.plot([np.min(CYout), np.max(CYout)], [np.min(CYout), np.max(CYout)], 
        color = 'k')
plt.xlabel("experimental response (M)")
plt.ylabel("model-predicted response (M)")
plt.savefig('./reb_9_1/python/reb_9_1_parity.png')
plt.show()

# create, display and save the residuals plots
# residuals vs. VFR
plt.figure(2) 
plt.plot(VFR, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Volumetric Flow Rate (cm$^3$ min$^{-1}$)")
plt.ylabel("Residual (M)")
plt.savefig('./reb_9_1/python/reb_9_1_VFR_residuals.png')
plt.show()

# residuals vs. CAin
plt.figure(2) 
plt.plot(CAin, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Inlet Concentration of A (M)")
plt.ylabel("Residual (M)")
plt.savefig('./reb_9_1/python/reb_9_1_CA_residuals.png')
plt.show()

# residuals vs. CBin
plt.figure(3) 
plt.plot(CBin, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Inlet Concentration of B (M)")
plt.ylabel("Residual (M)")
plt.savefig('./reb_9_1/python/reb_9_1_CB_residuals.png')
plt.show()

# residuals vs. CYin
plt.figure(4) 
plt.plot(CYin, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Inlet Concentration of Y (M)")
plt.ylabel("Residual (M)")
plt.savefig('./reb_9_1/python/reb_9_1_CY_residuals.png')
plt.show()

# residuals vs. CZin
plt.figure(5) 
plt.plot(CZin, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Inlet Concentration of Z (M)")
plt.ylabel("Residual (M)")
plt.savefig('./reb_9_1/python/reb_9_1_CZ_residuals.png')
plt.show()
