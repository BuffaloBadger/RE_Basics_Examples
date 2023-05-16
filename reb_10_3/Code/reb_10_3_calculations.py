import pandas as pd
import numpy as np
import scipy as sp
import rebutils as reb
import matplotlib.pyplot as plt

# Given and known constants
V = 3.0 # gal

# Read the data file
df = pd.read_csv('./reb_10_3/Data/reb_10_3_data.csv')
        # Columns: VFR, CAin, CYin, CZin, CA

# Extract the data into arrays and get the number of experiments
VFR = df['VFR'].to_numpy()
CAin = df['CAin'].to_numpy()
CYin = df['CYin'].to_numpy()
CZin = df['CZin'].to_numpy()
CAout = df['CA'].to_numpy()
n_expts = len(VFR)

# create a matrix containing the adjusted inputs
adj_inputs = np.transpose(np.array([VFR, CAin, CYin, CZin]))

# allocate storage for the responses
resp = np.zeros(n_expts)

# define the response function
def response_function(inputs, log_kf_guess, log_kr_guess):
    kf_guess = 10.0**log_kf_guess
    kr_guess = 10.0**log_kr_guess
    for i in range(0,n_expts):
        # extract the inputs
        Vdot = inputs[i,0]
        CA0 = inputs[i,1]
        CY0 = inputs[i,2]
        CZ0 = inputs[i,3]
        # calculate inlet molar flow rates
        nAin = CA0*Vdot
        nYin = CY0*Vdot
        nZin = CZ0*Vdot

        # define the reactor model residuals
        def reactor_model(x):
            # extract the test values for the unknowns
            nA = x[0]
            nY = x[1]
            nZ = x[2]

            # calculate the rate
            CA = nA/Vdot
            CY = nY/Vdot
            CZ = nZ/Vdot
            r = kf_guess*CA**2 - kr_guess*CY*CZ

            # evaluate the reactor model residuals
            r1 = nAin - nA - r*V
            r2 = nYin - nY + r*V
            r3 = nZin - nZ + r*V
            return [r1, r2, r3]
        
        # guess the solution to the reactor model
        nA_guess = CAout[i]*Vdot
        extent = nAin - nA_guess
        nY_guess = nYin + extent
        nZ_guess = nZin + extent
        guess = [nA_guess, nY_guess, nZ_guess]

        # solve the reactor model equations
        soln = sp.optimize.root(reactor_model,guess)
        if not(soln.success):
            print("A solution was NOT obtained:")
            print(soln.message)
        
        # calculate the response, CYout
        resp[i] = soln.x[0]/Vdot
    return resp

# Estimate the parameters (here the base 10 log of k)
par_guess = [0.0, 0.0]
beta, beta_ci, r_squared = reb.fitNLSR(par_guess, adj_inputs, CAout, 
        response_function, False)

# Extract the results
kf = 10.0**beta[0]
kf_ll = 10.0**beta_ci[0,0]
kf_ul = 10.0**beta_ci[0,1]
kr = 10.0**beta[1]
kr_ll = 10.0**beta_ci[1,0]
kr_ul = 10.0**beta_ci[1,1]

# Display the results
print(' ')
print(f'kf: {kf:.3g} gal/lbmol/min, 95% CI [{kf_ll:.3g}, {kf_ul:.3g}]')
print(f'kr: {kr:.3g} gal/lbmol/min, 95% CI [{kr_ll:.3g}, {kr_ul:.3g}]')
print(f'R-squared: {r_squared:.3g}')
print(' ')

# Save the results to a .csv file
data = [['kf', f'{kf:.3g}', 'gal lbmol^-1^ min^-1^'],
    ['kf_lower_limit', f'{kf_ll:.3g}', 'gal lbmol^-1^ min^-1^'],
    ['kf_upper_limit', f'{kf_ul:.3g}', 'gal lbmol^-1^ min^-1^'],
    ['kr', f'{kr:.3g}', 'gal lbmol^-1^ min^-1^'],
    ['kr_lower_limit', f'{kr_ll:.3g}', 'gal lbmol^-1^ min^-1^'],
    ['kr_upper_limit', f'{kr_ul:.3g}', 'gal lbmol^-1^ min^-1^'],
    ['R_squared', f'{r_squared:.3g}', '']]
result = pd.DataFrame(data, columns=['item','value','units'])
result.to_csv('./reb_10_3/Results/reb_10_3_results.csv', index=False)

# calculate the model-predicted responses and the residuals
y_model = response_function(adj_inputs, beta[0], beta[1])
residuals = CAout - y_model

# create, display and save the parity plot
plt.figure(1) 
plt.plot(CAout, y_model, color = 'r', marker='o', ls='')
plt.plot([np.min(CAout), np.max(CAout)], [np.min(CAout), np.max(CAout)], 
        color = 'k')
plt.xlabel("experimental response (lbmol gal$^{-1}$)")
plt.ylabel("model-predicted response (lbmol gal$^{-1}$)")
plt.savefig('./reb_10_3/Results/reb_10_3_parity.png')
plt.savefig('../RE_Basics/Graphics/reb_10_3_parity.png')
plt.show()

# create, display and save the residuals plots
# residuals vs. VFR
plt.figure(2) 
plt.plot(VFR, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Volumetric Flow Rate (gal min$^{-1}$)")
plt.ylabel("Residual (lbmol gal$^{-1}$)")
plt.savefig('./reb_10_3/Results/reb_10_3_VFR_residuals.png')
plt.savefig('../RE_Basics/Graphics/reb_10_3_VFR_residuals.png')
plt.show()

# residuals vs. CAin
plt.figure(2) 
plt.plot(CAin, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Inlet Concentration of A (lbmol gal$^{-1}$)")
plt.ylabel("Residual (lbmol gal$^{-1}$)")
plt.savefig('./reb_10_3/Results/reb_10_3_CA_residuals.png')
plt.savefig('../RE_Basics/Graphics/reb_10_3_CA_residuals.png')
plt.show()

# residuals vs. CYin
plt.figure(4) 
plt.plot(CYin, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Inlet Concentration of Y (lbmol gal$^{-1}$)")
plt.ylabel("Residual (lbmol gal$^{-1}$)")
plt.savefig('./reb_10_3/Results/reb_10_3_CY_residuals.png')
plt.savefig('../RE_Basics/Graphics/reb_10_3_CY_residuals.png')
plt.show()

# residuals vs. CZin
plt.figure(5) 
plt.plot(CZin, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Inlet Concentration of Z (lbmol gal$^{-1}$)")
plt.ylabel("Residual (lbmol gal$^{-1}$)")
plt.savefig('./reb_10_3/Results/reb_10_3_CZ_residuals.png')
plt.savefig('../RE_Basics/Graphics/reb_10_3_CZ_residuals.png')
plt.show()
