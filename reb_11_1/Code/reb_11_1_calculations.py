import pandas as pd
import numpy as np
import math
import rebutils as reb
import matplotlib.pyplot as plt

# Given and known constants
P = 1.0 # atm
D = 1.0 # cm
L = 10.0 # cm
T = 1500. # K
mol_per_cc_at_stp = 1.0E-3/22.4

# path for saving files
filepath1 = './reb_11_1/Results/'
filepath2 = '../RE_Basics/Graphics/'

# Read the experimental data into a dataframe
df = pd.read_csv("reb_11_1/Data/reb_11_1_data.csv")
        # columns: A_in, Y_in, Z_in, fA

# Extract the data as arrays
Vdot_Ain = df['A_in'].to_numpy()
Vdot_Yin = df['Y_in'].to_numpy()
Vdot_Zin = df['Z_in'].to_numpy()
fA = df['fA'].to_numpy()

# combine the adjusted inputs into a matrix
adjusted_inputs = np.transpose(np.array([Vdot_Ain, Vdot_Yin, Vdot_Zin]))

# Allocate storage for the responses
resp = np.ones(len(Vdot_Ain)) * -1

# Define the response function
def response_function(inputs, par):
    # extract the rate coefficient
    k_test = 10.0**par

    # loop through the data points
    for i, input in enumerate(inputs):
        # define the reactor model
        def reactor_model(z,nDot):
            nDot_A = nDot[0]
            nDot_Y = nDot[1]
            nDot_Z = nDot[2]
            PA = nDot_A/(nDot_A + nDot_Y + nDot_Z)*P
            r = k_test*PA
            ddz = np.array([-r, r, r]) * math.pi *D**2/4
            return ddz

        # Solve the reactor model equations
        z0 = 0.0
        n0 = np.array([input[0], input[1], input[2]])*mol_per_cc_at_stp
        f_var = 0
        f_val = L
        soln = reb.solveIVODEs(z0, n0, f_var, f_val, reactor_model)
        if not(soln.success):
            print("A solution was NOT obtained:")
            print(soln.message)
        # Calculate the response
        nA = soln.y[0,-1]
        nAin = input[0]*mol_per_cc_at_stp
        resp[i] = 100*(nAin - nA)/nAin
    return resp

# Estimate the kinetics parameters
par_guess = [-2.0]
beta, beta_ci, r_squared = reb.fitNLSR(par_guess, adjusted_inputs, fA, 
        response_function, False)

# extract the results
k = 10**beta[0]
k_ll = 10**beta_ci[0,0]
k_ul = 10**beta_ci[0,1]

# report the results
print(' ')
print(f'k: {k:.3g} mol/cc/min/atm, 95% CI ['\
      + f'{k_ll:.3g}, {k_ul:.3g}]')
print(f'R-squared: {r_squared:.3g}')
print(' ')

# save the results to a .csv file
data = [['k', f'{k:.3g}', 'mol cm^-3^ min^-1^ atm^-1^'],
    ['k_ll', f'{k_ll:.3g}', 'mol cm^-3^ min^-1^ atm^-1^'],
    ['k_ul', f'{k_ul:.3g}', 'mol cm^-3^ min^-1^ atm^-1^'],
    ['R_squared', f'{r_squared:.3g}', '']]
result = pd.DataFrame(data, columns=['item','value','units'])
result.to_csv(filepath1 + "reb_11_1_results.csv", index=False)

# calculate the model-predicted responses
y_model = response_function(adjusted_inputs,beta[0])

# calculate the residuals
residuals = fA - y_model

# create a parity plot
plt.figure(1) 
plt.plot(fA, y_model, color = 'r', marker='o', ls='')
plt.plot([np.min(fA), np.max(fA)], [np.min(fA), np.max(fA)], color = 'k')
plt.xlabel("experimental response (%)")
plt.ylabel("model-predicted response (%)")

# save and show the parity plot
filename = filepath1 + 'reb_11_1_parity.png'
plt.savefig(filename)
filename = filepath2 + 'reb_11_1_parity.png'
plt.savefig(filename)
plt.show()

# create a residuals plot for the inlet A flow rate
plt.figure(2) 
plt.plot(Vdot_Ain, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Reagent A Inlet Volumetric Flow Rate (sccm)")
plt.ylabel("Residual (%)")

# save and show the residuals plot for the reaction time
filename = filepath1 + 'reb_11_1_residuals_vs_VFR_A.png'
plt.savefig(filename)
filename = filepath2 + 'reb_11_1_residuals_vs_VFR_A.png'
plt.savefig(filename)
plt.show()

# create a residuals plot for the inlet Y flow rate
plt.figure(2) 
plt.plot(Vdot_Yin, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Reagent Y Inlet Volumetric Flow Rate (sccm)")
plt.ylabel("Residual (%)")

# save and show the residuals plot for the reaction time
filename = filepath1 + 'reb_11_1_residuals_vs_VFR_Y.png'
plt.savefig(filename)
filename = filepath2 + 'reb_11_1_residuals_vs_VFR_Y.png'
plt.savefig(filename)
plt.show()

# create a residuals plot for the inlet Z flow rate
plt.figure(2) 
plt.plot(Vdot_Zin, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Reagent Z Inlet Volumetric Flow Rate (sccm)")
plt.ylabel("Residual (%)")

# save and show the residuals plot for the reaction time
filename = filepath1 + 'reb_11_1_residuals_vs_VFR_Z.png'
plt.savefig(filename)
filename = filepath2 + 'reb_11_1_residuals_vs_VFR_Z.png'
plt.savefig(filename)
plt.show()