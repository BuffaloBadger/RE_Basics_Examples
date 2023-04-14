import pandas as pd
import numpy as np
import rebutils as reb
import matplotlib.pyplot as plt

# Given and known constants
P = 1.0 # atm
m = 3.0 # g
VFR_in = 0.85 # L/min
T = 400. + 273.15 # K
K = 12.2

# Parameter guess
#par_guess = [-2.0, -3.0, -3.0, -3.0, -3.0]
par_guess = [ 10.0**5.09352415, 10.0**-4.30658933,  10.0**7.56602674,  \
             10.0**7.3611074,  10.0**-2.85303979]

# path for saving files
filepath = './reb_10_3/python/'

# Read the experimental data into a dataframe
df = pd.read_csv("reb_10_2/reb_10_2_data.csv")
        # columns: yA, yB, yY, yZ, PA

# Extract the data as arrays
yAin = df['yA'].to_numpy()
yBin = df['yB'].to_numpy()
yYin = df['yY'].to_numpy()
yZin = df['yZ'].to_numpy()
PAexpt = df['PA'].to_numpy()

# combine the adjusted inputs into a matrix
adjusted_inputs = np.transpose(np.array([yAin, yBin, yYin, yZin]))

# Allocate storage for the responses
resp = np.ones(len(yAin)) * -1

# Define the response function
def response_function(inputs, par1, par2, par3, par4, par5):
    # extract the rate coefficient
    k_test = par1
    KA = par2
    KB = par3
    KY = par4
    KZ = par5

    # loop through the data points
    for i, input in enumerate(inputs):
        # define the reactor model
        def reactor_model(z,n):
            nA = n[0]
            nB = n[1]
            nY = n[2]
            nZ = n[3]
            ntot = nA + nB + nY + nZ
            PA = nA/ntot*P
            PB = nB/ntot*P
            PY = nY/ntot*P
            PZ = nZ/ntot*P
            r = k_test*PA*PB/(1+KA*PA+KB*PB+KY*PY+KZ*PZ)*(1 - PY*PZ/K/PA/PB)
            ddm = np.array([-r, -r, r, r])
            return ddm

        # Solve the reactor model equations
        m0 = 0.0
        n0 = np.array([input[0], input[1], input[2], input[3]])*VFR_in/22.4
        f_var = 0
        f_val = m
        soln = reb.solveIVODEs(m0, n0, f_var, f_val, reactor_model)
        if not(soln.success):
            print("A solution was NOT obtained:")
            print(soln.message)
        # Calculate the response
        n_A = soln.y[0,-1]
        n_B = soln.y[1,-1]
        n_Y = soln.y[2,-1]
        n_Z = soln.y[3,-1]
        n_tot = n_A + n_B + n_Y + n_Z
        resp[i] = n_A/n_tot*P
    return resp

# Estimate the kinetics parameters

beta, beta_ci, r_squared = reb.fitNLSR(par_guess, adjusted_inputs, PAexpt, 
        response_function, False)

# extract the results
k = beta[0]
k_ll = beta_ci[0,0]
k_ul = beta_ci[0,1]
K_A = beta[1]
K_A_ll = beta_ci[1,0]
K_A_ul = beta_ci[1,1]
K_B = beta[2]
K_B_ll = beta_ci[2,0]
K_B_ul = beta_ci[2,1]
K_Y = beta[3]
K_Y_ll = beta_ci[3,0]
K_Y_ul = beta_ci[3,1]
K_Z = beta[4]
K_Z_ll = beta_ci[4,0]
K_Z_ul = beta_ci[4,1]

# report the results
print(' ')
print(f'k: {k:.3g} mol/L/min, 95% CI [{k_ll:.3g}, {k_ul:.3g}]')
print(f'K_A: {K_A:.3g}, 95% CI [{K_A_ll:.3g}, {K_A_ul:.3g}]')
print(f'K_B: {K_B:.3g}, 95% CI [{K_B_ll:.3g}, {K_B_ul:.3g}]')
print(f'K_Y: {K_Y:.3g}, 95% CI [{K_Y_ll:.3g}, {K_Y_ul:.3g}]')
print(f'K_Z: {K_Z:.3g}, 95% CI [{K_Z_ll:.3g}, {K_Z_ul:.3g}]')
print(f'R-squared: {r_squared:.3g}')
print(' ')

# save the results to a .csv file
data = [['k', f'{k:.3g}', 'mol L^-1^ min^-1^'],
    ['k_ll', f'{k_ll:.3g}', 'mol L^-1^ min^-1^'],
    ['k_ul', f'{k_ul:.3g}', 'mol L^-1^ min^-1^'],
    ['K_A', f'{K_A:.3g}',''],
    ['K_A_ll', f'{K_A_ll:.3g}',''],
    ['K_A_ul', f'{K_A_ul:.3g}',''],
    ['K_B', f'{K_B:.3g}',''],
    ['K_B_ll', f'{K_B_ll:.3g}',''],
    ['K_B_ul', f'{K_B_ul:.3g}',''],
    ['K_Y', f'{K_Y:.3g}',''],
    ['K_Y_ll', f'{K_Y_ll:.3g}',''],
    ['K_Y_ul', f'{K_Y_ul:.3g}',''],
    ['K_Z', f'{K_Z:.3g}',''],
    ['K_Z_ll', f'{K_Z_ll:.3g}',''],
    ['K_Z_ul', f'{K_Z_ul:.3g}',''],
    ['R_squared', f'{r_squared:.3g}', '']]
result = pd.DataFrame(data, columns=['item','value','units'])
result.to_csv(filepath + "reb_10_3_results.csv", index=False)

# calculate the model-predicted responses
y_model = response_function(adjusted_inputs,beta[0],beta[1],beta[2],beta[3],
        beta[4])

# calculate the residuals
residuals = PAexpt - y_model

# create a parity plot
plt.figure(1) 
plt.plot(PAexpt, y_model, color = 'r', marker='o', ls='')
plt.plot([np.min(PAexpt), np.max(PAexpt)], [np.min(PAexpt), np.max(PAexpt)], color = 'k')
plt.xlabel("experimental response (atm)")
plt.ylabel("model-predicted response (atm)")

# save and show the parity plot
filename = filepath + 'reb_10_3_parity.png'
plt.savefig(filename)
plt.show()

# create a residuals plot for the inlet yA
plt.figure(2) 
plt.plot(yAin, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Inlet A Mole Fraction")
plt.ylabel("Residual (%)")

# save and show the residuals plot for the reaction time
filename = filepath + 'reb_10_3_residuals_vs_yAin.png'
plt.savefig(filename)
plt.show()

# create a residuals plot for the inlet yB
plt.figure(3) 
plt.plot(yBin, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Inlet B Mole Fraction")
plt.ylabel("Residual (%)")

# save and show the residuals plot for the reaction time
filename = filepath + 'reb_10_3_residuals_vs_yBin.png'
plt.savefig(filename)
plt.show()

# create a residuals plot for the inlet yY
plt.figure(4) 
plt.plot(yYin, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Inlet Y Mole Fraction")
plt.ylabel("Residual (%)")

# save and show the residuals plot for the reaction time
filename = filepath + 'reb_10_3_residuals_vs_yYin.png'
plt.savefig(filename)
plt.show()

# create a residuals plot for the inlet yZ
plt.figure(5) 
plt.plot(yZin, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Inlet Z Mole Fraction")
plt.ylabel("Residual (%)")

# save and show the residuals plot for the reaction time
filename = filepath + 'reb_10_3_residuals_vs_yZin.png'
plt.savefig(filename)
plt.show()