import pandas as pd
import numpy as np
import rebutils as reb
import matplotlib.pyplot as plt
import response_function as rf

# Given and known constants
P = 1.0 # atm
m = 3.0 # g
VFR_in = 0.85 # L/min
T = 400. + 273.15 # K
K = 12.2

# path for saving files
filepath1 = './reb_11_2/Results/'
filepath2 = '../RE_Basics/Graphics/'

# Read the experimental data into a dataframe
df = pd.read_csv("reb_11_2/Data/reb_11_2_data.csv")
        # columns: yA, yB, yY, yZ, PA

# Extract the data as arrays
yAin = df['yA'].to_numpy()
yBin = df['yB'].to_numpy()
yYin = df['yY'].to_numpy()
yZin = df['yZ'].to_numpy()
PA = df['PA'].to_numpy()

# combine the adjusted inputs into a matrix
adjusted_inputs = np.transpose(np.array([yAin, yBin, yYin, yZin]))


# Estimate the kinetics parameters
#par_guess = [0.0, 1.0, 1.0, -0.5, -0.5]
par_guess = [-3.0, 1.0, 0.3, -0.7, 0.1]
beta, beta_ci, r_squared = reb.fitNLSR(par_guess, adjusted_inputs, PA, 
        rf.response_function, False)

# extract the results
k = 10**beta[0]
k_ll = 10**beta_ci[0,0]
k_ul = 10**beta_ci[0,1]
alphaA = beta[1]
alphaA_ll = beta_ci[1,0]
alphaA_ul = beta_ci[1,1]
alphaB = beta[2]
alphaB_ll = beta_ci[2,0]
alphaB_ul = beta_ci[2,1]
alphaY = beta[3]
alphaY_ll = beta_ci[3,0]
alphaY_ul = beta_ci[3,1]
alphaZ = beta[4]
alphaZ_ll = beta_ci[4,0]
alphaZ_ul = beta_ci[4,1]

# report the results
print(' ')
print(f'k: {k:.3g} mol/L/min, 95% CI [{k_ll:.3g}, {k_ul:.3g}]')
print(f'alphaA: {alphaA:.3g}, 95% CI [{alphaA_ll:.3g}, {alphaA_ul:.3g}]')
print(f'alphaB: {alphaB:.3g}, 95% CI [{alphaB_ll:.3g}, {alphaB_ul:.3g}]')
print(f'alphaY: {alphaY:.3g}, 95% CI [{alphaY_ll:.3g}, {alphaY_ul:.3g}]')
print(f'alphaZ: {alphaZ:.3g}, 95% CI [{alphaZ_ll:.3g}, {alphaZ_ul:.3g}]')
print(f'R-squared: {r_squared:.3g}')
print(' ')

# save the results to a .csv file
data = [['k', f'{k:.3g}', 'mol L^-1^ min^-1^'],
    ['k_ll', f'{k_ll:.3g}', 'mol L^-1^ min^-1^'],
    ['k_ul', f'{k_ul:.3g}', 'mol L^-1^ min^-1^'],
    ['alpha_A', f'{alphaA:.3g}',''],
    ['alpha_A_ll', f'{alphaA_ll:.3g}',''],
    ['alpha_A_ul', f'{alphaA_ul:.3g}',''],
    ['alpha_B', f'{alphaB:.3g}',''],
    ['alpha_B_ll', f'{alphaB_ll:.3g}',''],
    ['alpha_B_ul', f'{alphaB_ul:.3g}',''],
    ['alpha_Y', f'{alphaY:.3g}',''],
    ['alpha_Y_ll', f'{alphaY_ll:.3g}',''],
    ['alpha_Y_ul', f'{alphaY_ul:.3g}',''],
    ['alpha_Z', f'{alphaZ:.3g}',''],
    ['alpha_Z_ll', f'{alphaZ_ll:.3g}',''],
    ['alpha_Z_ul', f'{alphaZ_ul:.3g}',''],
    ['R_squared', f'{r_squared:.3g}', '']]
result = pd.DataFrame(data, columns=['item','value','units'])
result.to_csv(filepath1 + "reb_11_2_results.csv", index=False)

# calculate the model-predicted responses
y_model = rf.response_function(adjusted_inputs,beta[0],beta[1],beta[2],beta[3],
        beta[4])

# calculate the residuals
residuals = PA - y_model

# create a parity plot
plt.figure(1) 
plt.plot(PA, y_model, color = 'r', marker='o', ls='')
plt.plot([np.min(PA), np.max(PA)], [np.min(PA), np.max(PA)], color = 'k')
plt.xlabel("experimental response (atm)")
plt.ylabel("model-predicted response (atm)")

# save and show the parity plot
filename = filepath1 + 'reb_11_2_parity.png'
plt.savefig(filename)
filename = filepath2 + 'reb_11_2_parity.png'
plt.savefig(filename)
plt.show()

# create a residuals plot for the inlet yA
plt.figure(2) 
plt.plot(yAin, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Inlet A Mole Fraction")
plt.ylabel("Residual (%)")

# save and show the residuals plot for the reaction time
filename = filepath1 + 'reb_11_2_residuals_vs_yAin.png'
plt.savefig(filename)
filename = filepath2 + 'reb_11_2_residuals_vs_yAin.png'
plt.savefig(filename)
plt.show()

# create a residuals plot for the inlet yB
plt.figure(2) 
plt.plot(yBin, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Inlet B Mole Fraction")
plt.ylabel("Residual (%)")

# save and show the residuals plot for the reaction time
filename = filepath1 + 'reb_11_2_residuals_vs_yBin.png'
plt.savefig(filename)
filename = filepath2 + 'reb_11_2_residuals_vs_yBin.png'
plt.savefig(filename)
plt.show()

# create a residuals plot for the inlet yY
plt.figure(2) 
plt.plot(yYin, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Inlet Y Mole Fraction")
plt.ylabel("Residual (%)")

# save and show the residuals plot for the reaction time
filename = filepath1 + 'reb_11_2_residuals_vs_yYin.png'
plt.savefig(filename)
filename = filepath2 + 'reb_11_2_residuals_vs_yYin.png'
plt.savefig(filename)
plt.show()

# create a residuals plot for the inlet yZ
plt.figure(2) 
plt.plot(yZin, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Inlet Z Mole Fraction")
plt.ylabel("Residual (%)")

# save and show the residuals plot for the reaction time
filename = filepath1 + 'reb_11_2_residuals_vs_yZin.png'
plt.savefig(filename)
filename = filepath2 + 'reb_11_2_residuals_vs_yZin.png'
plt.savefig(filename)
plt.show()