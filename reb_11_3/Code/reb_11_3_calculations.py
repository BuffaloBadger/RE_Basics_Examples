import pandas as pd
import numpy as np
import response_function as rf
import rebutils as reb
import matplotlib.pyplot as plt

# Given and known constants
V = 3.0 # gal

# Read the data file
df = pd.read_csv('./reb_11_3/Data/reb_11_3_data.csv')
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

# Estimate the parameters (here the base 10 log of k)
par_guess = [0.0, 0.0]
beta, beta_ci, r_squared = reb.fitNLSR(par_guess, adj_inputs, CAout, 
        rf.response_function, False)

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
result.to_csv('./reb_11_3/Results/reb_11_3_results.csv', index=False)

# calculate the model-predicted responses and the residuals
y_model = rf.response_function(adj_inputs, beta[0], beta[1])
residuals = CAout - y_model

# create, display and save the parity plot
plt.figure(1) 
plt.plot(CAout, y_model, color = 'r', marker='o', ls='')
plt.plot([np.min(CAout), np.max(CAout)], [np.min(CAout), np.max(CAout)], 
        color = 'k')
plt.xlabel("experimental response (lbmol gal$^{-1}$)")
plt.ylabel("model-predicted response (lbmol gal$^{-1}$)")
plt.tight_layout()
plt.savefig('./reb_11_3/Results/reb_11_3_parity.png')
plt.savefig('../RE_Basics/Graphics/reb_11_3_parity.png')
plt.show()

# create, display and save the residuals plots
# residuals vs. VFR
plt.figure(2) 
plt.plot(VFR, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Volumetric Flow Rate (gal min$^{-1}$)")
plt.ylabel("Residual (lbmol gal$^{-1}$)")
plt.tight_layout()
plt.savefig('./reb_11_3/Results/reb_11_3_VFR_residuals.png')
plt.savefig('../RE_Basics/Graphics/reb_11_3_VFR_residuals.png')
plt.show()

# residuals vs. CAin
plt.figure(2) 
plt.plot(CAin, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Inlet Concentration of A (lbmol gal$^{-1}$)")
plt.ylabel("Residual (lbmol gal$^{-1}$)")
plt.tight_layout()
plt.savefig('./reb_11_3/Results/reb_11_3_CA_residuals.png')
plt.savefig('../RE_Basics/Graphics/reb_11_3_CA_residuals.png')
plt.show()

# residuals vs. CYin
plt.figure(4) 
plt.plot(CYin, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Inlet Concentration of Y (lbmol gal$^{-1}$)")
plt.ylabel("Residual (lbmol gal$^{-1}$)")
plt.tight_layout()
plt.savefig('./reb_11_3/Results/reb_11_3_CY_residuals.png')
plt.savefig('../RE_Basics/Graphics/reb_11_3_CY_residuals.png')
plt.show()

# residuals vs. CZin
plt.figure(5) 
plt.plot(CZin, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Inlet Concentration of Z (lbmol gal$^{-1}$)")
plt.ylabel("Residual (lbmol gal$^{-1}$)")
plt.tight_layout()
plt.savefig('./reb_11_3/Results/reb_11_3_CZ_residuals.png')
plt.savefig('../RE_Basics/Graphics/reb_11_3_CZ_residuals.png')
plt.show()
