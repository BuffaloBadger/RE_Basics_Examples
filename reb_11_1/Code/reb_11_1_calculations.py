import pandas as pd
import numpy as np
import response_function as rf
import rebutils as reb
import matplotlib.pyplot as plt

# Given and known constants
V = 0.1 # L
T = 35 + 273.15 # K

# Read the data file
df = pd.read_csv('./reb_11_1/Data/reb_11_1_data.csv')
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

# Estimate the parameters (here the base 10 log of k)
par_guess = [0.0]
beta, beta_ci, r_squared = reb.fitNLSR(par_guess, adj_inputs, CYout, 
        rf.response_function, False)

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
result.to_csv('./reb_11_1/Results/reb_11_1_results.csv', index=False)

# calculate the model-predicted responses and the residuals
y_model = rf.response_function(adj_inputs, beta[0])
residuals = CYout - y_model

# create, display and save the parity plot
plt.figure(1) 
plt.plot(CYout, y_model, color = 'r', marker='o', ls='')
plt.plot([np.min(CYout), np.max(CYout)], [np.min(CYout), np.max(CYout)], 
        color = 'k')
plt.xlabel("experimental response (M)")
plt.ylabel("model-predicted response (M)")
plt.tight_layout()
plt.savefig('./reb_11_1/Results/reb_11_1_parity.png')
plt.savefig('../RE_Basics/Graphics/reb_11_1_parity.png')
plt.show()

# create, display and save the residuals plots
# residuals vs. VFR
plt.figure(2) 
plt.plot(VFR, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Volumetric Flow Rate (cm$^3$ min$^{-1}$)")
plt.ylabel("Residual (M)")
plt.tight_layout()
plt.savefig('./reb_11_1/Results/reb_11_1_VFR_residuals.png')
plt.savefig('../RE_Basics/Graphics/reb_11_1_VFR_residuals.png')
plt.show()

# residuals vs. CAin
plt.figure(2) 
plt.plot(CAin, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Inlet Concentration of A (M)")
plt.ylabel("Residual (M)")
plt.tight_layout()
plt.savefig('./reb_11_1/Results/reb_11_1_CA_residuals.png')
plt.savefig('../RE_Basics/Graphics/reb_11_1_CA_residuals.png')
plt.show()

# residuals vs. CBin
plt.figure(3) 
plt.plot(CBin, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Inlet Concentration of B (M)")
plt.ylabel("Residual (M)")
plt.tight_layout()
plt.savefig('./reb_11_1/Results/reb_11_1_CB_residuals.png')
plt.savefig('../RE_Basics/Graphics/reb_11_1_CB_residuals.png')
plt.show()

# residuals vs. CYin
plt.figure(4) 
plt.plot(CYin, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Inlet Concentration of Y (M)")
plt.ylabel("Residual (M)")
plt.tight_layout()
plt.savefig('./reb_11_1/Results/reb_11_1_CY_residuals.png')
plt.savefig('../RE_Basics/Graphics/reb_11_1_CY_residuals.png')
plt.show()

# residuals vs. CZin
plt.figure(5) 
plt.plot(CZin, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Inlet Concentration of Z (M)")
plt.ylabel("Residual (M)")
plt.tight_layout()
plt.savefig('./reb_11_1/Results/reb_11_1_CZ_residuals.png')
plt.savefig('../RE_Basics/Graphics/reb_11_1_CZ_residuals.png')
plt.show()
