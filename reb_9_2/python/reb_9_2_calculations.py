import pandas as pd
import numpy as np
import scipy as sp
import rebutils as reb
import matplotlib.pyplot as plt

# Given and known constants
P = 3.0 # atm
T = 450 + 273.15 # K
V = 1.0 # L (basis)
R = 0.08206 # L*atm/mol/K

# Read the data file
df = pd.read_csv('./reb_9_2/reb_9_2_data.csv')
        # Columns: tau, yAin, CZout

# Extract the data into arrays and get the number of experiments
tau = df['tau'].to_numpy()
yAin = df['yAin'].to_numpy()
CZout = df['CZout'].to_numpy()/1000.0
n_expts = len(tau)

# create a matrix containing the adjusted inputs
adj_inputs = np.transpose(np.array([tau, yAin]))

# allocate storage for the responses
resp = np.zeros(n_expts)

# define the response function
def response_function(inputs, log_k_guess, alphaA_guess, alphaB_guess):
    k_guess = 10**log_k_guess
    for i in range(0,n_expts):
        # extract the inputs
        space_time = inputs[i,0]
        yA0 = inputs[i,1]
        # calculate inlet molar flow rates
        VFRin = V/space_time
        nAin = yA0*P*VFRin/R/T
        nBin = (1-yA0)*P*VFRin/R/T

        # define the reactor model residuals
        def reactor_model(x):
            # extract the test values for the unknowns
            nA = x[0]
            nB = x[1]
            nZ = x[2]

            # calculate the rate
            PA = nA*P/(nA+nB+nZ)
            PB = nB*P/(nA+nB+nZ)
            r = k_guess*PA**alphaA_guess*PB**alphaB_guess

            # evaluate the reactor model residuals
            r1 = nAin - nA - r*V
            r2 = nBin - nB - r*V
            r3 = -nZ + r*V
            return [r1, r2, r3]
        
        # guess the solution to the reactor model
        nZ_guess = CZout[i]*VFRin
        nA_guess = nAin - nZ_guess
        nB_guess = nBin - nZ_guess
        guess = [nA_guess, nB_guess, nZ_guess]

        # solve the reactor model equations
        soln = sp.optimize.root(reactor_model,guess)
        if not(soln.success):
            print("A solution was NOT obtained:")
            print(soln.message)
        
        # calculate the response, CYout
        nAf = soln.x[0]
        nBf = soln.x[1]
        nZf = soln.x[2]
        resp[i] = nZf*P/(nAf + nBf + nZf)/R/T
    return resp

# Estimate the parameters (here the base 10 log of k)
#par_guess = [np.log10(4.4e-5), 1.4, 0.6]
par_guess = [-5.0, 1.3, 0.5]
beta, beta_ci, r_squared = reb.fitNLSR(par_guess, adj_inputs, CZout, 
        response_function, False)

# Extract the results
k = 10.0**beta[0]
k_ll = 10.0**beta_ci[0,0]
k_ul = 10.0**beta_ci[0,1]
alphaA = beta[1]
alphaA_ll = beta_ci[1,0]
alphaA_ul = beta_ci[1,1]
alphaB = beta[2]
alphaB_ll = beta_ci[2,0]
alphaB_ul = beta_ci[2,1]

# Display the results
print(' ')
print(f'k: {k:.3g} L/mol/min, 95% CI [{k_ll:.3g}, {k_ul:.3g}]')
print(f'alpha_A: {alphaA:.3g}, 95% CI [{alphaA_ll:.3g}, {alphaA_ul:.3g}]')
print(f'alpha_B: {alphaB:.3g}, 95% CI [{alphaB_ll:.3g}, {alphaB_ul:.3g}]')
print(f'R-squared: {r_squared:.3g}')
print(' ')

# Save the results to a .csv file
order = alphaA + alphaB
order_as_text = format(order,'.3f')
data = [['k', f'{k:.3g}', 'mol L^-1^ s^-1^ atm^' + order_as_text + '^'],
    ['k_lower_limit', f'{k_ll:.3g}', 'mol L^-1^ s^-1^ atm^' + order_as_text + '^'],
    ['k0_upper_limit', f'{k_ul:.3g}', 'mol L^-1^ s^-1^ atm^' + order_as_text + '^'],
    ['alphaA', f'{alphaA:.3g}', ''],
    ['alphaA_ll', f'{alphaA_ll:.3g}', ''],
    ['alphaA_ul', f'{alphaA_ul:.3g}', ''],
    ['alphaB', f'{alphaB:.3g}', ''],
    ['alphaB_ll', f'{alphaB_ll:.3g}', ''],
    ['alphaB_ul', f'{alphaB_ul:.3g}', ''],
    ['R_squared', f'{r_squared:.3g}', '']]
result = pd.DataFrame(data, columns=['item','value','units'])
result.to_csv('./reb_9_2/python/reb_9_2_results.csv', index=False)

# calculate the model-predicted responses and the residuals
y_model = response_function(adj_inputs, beta[0], beta[1], beta[2])
residuals = CZout - y_model

# create, display and save the parity plot
plt.figure(1) 
plt.plot(CZout, y_model, color = 'r', marker='o', ls='')
plt.plot([np.min(CZout), np.max(CZout)], [np.min(CZout), np.max(CZout)], 
        color = 'k')
plt.xlabel("experimental response (M)")
plt.ylabel("model-predicted response (M)")
plt.savefig('./reb_9_2/python/reb_9_2_parity.png')
plt.show()

# create, display and save the residuals plots
# residuals vs. tau
plt.figure(2) 
plt.plot(tau, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Space Time (s)")
plt.ylabel("Residual (M)")
plt.savefig('./reb_9_2/python/reb_9_2_tau_residuals.png')
plt.show()

# residuals vs. yAin
plt.figure(2) 
plt.plot(yAin, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Inlet Mole Fraction of A")
plt.ylabel("Residual (M)")
plt.savefig('./reb_9_2/python/reb_9_2_yA_residuals.png')
plt.show()