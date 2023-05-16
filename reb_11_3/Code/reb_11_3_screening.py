import pandas as pd
import numpy as np
import rebutils as reb
import matplotlib.pyplot as plt

par_guess = [-2.0, -1.0, 0.0, 0.0, -1.0]
#par_guess = [5.58, -9.75, 8.05, 7.84, -1.32]
#par_guess = [5.58, 8.05, -9.75, 7.84, -1.32]

# Given and known constants
P = 1.0 # atm
m = 3.0 # g
VFR_in = 0.85 # L/min
T = 400. + 273.15 # K
K = 12.2

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

# Allocate storage for the responses
resp = np.ones(len(yAin)) * -1

# Define the response function
def response_function(inputs, par1, par2, par3, par4, par5):
    # extract the rate coefficient
    k_test = 10.0**par1
    KA = 10.0**par2
    KB = 10.0**par3
    KY = 10.0**par4
    KZ = 10.0**par5

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
#par_guess = [0.0, -10.0, -10.0, -10.0, -10.0] # gave r_squared = 0.615
#par_guess = [2.0, 3.0, -10.0, -1.0, -10.0] # gave r_squared = 0.981
#par_guess = [3.1, -5.0, -5.0, -5.0, -5.0] # gave r_sauared = 0.995
#par_guess = [5.58, -9.75, 8.05, 7.84, -1.32] # gave r_sauared = 0.995

#par_guess = [1.0, 0.0, 0.0, 0.0, 0.0] # gave r_sauared = -0.138
#par_guess = [3.0, 0.0, 0.0, 0.0, 0.0] # gave r_sauared = -0.138
#par_guess = [5.0, 2.66, 3.67, 2.97, -2.5] # gave r_sauared = -0.138
#par_guess = [6.0, 2.66, 3.67, 2.97, -2.5] # gave r_sauared = -0.138

# calculate the model-predicted responses
y_model = response_function(adjusted_inputs,par_guess[0],par_guess[1],\
        par_guess[2],par_guess[3],par_guess[4])

# calculate the residuals
residuals = PA - y_model

# create a parity plot
plt.figure(1) 
plt.plot(PA, y_model, color = 'r', marker='o', ls='')
plt.plot([np.min(PA), np.max(PA)], [np.min(PA), np.max(PA)], color = 'k')
plt.xlabel("experimental response (atm)")
plt.ylabel("model-predicted response (atm)")

# show the parity plot
plt.show()

# create a residuals plot for the inlet yA
plt.figure(2) 
plt.plot(yAin, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Inlet A Mole Fraction")
plt.ylabel("Residual (%)")

# show the residuals plot for the reaction time
#plt.show()

# create a residuals plot for the inlet yB
plt.figure(3) 
plt.plot(yBin, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Inlet B Mole Fraction")
plt.ylabel("Residual (%)")

# show the residuals plot for the reaction time
#plt.show()

# create a residuals plot for the inlet yY
plt.figure(4) 
plt.plot(yYin, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Inlet Y Mole Fraction")
plt.ylabel("Residual (%)")

# show the residuals plot for the reaction time
#plt.show()

# create a residuals plot for the inlet yZ
plt.figure(5) 
plt.plot(yZin, residuals, color = 'r', marker='o', ls='')
plt.axhline(y=0, color = 'k')
plt.xlabel("Inlet Z Mole Fraction")
plt.ylabel("Residual (%)")

# show the residuals plot for the reaction time
#plt.show()