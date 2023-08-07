import pandas as pd
import numpy as np
import rebutils as reb
import matplotlib.pyplot as plt
import response_function as rf

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
y_model = rf.response_function(adjusted_inputs,par_guess[0],par_guess[1],\
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