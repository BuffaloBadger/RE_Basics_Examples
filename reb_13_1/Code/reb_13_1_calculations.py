import rebutils as reb
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Path to directory for saving files
filepath = './reb_12_1/python/'

# Given and known in consistent units
dH = -41200.0 # J/mol
k0 = 51100.0 * 3600 # L/mol/h
E = 74800.0 # J/mol
T0 = 180 + 273.15 # K
CA0 = 2.9 # mol/L
CB0 = 3.2 # mol/L
CY0 = 0.0 # mol/L
CZ0 = 0.0 # mol/L
Cp = 1.23 * 4.184 # J/g/K
rho = 1.02E3 # g/L
t_rxn = 2 # h
V = 1.0 # L (basis)
R = 8.314 # J/mol/K

# Define the response function
def response():
    # Define the reactor model
    def reactor_model(time,dep_vars):
        # evaluate the derivatives at the current values of the variables
        nA_cur = dep_vars[0]
        nB_cur = dep_vars[1]
        T_cur = dep_vars[4]
        k = k0*np.exp(-E/R/T_cur)
        CA_cur = nA_cur/V
        CB_cur = nB_cur/V
        r = k*CA_cur*CB_cur
        ddt = np.array([-V*r, -V*r, V*r, V*r, -r*dH/rho/Cp])
        # Return the derivatives
        return ddt

    # Solve the reactor model equations
    t0 = 0
    dep0 = np.array([CA0*V, CB0*V, CY0*V, CZ0*V, T0])
    f_var = 0
    f_val = t_rxn
    soln = reb.solveIVODEs(t0, dep0, f_var, f_val, reactor_model, \
            odes_are_stiff=False)
    if not(soln.success):
        print("A solution was NOT obtained:")
        print(soln.message)

    # Extract the results
    t = soln.t
    nA = soln.y[0,]
    nB = soln.y[1,]
    nY = soln.y[2,]
    nZ = soln.y[3,]
    T = soln.y[4,]

    # Calculate the response
    CA = nA/V
    CB = nB/V
    CY = nY/V
    CZ = nZ/V
    T = T - 273.15

    # Return the response
    return t, CA, CB, CY, CZ, T

# Calculate missing input to the response function

# Calculate the response
t, CA, CB, CY, CZ, T = response()

# Save the response
df = pd.DataFrame()
df['t'] = t.tolist()
df['CA'] = CA.tolist()
df['CB'] = CB.tolist()
df['CY'] = CY.tolist()
df['CZ'] = CZ.tolist()
df['T'] = T.tolist()
df.to_csv(filepath+'reb_12_1_response.csv')

# Process the resulst
plt.figure(1) # concentrations vs time
plt.plot(t, CA, color = 'tab:orange', label = 'A')
plt.plot(t, CB, color = 'tab:blue', label = 'B')
plt.plot(t, CY, color = 'tab:red', lw = 4.0, ls = 'dashed', label = 'Y')
plt.plot(t, CZ, color = 'tab:green', lw = 2.0, label = 'Z')
plt.xlabel("Time (h)")
plt.ylabel("Cooncentration (M)")
plt.legend()

# save and show the figure
plt.savefig(filepath + 'reb_12_1_C_vs_t.png')
plt.show()

plt.figure(2) # temperature vs time
plt.plot(t, T, color = 'k')
plt.xlabel("Time (h)")
plt.ylabel("Temperature (°C)")

# save and show the figure
plt.savefig(filepath + 'reb_12_1_T_vs_t.png')
plt.show()

# plot rate vs. t
rate = k0*np.exp(-E/R/T)*CA*CB
plt.figure(3) # rate vs time
plt.plot(t, rate, color = 'k')
plt.xlabel("Time (h)")
plt.ylabel("Rate (mol L$^{-1}$ h$^{-1}$)")

# save and show the figure
plt.savefig(filepath + 'reb_12_1_r_vs_t.png')
plt.show()
