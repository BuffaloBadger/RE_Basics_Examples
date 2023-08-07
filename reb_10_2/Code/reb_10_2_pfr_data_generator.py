import numpy as np
import pandas as pd
import scipy as sp
import rebutils as reb
import random

# Constant inputs
m = 3.0 # g
VFRin = 0.85 # L/min @stp
P = 1.0 # atm
T = 400 + 273.15 # K
K = 12.2

# kinetics parameters
k = 0.058/60.0 # mol/g/min

# Define rate expression
def rate(P_A, P_B, P_Y, P_Z):
    reaction_rate = k*P_A**0.9*P_B**0.25*P_Y**-0.6*(1 - P_Y*P_Z/K/P_A/P_B)
    return reaction_rate

# Adjusted inputs
y1 = np.array([0.3, 0.4, 0.5, 0.6, 0.7])
y2 = np.array([0.33, 0.5, 0.66])

# Create empty dataframe for the results
df = pd.DataFrame(columns=['yA','yB','yY','yZ','PA'])

count = 0

# Calculate the responses
for yA in y1:
    for yb in y2:
        yB = yb*(1-yA)
        for yc in y2:
            yY = yc*(1-yA-yB)
            yZ = 1 - yA - yB - yY

            equil_factor = yY*yZ/K/yA/yB
            if equil_factor > 1.0:
                count += 1
    
            # Define the mole balances
            def mole_balances(t,n):
                n_total = np.sum(n)
                PA = n[0]/n_total*P
                PB = n[1]/n_total*P
                PY = n[2]/n_total*P
                PZ = n[3]/n_total*P
                r = rate(PA, PB, PY, PZ)
                ddm = np.array([-r, -r, r, r])
                return ddm
            
            # Solve the mole balances
            m0 = 0
            n0 = np.array([yA, yB, yY, yZ]) * VFRin/22.4
            f_var = 0
            f_val = m
            soln = reb.solveIVODEs(m0, n0, f_var, f_val, mole_balances)

            # calculate the response
            nA = soln.y[0,-1]
            nB = soln.y[1,-1]
            nY = soln.y[2,-1]
            nZ = soln.y[3,-1]
            n_total = nA + nB + nY + nZ
            PA = nA/n_total*P

            # add +/- 0.01 random "error"
            random_error = (2*random.random() - 1.0)*0.01
            PA += random_error

            # round to 3 decimal places
            yA = round(yA,3)
            yB = round(yB,3)
            yY = round(yY,3)
            yZ = round(yZ,3)
            PA = round(PA,3)

            # append the result to the dataframe
            df.loc[len(df.index)] = [yA, yB, yY, yZ, PA]

print('There were %i experiments that started beyond equilibrium' %(count))

# display the results
print("\n")
print(df)

# save the results
filename = "reb_10_2/Data/reb_10_2_data.csv"
print("\nSaving results to " + filename + "\n")
df.to_csv(filename,index=False)

filename = "../RE_Basics/Data/reb_10_2_data.csv"
print("\nSaving results to " + filename + "\n")
df.to_csv(filename,index=False)
