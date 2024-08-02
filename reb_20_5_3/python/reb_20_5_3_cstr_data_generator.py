import numpy as np
import pandas as pd
import scipy as sp
import random

# Constant inputs
V = 3.0 # gal

# kinetics parameters
kf = 15.98 # gal/lbmol/min
kr = 5.02 # gal/lbmol/min

# Define rate expression
def rate(concA, concY, concZ):
    return kf*concA**2 - kr*concY*concZ

# Adjusted inputs
VFR = np.array([0.5, 1.0, 1.5, 2.0, 2.5]) # gal/min
CA0 = np.array([0.01, 0.02, 0.03]) # lbmol/gal
CY0 = np.array([0.0, 0.01, 0.02, 0.03]) # lbmol/gal
CZ0 = np.array([0.0, 0.01, 0.02, 0.03]) # lbmol/gal

# Create empty dataframe for the results
df = pd.DataFrame(columns=['Vdot', 'CA_0', 'CY_0', 'CZ_0', 'CA_1'])

# Calculate the responses
for Vdot in VFR:
    for CAin in CA0:
        for CYin in CY0:
            for CZin in CZ0:
                # define the mole balance
                def mole_balance(fA):
                    nAin = CAin*Vdot
                    nYin = CYin*Vdot
                    nZin = CZin*Vdot
                    nA = nAin*(1-fA)
                    nY = nYin + nAin*fA
                    nZ = nZin + nAin*fA
                    CA = nA/Vdot
                    CY = nY/Vdot
                    CZ = nZ/Vdot
                    r = rate(CA, CY, CZ)
                    return nAin - nA - V*r
                
                # solve the mole balance
                guess = [0.1]
                soln = sp.optimize.root(mole_balance, guess)
                if not(soln.success):
                    print("A solution was NOT obtained:")
                    print(soln.message)
                fA = soln.x[0]

                # add +/- 0.01 random "error"
                random_error = (2*random.random() - 1.0)*0.05
                fA += random_error

                CAout = CAin*(1-fA)

                # round to 2 decimal places
                CAout = round(CAout,4)

                # append the result to the dataframe
                df.loc[len(df.index)] = [Vdot,CAin,CYin,CZin,CAout]

# display the results
print('\n')
print(df)

# save the results
df.to_csv('reb_20_5_3/reb_20_5_3_data.csv',index=False)