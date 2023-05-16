import numpy as np
import pandas as pd
import scipy as sp
import random

# set the filename for the results
filename1 = 'reb_10_1/Data/reb_10_1_data.csv'
filename2 = '../RE_Basics/Data/reb_10_1_data.csv'

# Constant inputs
V = 0.1 # L

# kinetics parameters
k = 5.29E9*np.exp(-12100/1.987/(35 + 273.15))
print('The rate coefficient used to generate these data was %.2g '%(k) + \
        'M^-0.5^ min')

# Define rate expression
def rate(concA, concB, concZ):
    return k*concA*concB/concZ**2

# Adjusted inputs
VFR = np.array([50, 100]) # cc/min
VFR = VFR/1000.0 # L/min
CA0 = np.array([0.5, 1.0, 2.5, 5]) # M
CB0 = np.array([0.5, 1.0, 2.5, 5]) # M
CY0 = np.array([0.5, 1.0, 2.5, 5]) # M
CZ0 = np.array([0.5, 1.0, 2.5, 5]) # M

# Create empty dataframe for the results
df = pd.DataFrame(columns=['VFR', 'CAin', 'CBin', 'CYin', 'CZin', 'CYout'])

# Calculate the responses
for Vdot in VFR:
    for CAin in CA0:
        for CBin in CB0:
            for CYin in CY0:
                for CZin in CZ0:
                    # define the mole balance
                    def mole_balance(CYout):
                        extent = CYout - CYin
                        CA = CAin - extent
                        CB = CBin - extent
                        CZ = CZin + extent
                        r = rate(CA, CB, CZ)
                        return (CAin - CA)*Vdot - V*r
                    
                    # solve the mole balance
                    guess = [CYin + 0.5]
                    soln = sp.optimize.root(mole_balance, guess)
                    if not(soln.success):
                        print("A solution was NOT obtained:")
                        print(soln.message)
                    CY = soln.x[0]

                    # add +/- 0.01 random "error"
                    random_error = (2*random.random() - 1.0)*0.01
                    CY += random_error

                    # round to 2 decimal places
                    CY = round(CY,2)

                    # append the result to the dataframe
                    df.loc[len(df.index)] = [1000*Vdot,CAin,CBin,CYin,CZin,CY]

# display the results
print("\n")
print(df)

# save the results
print("\nSaving results to " + filename1 + "\n")
df.to_csv(filename1,index=False)

print("\nSaving results to " + filename2 + "\n")
df.to_csv(filename2,index=False)
