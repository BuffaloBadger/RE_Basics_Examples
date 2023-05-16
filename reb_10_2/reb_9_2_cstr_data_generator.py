import numpy as np
import pandas as pd
import scipy as sp
import random

# set the filename for the results
filename = 'reb_9_2/reb_9_2_data.csv'

# Constant inputs
V = 1.0 # L (basis)
P = 3.0 # atm
R = 0.08206 # L*atm/mol/K
T = 450 + 273.15

# kinetics parameters
k = 0.1556/(R*T)**2
print('The rate coefficient used to generate these data was %.2g '%(k) + \
        'mol/L/s/atm**2')

# Define rate expression
def rate(Apressure,Bpressure):
    return k*Apressure**1.4*Bpressure**0.6

# Adjusted inputs
tau = np.array([30, 40, 50, 60]) # s
yA0 = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])

# Create empty dataframe for the results
df = pd.DataFrame(columns=['tau', 'yAin', 'CZout'])

# Calculate the responses
for space_time in tau:
    VFRin = V/space_time
    for yAin in yA0:
        nAin = yAin*P*VFRin/R/T
        nBin = (1.0 - yAin)*P*VFRin/R/T
        # define the mole balance
        def mole_balance(nZ):
            nA = nAin - nZ
            nB = nBin - nZ
            n_total = nA + nB + nZ
            PA = nA/n_total*P
            PB = nB/n_total*P
            r = rate(PA,PB)
            return nAin - nA - V*r
        
        # solve the mole balance
        guess = [0.01*nAin]
        soln = sp.optimize.root(mole_balance, guess)
        if not(soln.success):
            print("A solution was NOT obtained:")
            print(soln.message)
        nZ = soln.x[0]

        # calculate the response
        n_tot = nAin + nBin - nZ
        VFR = n_tot*R*T/P
        CZ = nZ/VFR*1000 # mmol/L

        # add +/- 0.01 random "error"
        random_error = (2*random.random() - 1.0)*0.2
        CZ += random_error

        # round to 2 decimal places
        CZ = round(CZ,2)

        # append the result to the dataframe
        df.loc[len(df.index)] = [space_time, yAin, CZ]

# display the results
print("\n")
print(df)

# save the results
print("\nSaving results to " + filename + "\n")
df.to_csv(filename,index=False)
