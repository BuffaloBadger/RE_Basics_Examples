import numpy as np
import pandas as pd
import scipy as sp
import random

# Constant inputs
V = 1.0 # L (basis)
P = 3.0 # atm
R = 0.08206 # L*atm/mol/K
Tbase = 450 + 273.15

# kinetics parameters
k = 4.4e-5 # mol/L/atm^2/s
E = 27500. # cal /mol
k0 = k*np.exp(E/1.987/Tbase)
EkJ = E/1000*4.184
print(' ')
print(f'pre-exponential: {k0:.3g}')
print(f'activation energy: {EkJ:.3g}')

# Define rate expression
def rate(k, Apressure,Bpressure):
    return k*Apressure**1.4*Bpressure**0.6

# Adjusted inputs
T_levels = np.array([440, 450, 460]) + 273.15
tau = np.array([30, 40, 50, 60]) # s
yA0 = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])

# Create empty dataframe for the results
df = pd.DataFrame(columns=['T (°C)', 'tau (s)', 'yA_0', 'CZ_1 (mmol/L)'])

# Calculate the responses
for T in T_levels:
    k = k0*np.exp(-E/1.987/T)
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
                r = rate(k,PA,PB)
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
            df.loc[len(df.index)] = [T-273.15,space_time, yAin, CZ]

# display the results
print("\n")
print(df)

# save the results
df.to_csv('reb_20_5_2/reb_20_5_2_data.csv',index=False)
