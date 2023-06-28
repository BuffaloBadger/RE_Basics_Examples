import response_function as resp
import pandas as pd
import matplotlib.pyplot as plt

# Path to directory for saving files
filepath = './reb_13_1/Results/'

# Calculate the response
t, CA, CB, CY, CZ, T, r = resp.response()

# Save the response
df = pd.DataFrame()
df['t'] = t.tolist()
df['CA'] = CA.tolist()
df['CB'] = CB.tolist()
df['CY'] = CY.tolist()
df['CZ'] = CZ.tolist()
df['T'] = T.tolist()
df['r'] = r.tolist()
df.to_csv(filepath+'reb_13_1_response.csv')

# Process the results
plt.figure(1) # concentrations vs time
plt.plot(t, CA, color = 'tab:orange', label = 'A')
plt.plot(t, CB, color = 'tab:blue', label = 'B')
plt.plot(t, CY, color = 'tab:red', lw = 4.0, ls = 'dashed', label = 'Y')
plt.plot(t, CZ, color = 'tab:green', lw = 2.0, label = 'Z')
plt.xlabel("Time (h)")
plt.ylabel("Cooncentration (M)")
plt.legend()

# save and show the figure
plt.savefig(filepath + 'reb_13_1_C_vs_t.png')
plt.show()

plt.figure(2) # temperature vs time
plt.plot(t, T, color = 'k')
plt.xlabel("Time (h)")
plt.ylabel("Temperature (°C)")

# save and show the figure
plt.savefig(filepath + 'reb_13_1_T_vs_t.png')
plt.show()

# plot rate vs. t
plt.figure(3) # rate vs time
plt.plot(t, r, color = 'k')
plt.xlabel("Time (h)")
plt.ylabel("Rate (mol L$^{-1}$ h$^{-1}$)")

# save and show the figure
plt.savefig(filepath + 'reb_13_1_r_vs_t.png')
plt.show()

plt.figure(1) # temperature at early time
plt.plot(t[0:20], T[0:20], color = 'k')
plt.xlabel("Time (h)")
plt.ylabel("Temperature (°C)")

# save and show the figure
#plt.savefig(filepath + 'reb_13_1_C_vs_t_early.png')
plt.show()
