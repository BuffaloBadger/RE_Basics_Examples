import matplotlib.pyplot as plt
import numpy as np

# set filepaths
path_to_data = './reb_4_3/Data/'
path_to_results = './reb_4_3/Results/'
path_to_figures = '../RE_Basics/Graphics/'

# Given or known
mu_max = 1.0 # /h
Ks = 0.2 # g/L

# Create vectors for plot data
C_S = np.linspace(0.0, 5.0, 100)
mu = np.zeros(100, dtype=float, order='C')

# Calculate rates for each conversion
for i in range(0, np.size(C_S)) :
    # calculate the rate
    mu[i] = mu_max*C_S[i]/(Ks + C_S[i])

# Plot the results
plt.figure(1)
plt.plot(C_S, mu, color = 'b')
plt.xlabel("Substrate Mass Concentration (g L$^{-1}$)")
plt.ylabel("Specific Growth Rate (h$^{-1}$)")
# save and show the figure
plt.savefig(path_to_results + 'reb_4_3_fig_1.png')
plt.savefig(path_to_figures + 'reb_4_3_fig_1.png')
plt.show()
