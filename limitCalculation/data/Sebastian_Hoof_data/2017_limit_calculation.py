import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from scipy.interpolate import interp1d

# Load the data from the file
lgm, lgg, ll = np.genfromtxt("cast2017_loglike.dat", unpack=True)

# Define the target mass and find the closest available mass in the dataset
target_mass = 0.0001
closest_mass_index = np.argmin(np.abs(lgm - target_mass))
closest_mass_value = lgm[closest_mass_index]

# Extract the corresponding likelihood and g values for this closest mass
mask_closest_mass = np.isclose(lgm, closest_mass_value)
selected_likelihood = np.exp(ll[mask_closest_mass] - ll[mask_closest_mass].max())
selected_g = lgg[mask_closest_mass]
selected_g4 = selected_g**4

# Compute the total area under the likelihood curve using Simpson's rule
area_g = simps(selected_likelihood, selected_g)
area_g4 = simps(selected_likelihood, selected_g4)

# Compute the cumulative integral using Simpson's rule for g and g^4 values
cumulative_integral_g = np.array([simps(selected_likelihood[:i+1], selected_g[:i+1]) for i in range(selected_g.shape[0])])
cumulative_integral_g4 = np.array([simps(selected_likelihood[:i+1], selected_g4[:i+1]) for i in range(selected_g4.shape[0])])

# Total areas using Simpson's rule
total_area_g = cumulative_integral_g[-1]
total_area_g4 = cumulative_integral_g4[-1]

# Calculate 95% of the total area
area_95_percent_g = 0.95 * total_area_g
area_95_percent_g4 = 0.95 * total_area_g4

# Create an interpolation function for the cumulative integral
interp_cumulative_integral_g = interp1d(cumulative_integral_g, selected_g, kind='linear', fill_value="extrapolate")
interp_cumulative_integral_g4 = interp1d(cumulative_integral_g4, selected_g4, kind='linear', fill_value="extrapolate")

# Find the value of g and g^4 corresponding to 95% of the total area
g_value_95_percent = interp_cumulative_integral_g(area_95_percent_g)
g4_value_95_percent = interp_cumulative_integral_g4(area_95_percent_g4)
g_value_based_on_g4 = g4_value_95_percent**0.25
"""
# Plotting for g
plt.figure(figsize=(12, 7))
plt.plot(selected_g, selected_likelihood, '-o', markersize=3, label='Using $g$')
plt.axvline(g_value_95_percent, color='red', linestyle='--', label=f'95% value: $g = {g_value_95_percent:.2e}$')
plt.xlabel('ALP-photon coupling $g$')
plt.ylabel('Likelihood')
plt.title(f'Likelihood vs ALP-photon coupling $g$ for ALP mass {closest_mass_value}')
plt.legend()
plt.xlim([0, 1e-10])
plt.grid(True)
plt.tight_layout()
plt.show()
"""
# Plotting for g^4
plt.figure(figsize=(12, 7))
g4_xlim_min = 0
g4_xlim_max = 1e-40
plt.plot(selected_g4, selected_likelihood, '-o', markersize=3)
plt.axvline(g4_value_95_percent, color='red', linestyle='--', label=f'95% value: $g^4 = {g4_value_95_percent:.2e}$ (corresponding $g$ value: {g4_value_95_percent**0.25:.2e})')
plt.xlabel('ALP-photon coupling $g^4$')
plt.ylabel('Likelihood')
plt.title(f'Likelihood vs ALP-photon coupling $g^4$ for ALP mass {closest_mass_value}')
plt.legend(fontsize=12)
plt.xlim([g4_xlim_min, g4_xlim_max])
plt.grid(True)
plt.savefig("Hoof_data_limit.pdf")
plt.show()


# Print the results
print(f"For the ALP mass value closest to {target_mass} eV:")
print(f"Total area under the likelihood curve for g: {area_g:.2e}")
print(f"Total area under the likelihood curve for g^4: {area_g4:.2e}")
print(f"Total area method 2 under the likelihood curve for g: {total_area_g:.2e}")
print(f"Total area method 2 under the likelihood curve for g^4: {total_area_g4:.2e}")
print(f"95% of the area under the likelihood curve for g: {area_95_percent_g:.2e}")
print(f"95% of the area under the likelihood curve for g^4: {area_95_percent_g4:.2e}")
print(f"g value leaving 95% of the area under the likelihood curve to the left: {g_value_95_percent:.2e}")
print(f"g4 value leaving 95% of the area under the likelihood curve to the left: {g4_value_95_percent:.2e}")
print(f"g value as the 4th root of g4 leaving 95% of the area under the likelihood curve to the left: {g_value_based_on_g4:.2e}")

