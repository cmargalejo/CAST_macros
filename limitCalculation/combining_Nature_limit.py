from unbinned_likelihood_CAST_g4_datasets import *

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from scipy.interpolate import interp1d

# Load the data from the file
lgm, lgg, ll = np.genfromtxt("data/Sebastian_Hoof_data/cast2017_loglike.dat", unpack=True)

# Define the target mass and find the closest available mass in the dataset
target_mass = 0.0001
closest_mass_index = np.argmin(np.abs(lgm - target_mass))
closest_mass_value = lgm[closest_mass_index]

# Extract the corresponding likelihood and g values for this closest mass
mask_closest_mass = np.isclose(lgm, closest_mass_value)
selected_likelihood = np.exp(ll[mask_closest_mass] - ll[mask_closest_mass].max())
selected_g = lgg[mask_closest_mass]
selected_g4 = selected_g**4

# Perform linear interpolation to be able to get the likelihood value at any g
interp_likelihood_g = interp1d(selected_g, selected_likelihood, kind='linear', fill_value="extrapolate")
interp_likelihood_g4 = interp1d(selected_g4, selected_likelihood, kind='linear', fill_value="extrapolate")

# Test the interpolation at a specific value of g and g^4
test_g_value = 0.66e-10
test_g4_value = test_g_value**4
likelihood_at_test_g = interp_likelihood_g(test_g_value)
likelihood_at_test_g4 = interp_likelihood_g4(test_g4_value)

print("Likelihood for g = ", test_g_value, ", is ", likelihood_at_test_g)
print("Likelihood for g4 = ", test_g4_value, ", is ", likelihood_at_test_g4)
#test_g_value, likelihood_at_test_g, test_g4_value, likelihood_at_test_g4

# Plot the interpolation for g
plt.figure(figsize=(12, 6))
g_xlim_min = 0
g_xlim_max = 1e-10
plt.plot(selected_g, selected_likelihood, 'o', label='Original Data')
plt.plot(test_g_value, likelihood_at_test_g, 'r*', markersize=10, label=f'Interpolated Likelihood at $g={test_g_value}$')
plt.axvline(test_g_value, color='red', linestyle='--', label=f'Test g value = {test_g_value:.2e}$')
plt.title('Interpolation of Likelihood for $g$')
plt.xlabel('ALP-photon coupling $g$')
plt.ylabel('Likelihood')
plt.xlim([g_xlim_min, g_xlim_max])
plt.legend()
plt.grid(True)
plt.show()

# Plot the interpolation for g^4
plt.figure(figsize=(12, 6))
g4_xlim_min = 0
g4_xlim_max = 1e-40
plt.plot(selected_g4, selected_likelihood, 'o', label='Original Data')
plt.plot(test_g4_value, likelihood_at_test_g4, 'r*', markersize=10, label=f'Interpolated Likelihood at $g^4={test_g4_value}$')
plt.axvline(test_g4_value, color='red', linestyle='--', label=f'Test g4 value = {test_g4_value:.2e}$')
plt.axhline(likelihood_at_test_g4, color='red', linestyle='--', label=f'Likelihood at test g4 value = {likelihood_at_test_g4:.2e}$')
plt.title('Interpolation of Likelihood for $g^4$')
plt.xlabel('ALP-photon coupling $g^4$')
plt.ylabel('Likelihood')
plt.xlim([g4_xlim_min, g4_xlim_max])
plt.legend()
plt.grid(True)
plt.show()


# Define the range of g^4 values for which you want to calculate the combined likelihood
#g_aγ4_values = np.linspace(start=np.min(selected_g), stop=np.max(selected_g), num=1000)
g_aγ4_values = np.linspace(start=0, stop=1e-40, num=10)

# Calculate the combined likelihood for each g^4 value
combined_likelihoods = [totalLikelihood(dataset_Ar1, dataset_Ar2, dataset_Xe, g_aγ4, likelihood) * interp_likelihood_g4(g_aγ4) for g_aγ4 in g_aγ4_values]

# Plotting
plt.figure(figsize=(14, 7))
# Plot the total likelihood from the datasets
plt.plot(g_aγ4_values, totalLikelihood, label='Total Likelihood from Datasets', color='blue')
# Plot the interpolated likelihood
plt.plot(g_aγ4_values, interp_likelihood_g4, label='Interpolated Likelihood', color='orange')
# Plot the product of the two likelihoods
plt.plot(g_aγ4_values, combined_likelihoods, label='Product of Likelihoods', color='green')
plt.xlabel('g^4 values')
plt.ylabel('Likelihood')
plt.title('Likelihood Comparisons')
plt.legend()
plt.grid(True)
plt.show()
