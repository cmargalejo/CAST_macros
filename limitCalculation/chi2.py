#from unbinned_likelihood_CAST_g4_datasets import *
from limitLucaLista import *
"""
g_aγs = np.linspace(-3.0e-40, 3.0e-40, 1000)
logL2 = [totalLogLikelihood(dataset_Ar1, dataset_Ar2, dataset_Xe, g_aγ, chi2) for g_aγ in g_aγs]
plt.plot(g_aγs, logL2)
# Calculate chi2 at g=0 and find the chi2 target for the 95% confidence level
chi2_at_g0 = totalLogLikelihood(dataset_Ar1, dataset_Ar2, dataset_Xe, 0, chi2)
chi2_target = chi2_at_g0 + 3.84
plt.xlabel(' $g_{aγ}^4$ (GeV$^{-4}$)')
plt.ylabel('-2log L or Chi$^2$')
#plt.show()
#plt.savefig(f"ChiSquareg4_Cris_approach_all_datasets_zoom.pdf")
#plt.close()
#plt.savefig(f"ChiSquareg4_Nature_approach_all_datasets.pdf")
#findMinimumMinusLogLikelihood(chi2, dataset, g_aγs)

# Find the value of g_aγ at the chi2 target
index = (np.abs(logL2 - chi2_target)).argmin()
g_aγ_at_95CL = g_aγs[index]
g = g_aγ_at_95CL**0.25

# Plot the horizontal line at chi2_target and vertical line at g
plt.axhline(y=chi2_target, color='r', linestyle='--', label='95% CL Target')
plt.axvline(x=g_aγ_at_95CL, color='g', linestyle='--', label='g_aγ at 95% CL')
# Annotate the limit value on the plot
#plt.text(g_aγ_at_95CL, chi2_target, f'  g = {g:.2e}', verticalalignment='bottom', horizontalalignment='right')

# Adding labels and legend
plt.xlabel('$g_{a\gamma}^4$ (GeV$^{-4}$)')
plt.ylabel('-2log L or Chi$^2$')
plt.legend()
# Show the plot
plt.show()
# Print the results
print("chi2 at g=0 is: ", chi2_at_g0)
print("Limit at g4 = ", g_aγ_at_95CL, " or g = ", g)
print("chi2 at g=0+3.84 is: ", chi2_target)
print("Limit at g4 = ", g_aγ_at_95CL, " or g = ", g)
"""
"""
# Given your code snippet and assuming you have already defined the totalLogLikelihood function
# and the datasets (dataset_Ar1, dataset_Ar2, dataset_Xe), we can proceed to find the value on the x-axis 
# that corresponds to the chi2_target.

# Generate g_aγ values and corresponding log likelihoods (chi2 values)
g_aγs = np.linspace(-3.0e-40, 3.0e-40, 1000)
logL2 = [totalLogLikelihood(dataset_Ar1, dataset_Ar2, dataset_Xe, g_aγ, chi2) for g_aγ in g_aγs]

# Convert the list to a NumPy array for element-wise operations
logL2_array = np.array(logL2)

# Calculate chi2 at g=0 using the actual function and data
chi2_at_g0 = totalLogLikelihood(dataset_Ar1, dataset_Ar2, dataset_Xe, 0, chi2)

# Calculate the chi2 target for the 95% confidence level
chi2_target = chi2_at_g0 + 3.84

# Find the value of g_aγ at the chi2 target
# Ensure we only consider positive g_aγ values since we are looking for a physical solution
positive_g_aγs = g_aγs[g_aγs >= 0]
positive_logL2 = logL2_array[g_aγs >= 0]

# Find the index of the chi2 value closest to the chi2_target
index = (np.abs(positive_logL2 - chi2_target)).argmin()
g_aγ_at_95CL = positive_g_aγs[index]
g = g_aγ_at_95CL ** 0.25

# Now let's plot the -2 log likelihood (or chi2) versus g_aγ including the lines for chi2_target
plt.plot(g_aγs, logL2, label='Chi-squared')
plt.axhline(y=chi2_target, color='r', linestyle='--', label='95% CL Target')
plt.axvline(x=g_aγ_at_95CL, color='g', linestyle='--', label='g_aγ at 95% CL')
plt.xlabel('$g_{a\gamma}^4$ (GeV$^{-4}$)')
plt.ylabel('-2log L or Chi$^2$')
plt.legend()
plt.show()

# Print the g_aγ value corresponding to the chi2 target
print("The value of g_aγ4 corresponding to the chi2 target is:", g_aγ_at_95CL)
print("The value of g_aγ corresponding to the chi2 target is:", g)

# Find the chi^2 value for g = 1.89e-41
target_g = 1.89e-41 #6.6e-11**4
# Finding the index of the g value that is closest to the target_g
closest_g_index = np.argmin(np.abs(g_aγs - target_g))
# Extracting the corresponding chi^2 value
chi2_value_at_target_g = chi2[closest_g_index]
print("The value of chi2 corresponding to the g=6.6e-41 GeV-1 is:", chi2_value_at_target_g)
diff = chi2_value_at_target_g - chi2_at_g0
print("The difference wrt chi2 when g=0 corresponding to the g=6.6e-41 GeV-1 is:", diff)
"""

# Generate g_aγ values and corresponding log likelihoods (chi2 values)
g_aγs = np.linspace(-3.0e-40, 1.0e-40, 1000)
#logL2 = [totalLogLikelihood(dataset_Ar1, dataset_Ar2, dataset_Xe, g_aγ, chi2) for g_aγ in g_aγs] #chi2 chi2igor chi2basti
logL2 = totalLogLikelihood(dataset_Ar1, dataset_Ar2, dataset_Xe, g_aγs, chi2) # to use with the optimised code
#logL2 = [totalLogLikelihood(dataset_Ar1, dataset_Ar2, dataset_Xe, g_aγ, logLikelihood2) for g_aγ in g_aγs]

# Convert the list to a NumPy array for element-wise operations
logL2_array = np.array(logL2)
print(logL2_array)

# Calculate chi2 at g=0 using the actual function and data
chi2_at_g0 = totalLogLikelihood(dataset_Ar1, dataset_Ar2, dataset_Xe, 0, chi2)
# chi2_at_g0 = totalLogLikelihood(dataset_Ar1, dataset_Ar2, dataset_Xe, 0, logLikelihood2)
"""
# find the minimum chi2 value
chi2_min_index = np.argmin(logL2_array)
chi2_min = logL2_array[chi2_min_index]
g4_min = g_aγs[chi2_min_index]
print("chi2 min = ", chi2_min, ", and g4 at min chi2 = ", g4_min)
"""
# Ignore nan values and find the minimum chi2 value
chi2_min = np.nanmin(logL2_array)
# Find the index of this minimum value
chi2_min_index = np.where(logL2_array == chi2_min)[0][0]  # Get the first occurrence of the minimum value
# Corresponding g4 value at minimum chi2
g4_min = g_aγs[chi2_min_index]
print("Minimum chi2 value ignoring nan:", chi2_min, ", and g4 at this minimum chi2:", g4_min)


# Calculate the chi2 target for the 95% confidence level
chi2_target_0 = chi2_at_g0 + 3.84 # starting from chi2=0
# Filter out nan values from the array
valid_logL2_array = logL2_array[~np.isnan(logL2_array)]
chi2_target = chi2_min + 3.84 # starting from the minimum
# Find the index of the closest value to chi2_target in the valid_logL2_array
closest_index = np.abs(valid_logL2_array - chi2_target).argmin()


# Find the value of g_aγ at the chi2 target
# Ensure we only consider positive g_aγ values since we are looking for a physical solution
positive_g_aγs = g_aγs[g_aγs >= 0]
positive_logL2 = logL2_array[g_aγs >= 0]

# Find the index of the chi2 value closest to the chi2_target
index = (np.abs(positive_logL2 - chi2_target_0)).argmin()
g_aγ_at_95CL = positive_g_aγs[index]
g = g_aγ_at_95CL ** 0.25

# Find the index of the chi2 value closest to the chi2_target
index_at_95CL = (np.abs(logL2_array - chi2_target)).argmin()
print("min index = ", chi2_min_index, ", chi target index = ", index_at_95CL)
# Find the corresponding g4 value at this index
g4_at_95CL = g_aγs[index_at_95CL]
print("95% CL chi2 value ignoring nan:", chi2_target, ", and g4 at this minimum chi2:", g4_at_95CL)

# Find the corresponding g4 value in the original g_aγs array
# The number of nan values needs to be accounted for to get the correct index in g_aγs
num_nans = np.isnan(logL2_array).sum()
g4_at_chi2_target = g_aγs[closest_index + num_nans]

# Printing the found values
print("Closest chi2 to chi2_target:", valid_logL2_array[closest_index])
print("g4 value at this chi2:", g4_at_chi2_target)


# Now let's plot the -2 log likelihood (or chi2) versus g_aγ including the lines for chi2_target
plt.plot(g_aγs, logL2, label='Chi-squared')
plt.axhline(y=chi2_target, color='r', linestyle='--', label='95% CL Target from minimum')
#plt.axhline(y=chi2_target_0, color='r', linestyle='-.', label='95% CL Target from 0')
plt.axvline(x=g4_at_chi2_target, color='g', linestyle='--', label='g_aγ at 95% CL from minumum')
#plt.axvline(x=g_aγ_at_95CL, color='g', linestyle='-.', label='g_aγ at 95% CL from 0')
plt.xlabel('$g_{a\gamma}^4$ (GeV$^{-4}$)')
plt.ylabel('-2log L or Chi$^2$')
plt.legend()
plt.show()

# Print the g_aγ value corresponding to the chi2 target
print("The value of g_aγ4 corresponding to the chi2 95% CL target is:", g_aγ_at_95CL)
print("The value of g_aγ corresponding to the chi2 95% CL target is:", g)

# Compute the chi2 for g^4 = 1.89e-41 or g = 6.6e-11
target_g4 = 1.89e-41
# Calculate the chi2 value for the given g^4
chi2_value_at_target_g4 = totalLogLikelihood(dataset_Ar1, dataset_Ar2, dataset_Xe, target_g4, chi2) #logLikelihood2
print("The chi2 value for g^4 = 1.89e-41 is or g = 6.6e-11 is: ", chi2_value_at_target_g4)
diff = chi2_value_at_target_g4 - chi2_at_g0
print("The value of chi2 at g = 0 is: ", chi2_at_g0, )
print("The difference wrt chi2 when g=0 corresponding to the g=6.6e-11 GeV-1 is:", diff)

"""
Confidence Level (CL)    χ 2 value
90%                         2.71
95%                         3.84
99%                         6.63
99.9%                       10.83
"""
# Now I repeat for likelihood and see what I get...


# I plot all the likelihoods together

# Combining the plotting commands to make a single plot with all the curves

# Define the datasets
datasets = [dataset_Ar1, dataset_Ar2, dataset_Xe]
"""
datasets = [dataset_Ar1, dataset_Ar2, dataset_Xe]
# Plotting the likelihood for each individual dataset
for dataset in datasets:
    g4Lin_individual = np.linspace(0.0, 1e-40, 1000)
    likelihoodLin_individual = likelihood(dataset, g4Lin_individual)
    plt.plot(g4Lin_individual, likelihoodLin_individual, label=f"{dataset.name}")
    #lim_individual = limit(dataset, likelihood)
    #plt.axvline(x=lim_individual, linestyle='--', color='grey')
"""
# Plotting the total likelihood for all datasets combined
g4Lin_total = np.linspace(0, 1e-40, 1000)
likelihoodLin_total = totalLikelihood(dataset_Ar1, dataset_Ar2, dataset_Xe, g4Lin_total, likelihood2)
plt.plot(g4Lin_total, likelihoodLin_total, label='All Datasets Combined', color='black', linewidth=3)
#plt.axvline(x=totalLim, color='r', linewidth=2)
# Adding labels, title, legend, and saving the plot
plt.xlabel("Coupling constant (GeV$^-4$)")
plt.ylabel("Likelihood")
plt.title("Likelihood vs. Coupling Constant for Different Datasets using likelihood")
plt.show()
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig("combined_likelihood_plot_debug.pdf")
#plt.show()
"""
# Find the index and value of the minimum chi2
min_chi2_index = np.argmin(logL2_array)
min_chi2_value = logL2_array[min_chi2_index]
min_g_aγ = g_aγs[min_chi2_index]
print("Minimum chi2 value:", min_chi2_value)
print("g_aγ at minimum chi2:", min_g_aγ)

# Calculate chi2 target (minimum chi2 + 3.84)
chi2_target = min_chi2_value + 3.84

# Find the index of the chi2 value closest to the chi2_target
target_index = (np.abs(logL2_array - chi2_target)).argmin()
target_g_aγ = g_aγs[target_index]
target_g = target_g_aγ ** 0.25

# Plotting
plt.plot(g_aγs, logL2_array, label='Chi-squared')
plt.axhline(y=chi2_target, color='r', linestyle='--', label='Chi2 Target')
plt.axvline(x=target_g_aγ, color='g', linestyle='--', label='Target g_aγ')
plt.xlabel('$g_{a\gamma}^4$ (GeV$^{-4}$)')
plt.ylabel('Chi$^2$')
plt.legend()
plt.show()

# Print the results
print("Minimum chi2 value:", min_chi2_value)
print("g_aγ at minimum chi2:", min_g_aγ)
print("The value of g_aγ corresponding to the chi2 target is:", target_g_aγ)
print("The value of g corresponding to the chi2 target is:", target_g)
"""
for dataset in datasets:
    g4Lin_individual = np.linspace(-1e-39, 1e-40, 1000)
    likelihoodLin_individual = likelihood(dataset, g4Lin_individual)
    plt.plot(g4Lin_individual, likelihoodLin_individual, label=f"{dataset.name}")
    plt.xlabel("Coupling constant (GeV$^-4$)")
    plt.ylabel("Likelihood")
    plt.title("Likelihood vs. Coupling Constant for Different Datasets using likelihood")
    plt.legend()
    plt.show()

for dataset in datasets:
    g4Lin_individual = np.linspace(-1e-39, 1e-40, 1000)
    chi2_individual = chi2(dataset, g4Lin_individual)
    plt.plot(g4Lin_individual, chi2_individual, label=f"{dataset.name}")
    plt.xlabel("Coupling constant (GeV$^-4$)")
    plt.ylabel("-2log L or Chi$^2$")
    plt.title("Chi2 vs. Coupling Constant for Different Datasets using likelihood")
    plt.legend()
    plt.show()