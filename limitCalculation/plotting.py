from unbinned_likelihood_CAST_g4_datasets import *
from scipy.signal import convolve


def signalIntegrated(dataset, E, g_aγ4): # in counts keV⁻¹ cm^-2
    ## Returns the axion flux based on `g` and energy `E`
    return solarAxionFlux(E) * dataset.total_time * areaBore * conversionProbability() * telescopeEff(E) * detectorEff(dataset, E) * softwareEff(dataset, E) * g_aγ4 
# My signal, which is the expected solar axion flux corrected by my efficiencies
g_aγ4 = 1e-40
E = np.linspace(0, 10, 200) 
datasets = [dataset_Ar1, dataset_Ar2, dataset_Xe]
signal_total = np.zeros_like(E)
for dataset in datasets:
    # Calculate the signal for the current dataset across all energy values
    signal_individual = np.array([signalIntegrated(dataset, e, g_aγ4) for e in E])
    # Plot the individual dataset signal
    plt.plot(E, signal_individual, label=f"{dataset.name}")
    # Accumulate the total signal
    signal_total += signal_individual
plt.plot(E, signal_total, label='Total', color='black', linestyle='--')
plt.xlabel('Energy (keV)')
plt.ylabel('Received solar axion flux (keV⁻¹cm⁻²s⁻¹)')
plt.legend()
plt.show()
plt.savefig("plots/SolarAxionTheoreticalSignalFluxCombined.pdf")
plt.close()

"""
# Define a Gaussian function for the energy resolution
def gaussian(x, sigma):
    return 1/(sigma * np.sqrt(2 * np.pi)) * np.exp(-x**2 / (2 * sigma**2))
energy_resolution = 0.5  # in keV
# Create an array of energy values for the Gaussian
# It should be large enough to cover the range where the Gaussian is significantly non-zero
gaussian_x = np.linspace(-3 * energy_resolution, 3 * energy_resolution, 100)
gaussian_y = gaussian(gaussian_x, energy_resolution)
datasets = [dataset_Ar1, dataset_Ar2, dataset_Xe]
for dataset in datasets:
    # Calculate the signal for the current dataset across all energy values
    #signal_individual = np.array([signalIntegrated(dataset, g_aγ4, e) for e in E])    
    # Convolve the individual signal with the Gaussian kernel
    signal_blurred_individual = convolve(signal_individual, gaussian_y, mode='same')
    # Normalize the convolved signal
    signal_blurred_individual /= np.sum(gaussian_y)
    # Plot the blurred individual dataset signal
    plt.plot(E, signal_blurred_individual, label=f"Blurred {dataset.name}")
signal_blurred_total = convolve(signal_total, gaussian_y, mode='same')
# Normalize the convolved signal
signal_blurred_total /= np.sum(gaussian_y)
plt.plot(E, signal_blurred_total, label='Blurred Total', color='black', linestyle='--')
plt.xlabel('Energy (keV)')
plt.ylabel('Blurred solar axion flux (keV⁻¹cm⁻²s⁻¹)')
plt.legend()
plt.show()
plt.savefig("plots/SolarAxionTheoreticalSignalFluxBlurred.pdf")
plt.close()
"""

"""
# Calculate the step size
dx = gaussian_x[1] - gaussian_x[0]
# Print the Gaussian kernel sum for verification
print(f"Sum of Gaussian kernel: {np.sum(gaussian_y)}")
print(f"Gaussian kernel step size (dx): {dx}")
# Print the maximum value of the original and blurred signals for each dataset
datasets = [dataset_Ar1, dataset_Ar2, dataset_Xe]
for dataset in datasets:
    # Calculate the signal for the current dataset across all energy values
    signal_individual = np.array([signalIntegrated(dataset, g_aγ4, e) for e in E])
    # Convolve the individual signal with the Gaussian kernel
    signal_blurred_individual = convolve(signal_individual, gaussian_y, mode='same')
    # Normalize the convolved signal, accounting for step size
    signal_blurred_individual /= np.sum(gaussian_y) * dx

    # Print the max values for checking
    print(f"{dataset.name} - Max original signal: {np.max(signal_individual)}, Max blurred signal: {np.max(signal_blurred_individual)}")
    # Print diagnostic information
    print(f"{dataset.name} - Max original signal: {np.max(signal_individual)}")
    print(f"{dataset.name} - Sum of blurred signal before normalization: {np.sum(signal_blurred_individual * dx)}")
    print(f"{dataset.name} - Max blurred signal after normalization: {np.max(signal_blurred_individual)}")
"""   
"""
# Now convolve the signal with the Gaussian
# We use 'same' to ensure the convolved signal has the same size as the original signal
signal_blurred = convolve(signal_total, gaussian_y, mode='same')
# Normalization factor to ensure the area under the curve remains the same
normalization_factor = np.sum(gaussian_y) * (gaussian_x[1] - gaussian_x[0])
# Apply the normalization factor to the convolved signal
signal_blurred /= normalization_factor
# Plot the blurred signal
plt.plot(E, signal_blurred, label='Blurred Signal', color='red')
plt.xlabel('Energy (keV)')
plt.ylabel('Blurred solar axion flux (keV⁻¹cm⁻²s⁻¹)')
plt.legend()
plt.show()
plt.savefig("plots/SolarAxionTheoreticalSignalFluxBlurred.pdf")
plt.close()
"""
# conversion probability plot
g_aγ = np.linspace(1e-13, 1e-10, 1000)
plt.plot(g_aγ, conversionProbability()*(g_aγ**2)) # I multiply by g^2 because in the function I had excluded it by dividing it out.
plt.xlabel("$g_{aγ}$ (GeV⁻¹)")
plt.ylabel("Conversion probability")
#plt.title("Conversion Probability vs g_aγ")
plt.show()
plt.savefig("plots/ConversionProbability.pdf")
plt.close()

#Primakoff solar axion flux
E = np.linspace(0, 10, 200)
g_aγ = 1e-10
solarAxionFlux_values = np.array([solarAxionFlux(e) for e in E])
plt.plot(E, solarAxionFlux_values*(g_aγ**2))
plt.xlabel('Energy (keV)')
plt.ylabel('Solar axion flux for $g_{a\gamma}·10^{-10}$ (keV⁻¹ cm⁻² s⁻¹)')
plt.savefig("plots/SolarAxionPrimakoffFlux.pdf")
plt.close()

solarAxionFlux2019_values = np.array([solarAxionFlux2019(e) for e in E])
plt.plot(E, solarAxionFlux2019_values*(g_aγ**2))
plt.xlabel('Energy (keV)')
plt.ylabel('Solar axion flux 2019 for $g_{a\gamma}·10^{-10}$ (keV⁻¹ cm⁻² s⁻¹)')
plt.show()
plt.savefig("plots/SolarAxionPrimakoffFlux2019.pdf")
plt.close()

"""
# My signal, which is the expected solar axion flux corrected by my efficiencies
g_aγ4 = 1e-40
signal_values = [signal(e, g_aγ4) for e in E]
plt.plot(E, signal_values)
plt.xlabel('Energy (keV)')
plt.ylabel('Received solar axion flux (keV⁻¹cm⁻²s⁻¹)')
plt.show()
plt.savefig("plots/SolarAxionTheoreticalSignalFlux.pdf")
plt.close()
"""
"""
# Total signal Ar1
g_aγ4 = np.linspace(1e-52, 1e-40, 1000)
totalSignal_values = [totalSignal(dataset_Ar1, g4) for g4 in g_aγ4]
plt.plot(g_aγ4, totalSignal_values)
plt.xlabel('$g_{aγ}^4$ (GeV$^{-4}$)')
plt.ylabel('Total signal (counts)')
plt.show()
plt.savefig("plots/TotalSignal.pdf")
plt.close()
"""

# Total signal all datasets
g_aγ4 = np.linspace(1e-52, 1e-40, 1000)
datasets = [dataset_Ar1, dataset_Ar2, dataset_Xe]
signal_total = np.zeros_like(g_aγ4)
for dataset in datasets:
    signal_individual = np.array([totalSignal(dataset, g4) for g4 in g_aγ4])
    plt.plot(g_aγ4, signal_individual, label=f"{dataset.name}")
    signal_total += signal_individual
plt.plot(g_aγ4, signal_total, label='Total', color='black', linestyle='--')
plt.xlabel('$g_{aγ}^4$ (GeV$^{-4}$)')
plt.ylabel('Total signal (counts)')
plt.legend()
plt.show()
plt.savefig("plots/TotalSignalCombined.pdf")
plt.close()

# Telescope efficiency
plt.plot(dfTel["E(keV)"], dfTel["Efficiency"])
plt.xlabel("Energy (keV)")
plt.ylabel("Efficiency")
plt.show()
plt.savefig("plots/telescopeEfficiency.pdf")
plt.close()
