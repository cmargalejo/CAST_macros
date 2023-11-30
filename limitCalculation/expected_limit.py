from unbinned_likelihood_CAST_g4_datasets import *

#MC simulation to compute the expected limit.
x_min = -10 # mm
x_max = 10 # mm
y_min = -10 # mm
y_max = 10 # mm
area=2*2 #cm2 This area should correpond to what we get from x and y min max
#binwidth = 1 # keV
def drawCandidates(dataset, energies: np.ndarray, x_min, x_max, y_min, y_max):
    candidates = []
    binwidth = energies[1] - energies[0] #binwidth between first and second bin
    for energy in energies:
        b = background(dataset, energy) * area * binwidth # so that the units are counts keV⁻¹
        num_candidates = np.random.poisson(b)  # Number of candidates in this energy bin
        for _ in range(num_candidates):
            x = np.random.uniform(x_min, x_max)
            y = np.random.uniform(y_min, y_max)
            E = np.random.uniform(energy-binwidth/2.0, energy+binwidth/2.0) #if fixed binwidth: E+binwidth/2.0 and E-binwidth/2.0 \\ low_bin_edge, high_bin_edge
            candidates.append((E, x, y))
    #print("Generated candidates:", candidates)
    print("Number of generated candidates:", len(candidates))
    return candidates

def likelihood_with_candidates(dataset, g_aγ4, candidates):
    # Calculate the total signal and background over the entire dataset
    total_signal = totalSignal(dataset, g_aγ4)
    total_background = totalBackground(dataset)
    result = np.exp(-(total_signal + total_background))  # e^(-(s_tot + b_tot))    
    # Iterate over the candidates and calculate the contribution to the likelihood
    for candidate in candidates:  # candidates is a list of tuples (E, x, y)
        E, x, y = candidate
        s = signal(dataset, E, g_aγ4, x, y)
        #b = background(dataset, E, x, y)  # Assuming background also needs x, y
        b = background(dataset, E)  
        result *= s + b  # Multiply by the likelihood contribution of this candidate    
    return result

def MCLimit(dataset, energies, x_min, x_max, y_min, y_max, nmc, likelihoodFunction, dataset_name ):
    limits = np.zeros(nmc)  # Array to store the limit for each Monte Carlo iteration
    for i in range(nmc):
        # Draw candidates for this iteration
        candidates = drawCandidates(dataset, energies, x_min, x_max, y_min, y_max)
        #filename = f"drawn_candidates_{i+1}.csv"
        #np.savetxt(filename, candidates, delimiter=",")
        # Compute the likelihood line for this set of candidates
        g4Lin = np.linspace(0.0, 5e-40, 5000)
        #likelihoodLin = [likelihoodFunction(dataset, g4, candidates) for g4 in g4Lin]
        likelihoodLin = likelihoodFunction(dataset, g4Lin, candidates)
        LCumSum = np.cumsum(likelihoodLin)
        LMax = LCumSum.max()
        LCdf = LCumSum / LMax
        limitIdx = bisect_left(LCdf, 0.95)
        limits[i] = g4Lin[limitIdx]  # Store the computed limit for this iteration
        print(f"Toy (expected) limit i = {i} is = {pow(limits[i], 0.25)}")

    # Save the limits to a text file after all iterations are complete
    filename = f"unbinned_expected_limits_g4_fast_arrays_{dataset_name}.csv"
    np.savetxt(filename, limits, delimiter=",")
    #np.savetxt("unbinned_expected_limits_g4_fast_arrays_Ar1.csv", limits, delimiter=",")
    
    # Return the limits
    return limits


def MCLimit_combined(datasets, energies, x_min, x_max, y_min, y_max, nmc, likelihoodFunction):
    limits = np.zeros(nmc)  # Store the limit for each Monte Carlo iteration
    g4Lin = np.linspace(0.0, 5e-40, 5000)  # g4 values
    for i in range(nmc):
        combined_likelihoodLin = np.ones_like(g4Lin)
        for dataset in datasets:
            candidates = drawCandidates(dataset, energies, x_min, x_max, y_min, y_max)
            likelihoodLin = likelihoodFunction(dataset, g4Lin, candidates)
            combined_likelihoodLin *= likelihoodLin  # Combine likelihoods
        LCumSum = np.cumsum(combined_likelihoodLin)
        LMax = LCumSum.max()
        LCdf = LCumSum / LMax
        limitIdx = bisect_left(LCdf, 0.95)
        limits[i] = g4Lin[limitIdx]
        print(f"Toy (expected) limit i = {i} is = {pow(limits[i], 0.25)}")
    np.savetxt("combined_expected_limits_g4_fast_arrays.csv", limits, delimiter=",")
    return limits

# Example energy bins and number of simulations
energies = np.linspace(0.0, 12.0, 1000) #when I'm sure it all works fine.
#energies = np.linspace()
nmc = 50000  # Number of Monte Carlo simulations 1000 takes 330 minutes usign for loops.

# Perform the Monte Carlo limit calculations
datasets = [dataset_Ar1, dataset_Ar2, dataset_Xe]
MCLimits_combined = MCLimit_combined(datasets, energies, x_min, x_max, y_min, y_max, nmc, likelihood_with_candidates)
print(f"Combined expected limit via median at : {pow(np.median(MCLimits_combined), 0.25)}")




# Perform the Monte Carlo limit calculations
#MCLimits = MCLimit(dataset_Ar1, energies, x_min, x_max, y_min, y_max, nmc, likelihood_with_candidates)

datasets = [(dataset_Ar1, "Ar1"), (dataset_Ar2, "Ar2"), (dataset_Xe, "Xe")]
MCLimits = []

for dataset, dataset_name in datasets:
    limits = MCLimit(dataset, energies, x_min, x_max, y_min, y_max, nmc, likelihood_with_candidates, dataset_name)
    MCLimits.append(limits)

# Print out the median expected limit
#print(f"Expected limit via median at : {pow(np.median(MCLimits), 0.25)}")