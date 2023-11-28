from unbinned_likelihood_CAST_g4_datasets import *
import pickle

# The total time increases 3 h per day, so I can add 30 hours in each iteration to get it every 10 days, and do it 36 times to have a full year (360 days)
# dataset.totalTime

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

def minusLogLikelihoodCris_with_candidates(dataset, g_aγ4, candidates) -> float:
    # Calculate the total signal and background over the entire dataset
    total_signal = totalSignal(dataset, g_aγ4)
    total_background = totalBackground(dataset)
    result = total_signal + total_background  # ((s_tot + b_tot))
    # Iterate over the candidates and calculate the contribution to the likelihood
    for candidate in candidates:  # candidates is a list of tuples (E, x, y)
        E, x, y = candidate
        s = signal(dataset, E, g_aγ4, x, y)
        b = background(dataset, E)
        result -= np.log(s + b)  # Subtract the log of the likelihood contribution of this candidate
    return result


def likelihoodCris_with_candidates(dataset, g_aγ4, candidates) -> float:
    return np.exp(-minusLogLikelihoodCris_with_candidates(dataset, g_aγ4, candidates))


def limit_with_candidates(dataset, likelihoodFunction, candidates):
    # Compute the likelihood line for this set of candidates
    g4Lin = np.linspace(0.0, 5e-40, 5000)
    #likelihoodLin = [likelihoodFunction(dataset, g4, candidates) for g4 in g4Lin]
    likelihoodLin = likelihoodFunction(dataset, g4Lin, candidates)
    LCumSum = np.cumsum(likelihoodLin)
    LMax = LCumSum.max()
    LCdf = LCumSum / LMax
    limitIdx = bisect_left(LCdf, 0.95)
    return g4Lin[limitIdx]  # Store the computed limit for this iteration



def MCLimit(dataset, energies, x_min, x_max, y_min, y_max, nmc, likelihoodFunction):
    limits = np.zeros(nmc)  # Array to store the limit for each Monte Carlo iteration
    for i in range(nmc):
        # Draw candidates for this iteration
        candidates = drawCandidates(dataset, energies, x_min, x_max, y_min, y_max)
        #filename = f"drawn_candidates_{i+1}.csv"
        #np.savetxt(filename, candidates, delimiter=",")
        limits[i] = limit_with_candidates(dataset, likelihoodFunction, candidates)
        print(f"Toy (expected) limit i = {i} is = {pow(limits[i], 0.25)}")

    # Save the limits to a text file after all iterations are complete
    #np.savetxt("unbinned_expected_limits_g4_fast_arrays_Ar1.csv", limits, delimiter=",")
    
    # Return the limits
    return limits


# Example energy bins and number of simulations
energies = np.linspace(0.0, 12.0, 1000) # energy of the candidates
nmc = 1000  # Number of Monte Carlo simulations 1000 takes 330 minutes usign for loops.
results = []

timeStep = 30 * 60 * 60 # in seconds
for dataset in [dataset_Ar1, dataset_Xe]:
    resultDataset = []
    for i in range(36):
        dataset.total_time = dataset.total_time + i*timeStep # this modifies the datasets! As it is not the main code it is fine, but keep it in mind.
                                                            # it is also essential for this case, otherwise it wouldn't work.
        MCLimits = MCLimit(dataset, energies, x_min, x_max, y_min, y_max, nmc, likelihoodCris_with_candidates)
        noCandLim = limit_with_candidates(dataset, likelihoodCris_with_candidates, []) # this computes the limit for 0 candidates
        resultDataset.append([i, dataset.total_time, noCandLim, MCLimits])
    results.append(resultDataset)

#np.save("projected_sensitivity.dat", np.array(results))
with open("projected_sensitivity.dat", 'wb') as file:
    pickle.dump(results, file)