from unbinned_likelihood_CAST_g4_datasets import *
import pickle
import typer


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

def likelihood_with_candidates(dataset, g_aγ4, candidates): # it overflows!!!
    # Calculate the total signal and background over the entire dataset
    total_signal = totalSignal(dataset, g_aγ4)
    total_background = totalBackground(dataset)
    result = np.exp(-(total_signal + total_background))  # e^(-(s_tot + b_tot))    
    print("min exp -(stot+btot) = ", np.min(result))
    print("max exp -(stot+btot) = ", np.max(result))
    # Iterate over the candidates and calculate the contribution to the likelihood
    for candidate in candidates:  # candidates is a list of tuples (E, x, y)
        E, x, y = candidate
        s = signal(dataset, E, g_aγ4, x, y)
        #b = background(dataset, E, x, y)  # Assuming background also needs x, y
        b = background(dataset, E)  
        result *= s + b  # Multiply by the likelihood contribution of this candidate    
    return result

def minusLogLikelihoodCris_with_candidates(dataset, g_aγ4, candidates) -> float: # it also overflows!!!
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

def likelihood2_with_candidates(dataset, g_aγ4, candidates) -> float:
    # Calculate the total signal over the entire dataset
    total_signal = totalSignal(dataset, g_aγ4)
    result = np.exp(-total_signal)  # e^(-(s_tot))
    #print("exp -stot = ", result)
    #print("min exp -stot = ", np.min(result))
    #print("max exp -stot = ", np.max(result))
    # Iterate over the candidates and calculate the contribution to the likelihood
    for candidate in candidates:  # candidates is a list of tuples (E, x, y)
        E, x, y = candidate
        s = signal(dataset, E, g_aγ4, x, y)
        b = background(dataset, E)
        result *= (1.0 + s / b)  # Multiply by the likelihood contribution of this candidate
        # Debugging print statement (commented out)
        # print("result ", result, " for candidate energy ", E, "background = ", b, " and signal = ", s)
    return result

def logLikelihood2_with_candidates(dataset, g_aγ4, candidates) -> float:
    # Calculate the total signal over the entire dataset
    total_signal = totalSignal(dataset, g_aγ4)
    result = -total_signal  # -(s_tot))
    #print("exp -stot = ", result)
    print("min exp -stot = ", np.min(result))
    print("max exp -stot = ", np.max(result))
    # Iterate over the candidates and calculate the contribution to the likelihood
    for candidate in candidates:  # candidates is a list of tuples (E, x, y)
        E, x, y = candidate
        s = signal(dataset, E, g_aγ4, x, y)
        b = background(dataset, E)
        result += np.log(1.0 + s / b)  # Multiply by the likelihood contribution of this candidate
        # Debugging print statement (commented out)
        # print("result ", result, " for candidate energy ", E, "background = ", b, " and signal = ", s)
    print("min exp result = ", np.min(result))
    print("max exp result = ", np.max(result))
    return np.exp(result)

def limit_with_candidates(dataset, likelihoodFunction, candidates):
    # Compute the likelihood line for this set of candidates
    g4Lin = np.linspace(0.0, 5e-40, 5000)
    #likelihoodLin = [likelihoodFunction(dataset, g4, candidates) for g4 in g4Lin]
    likelihoodLin = likelihoodFunction(dataset, g4Lin, candidates)
    LCumSum = np.cumsum(likelihoodLin)
    LMax = LCumSum.max()
    LCdf = LCumSum / LMax
    limitIdx = bisect_left(LCdf, 0.95)
    #print("max likelihood = ", np.max(likelihoodLin))
    #print("min likelihood = ", np.min(likelihoodLin))
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
        if limits[i] == 0:
            print(candidates)
            quit()

    # Save the limits to a text file after all iterations are complete
    #np.savetxt("unbinned_expected_limits_g4_fast_arrays_Ar1.csv", limits, delimiter=",")
    
    # Return the limits
    return limits

np.random.seed(51090)
def projection(dataset, output): # output is the file.dat
    # Example energy bins and number of simulations
    energies = np.linspace(0.0, 12.0, 1000) # energy of the candidates
    nmc = 1000  # Number of Monte Carlo simulations 1000 takes 330 minutes usign for loops.
    results = []
    # The total time increases 3 h per day, so I can add 30 hours in each iteration to get it every 10 days, and do it 36 times to have a full year (360 days)
    # Finally I launch it 80 times to approach 30 months, like in the past plot. With Basti's version it should be fine and not overflow.
    # dataset.totalTime
    timeStep = 30 * 60 * 60 # 30 horus in seconds

    resultDataset = []
    dataset.total_time = 0 # to start the simulation including the data taking period
    for i in range(92): 
        dataset.total_time = dataset.total_time + timeStep # this modifies the datasets! As it is not the main code it is fine, but keep it in mind.
                                                            # it is also essential for this case, otherwise it wouldn't work.
        MCLimits = MCLimit(dataset, energies, x_min, x_max, y_min, y_max, nmc, likelihood2_with_candidates)
        noCandLim = limit_with_candidates(dataset, likelihood2_with_candidates, []) # this computes the limit for 0 candidates
        resultDataset.append([i, dataset.total_time, noCandLim, MCLimits])
        #print("Current time = ", dataset.total_time / 3600)
    results.append(resultDataset)

    #np.save("projected_sensitivity.dat", np.array(results))
    with open(output, 'wb') as file:
        pickle.dump(results, file)

def main(argon: bool = False, xenon: bool = False, output: str = ""): 
  if argon:
    projection(dataset_Ar1, output)
  if xenon:
    projection(dataset_Xe, output)
if __name__ == "__main__": 
  typer.run(main)

# to run it with one or other arguments:
# python script.py --argon --output "/path/to/argon.dat"
# python script.py --xenon --output "/path/to/xenon.dat"