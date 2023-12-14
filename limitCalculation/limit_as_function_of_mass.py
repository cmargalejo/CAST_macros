from unbinned_likelihood_CAST_g4_datasets import *
import matplotlib.pyplot as plt
import numpy as np

def conversionProbabilityWithMassAndG(g, E, ma):
    # the conversion probability in the CAST magnet ( does not depend on g_aγ explicitly because we have taken it out to be able to raise it
    # later to the power of 4 without turning it >0)
    B = 9.0 * 195.353 #to convert into natural units. I need the factor 195.353 and the resulting units for B turn into eV^2. See tesla_conversion.nim
    L = 9.26 * 5.06773e+06 # resulting units in eV^-1
    #print("L = ", L, " and B = ", B)
    # simplified vacuum conversion prob. for small masses
    # return (1 / 10**9 * B * L / 2.0) ** 2 # divided by 10**9 to convert from GeV^-1 to eV^-1
    # complete version
    q = ma**2 / (2*E)
    #print("q = ", q, ", ma = ", ma, ", E =", E)
    qL2 = q * L / 2 
    #print("qL2 = ", qL2)
    #P = (g * B * L / 2)**2 * (np.sin( qL2 ) / qL2)**2 # basti's thesis version
    P = (g * B *  np.sin(qL2) / q )**2
    #print("P = ", P)
    #print(P)
    return P

def conversionProbabilityWithMass(E, ma):
    # the conversion probability in the CAST magnet ( does not depend on g_aγ explicitly because we have taken it out to be able to raise it
    # later to the power of 4 without turning it >0)
    if E == 0:
        return 0 # if the energy is 0, we return a probability of 0 instead of nan.
    B = 9.0 * 195.353 #to convert into natural units. I need the factor 195.353 and the resulting units for B turn into eV^2. See tesla_conversion.nim
    L = 9.26 * 5.06773e+06 # resulting units in eV^-1
    # simplified vacuum conversion prob. for small masses
    # return (1 / 10**9 * B * L / 2.0) ** 2 # divided by 10**9 to convert from GeV^-1 to eV^-1
    # complete version
    q = ma**2 / (2*E * 1000) # to express energy in eV
    qL2 = q * L / 2 
    P = (1 / 10**9 * B * np.sin(qL2) / q )**2 #1e9 to convert GeV^-1 to eV^-1
    return P

def totalSignalWithMass(dataset, g_aγ4, ma):
    ## Flux integrated to total time, energy and area in counts of X-rays.
    # 1. integrate the solar flux
    ## NOTE: in practice this integration must not be done in this proc! Only perform once!
    xs = np.linspace(0.0, 10.0, 100)
    fl = np.array([solarAxionFlux(x) * telescopeEff(x) * detectorEff(dataset, x) * softwareEff(dataset, x) * conversionProbabilityWithMass(x, ma) for x in xs])
    #fl = np.array([solarAxionFlux(x, g_aγ) for x in xs])
    #print("fl = ", fl)
    integral = sc.simpson(fl, # convert units to float for compatibility
                          xs) # convert back to units (integrated out `keV⁻¹`!)
    # 2. compute final flux by "integrating" out the time and area
    #print("Integral = ", integral,  conversionProbabilityWithMass(3.0, ma))
    #print("Integral of axion flux * efficiency = ", integral, " time ", totalTime, " bore ", areaBore, " prob ", conversionProbability(g_aγ))
    #print("Total signal [counts] = ", integral * totalTime * areaBore * conversionProbability(g_aγ))
    #print(g_aγ4, " int ", integral)
    return integral * dataset.total_time * areaBore * g_aγ4

## NOTE: only important that signal and background have the same units!
# It gives the amount of signal in counts keV⁻¹ expected for each channel.
# In this case channels are energy bins, so correspond to energy E.
def signalWithMass(dataset, E, g_aγ4, x, y, ma): # in counts keV⁻¹ cm^-2
    ## Returns the axion flux based on `g` and energy `E`
    #print("Solar axion flux = ", solarAxionFlux(E), ", area of bore = ", areaBore, ", conversion probability = ", conversionProbability(), ", telescope eff = ", telescopeEff(E), ", detector eff = ", detectorEff(dataset, E), ", software eff = ", softwareEff(dataset, E), "gag = ", g_aγ4, "candidate weights = ", candidate_weights(dataset, x, y))
    return solarAxionFlux(E) * dataset.total_time * areaBore * conversionProbabilityWithMass(E, ma) * telescopeEff(E) * detectorEff(dataset, E) * softwareEff(dataset, E) * g_aγ4 * candidate_weights(dataset, x, y) #be careful because this area is that of the bore, not the spot on the readout

def likelihoodWithMass(dataset, g_aγ4, ma) -> float:
    result = np.exp(-(totalSignalWithMass(dataset, g_aγ4, ma)+totalBackground(dataset))) #e^(-(s_tot + b_tot))
    cEnergies = np.array(dataset.candidates["Es"]) # this will have to be modified when we have the position
    candidate_pos_x = np.array(dataset.candidates["xs"])
    candidate_pos_y = np.array(dataset.candidates["ys"])
    idx = 0
    #print("Result = ", result)
    for candidate in range(len(cEnergies)): # column names in df are xs, ys, Es
        E = cEnergies[candidate]
        x = candidate_pos_x[candidate]
        y = candidate_pos_y[candidate]
        s = signalWithMass(dataset, E, g_aγ4, x, y, ma)
        b = background(dataset, E)
        result *= s+b    
        #if result < 0.0:
        #    print("got less 0 ", result, " from ", s, " and ", b, " at ", E, " idx ", idx)
        #    quit()
        idx += 1
    return result

def totalLikelihoodWithMass(dataset1, dataset2, dataset3, g_aγ4, ma, likelihoodFunction=likelihood) -> float:
    likelihood_1 = likelihoodFunction(dataset1, g_aγ4, ma)
    likelihood_2 = likelihoodFunction(dataset2, g_aγ4, ma)
    likelihood_3 = likelihoodFunction(dataset3, g_aγ4, ma)
    
    return likelihood_1 * likelihood_2 * likelihood_3

#ma = np.linspace(1e-9, 0, 1000 )

# Now we are going to write a function to compute a limit.
def limitWithMass(dataset, ma, likelihoodFunction=likelihood):
    # Compute limit, CDF@95%
    # limit needs non logspace x & y data! (at least if computed in this simple way)
    g4Lin = np.linspace(0.0, 5e-40, 1000)
    likelihoodLin = likelihoodFunction(dataset, g4Lin, ma)
    LCumSum = np.cumsum(likelihoodLin)          # cumulative sum
    LMax = LCumSum.max()               # maximum of the cumulative sum
    LCdf = LCumSum / LMax              # normalize to get (empirical) CDF
    limitIdx = bisect_left(LCdf, 0.95) # limit at 95% of the CDF

    return g4Lin[limitIdx]

# And now we write a fucntion to compute the combined limit
def totalLimitWithMass(dataset1, dataset2, dataset3, ma, likelihoodFunction=likelihood):
    # Compute limit, CDF@95%
    # limit needs non logspace x & y data! (at least if computed in this simple way)
    g4Lin = np.linspace(0.0, 1e-40, 2000) # mayvbe I need more steps? Or a smaller g range?
    likelihoodLin = totalLikelihoodWithMass(dataset1, dataset2, dataset3, g4Lin, ma, likelihoodFunction)
    LCumSum = np.cumsum(likelihoodLin)          # cumulative sum
    LMax = LCumSum.max()               # maximum of the cumulative sum
    LCdf = LCumSum / LMax              # normalize to get (empirical) CDF
    limitIdx = bisect_left(LCdf, 0.95) # limit at 95% of the CDF

    return g4Lin[limitIdx]

# Example usage
E = 3e3  # in eV, which is 3 keV, max of Primakoff
mas = np.logspace(-4.0, -1.0, 10000)
g = 1e-11 / 10**9 # to express it in eV^-1 instead of GeV^-1
probabilities = conversionProbabilityWithMassAndG(g, E, mas)



print("P at mass 0.01 = ", conversionProbabilityWithMassAndG(g, E, 0.01))

# Plotting (if desired)
plt.plot(mas, probabilities)
plt.xlabel('Axion mass (eV)')
plt.ylabel('Conversion probability')
plt.yscale("log")
plt.xscale("log")
plt.xlabel("Axion mass (eV)")
plt.ylabel("Conversion probability")
#plt.title('Conversion Probability')
plt.grid()
plt.savefig("conversion_probability.pdf")
plt.close()
# plt.show()

limits = []
for ma in mas:
    limit = totalLimitWithMass(dataset_Ar1, dataset_Ar2, dataset_Xe, ma, likelihoodWithMass)
    #limit = limitWithMass(dataset_Ar1, ma, likelihoodWithMass)
    #print(f"m_a = {ma}, limit = {limit**0.25}")
    limits.append( limit**0.25 )

totalLim = totalLimitWithMass(dataset_Ar1, dataset_Ar2, dataset_Xe, 0.001, likelihoodWithMass)
print(f"\033[1;35;40m Combined limit Cris  at : {pow(totalLim, 0.25)}\033[0m")

plt.plot(mas, limits)
plt.yscale("log")
plt.xscale("log")
plt.xlabel("Axion mass (eV)")
plt.grid()
plt.ylabel("g_ag (GeV$^{-1}$)")
plt.tight_layout()
plt.savefig("limitMassesThingies.pdf") #, tight=True)
plt.close()
"""
# first two quick plots
g4Lin = np.linspace(0, 1e-40, 1000)
likelihoodLin = totalLikelihoodWithMass(dataset_Ar1, dataset_Ar2, dataset_Xe, g4Lin, ma, likelihood)
plt.plot(g4Lin, likelihoodLin)
plt.xlabel("Coupling constant (GeV$^{-4}$)")
plt.ylabel("Likelihood")
plt.show()
plt.axvline(x=totalLim, color='r')
plt.savefig(f"energy_bins_likelihood_g4_unbinned_all_datasets_optimized.pdf")
plt.close()
"""
