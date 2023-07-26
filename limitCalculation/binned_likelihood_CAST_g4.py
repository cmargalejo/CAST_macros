from codecs import backslashreplace_errors
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
from bisect import bisect_left
import scipy.integrate as sc # simpson numerical integration routine

#These are real data as of 14th july 2023
totalTime = 337.0 * 60 * 60 # 337 h of "tracking time"
areaBore = math.pi * (2.15 * 2.15) # cm² because bore diameter is 4.3 cm
chipArea = 0.16*np.pi # in cm. I assume all flux is focused into a circular area of radius 4 mm. Thus pi*r^2= pit* 0.4^2 mm² on the detector. Relevant area for background!

BinWidth = 1.0 # in keV, the bin width of the energy bins
Energies = np.array([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5]) # Energy of the center of each bin in keV
#Background = np.array([3.19e-8, 8.92e-7, 2.74e-6, 2.58e-6, 1.70e-6, 1.74e-6, 1.72e-6, 3.06e-6, 6.07e-6, 3.42e-6]) # for R = 10 mm converted to a rate in keV⁻¹cm⁻²s⁻¹
Background = np.array([0, 1.75e-6, 1.97e-6, 1.75e-6, 1.53e-6, 1.53e-6, 1.31e-6, 2.62e-6, 6.55e-6, 2.62e-6]) # for R = 4 mm
Candidates = np.array([0,      2,      0,     0,      0,      0,       0,      2,    8,      0]) #number of candidates in each energy bin or channel
EnergiesSoftware = np.array([1.75, 2.5,   3.75,    5.25,   7,   12]) #right edge of each energy bin
SoftwareEfficiencyAr1 = np.array([0.65, 0.68,   0.77,     0.71,   0.79,   0.74])
SoftwareEfficiencyAr2 = np.array([0.62, 0.79,   0.75,     0.68,   0.79,   0.71])
SoftwareEfficiencyXe = np.array([0.80, 0.94,   0.84,     0.81,   0.88,   0.87])

def solarAxionFlux(ω: float) -> float: # ω is axion energy in keV, in units 10^-12
    # From CAST constraints on the axion-electron coupling (2013)
    # axion flux produced by the Primakoff effect in solar core in units of keV⁻¹m⁻²yr⁻¹
    flux = (2.0 * (10 ** 18) * (1 / (10 ** -12)) **2 * ω**2.450 * np.exp(-0.829 * ω ))
    # convert flux to correct units --> keV⁻¹m⁻²yr⁻¹ = 3.17098e-12 keV⁻¹cm⁻²s⁻¹ !!!!
    return flux * 3.17098e-12 #in the correct units keV⁻¹cm⁻²s⁻¹

def solarAxionFlux2019(E: float) -> float:
    phi10 = 6.02 * 10 ** 10 # keV⁻¹cm⁻²s⁻¹
    flux2019 = phi10 * (1 / (10 ** -10)) ** 2 * E ** 2.481 / np.exp(E/1.205)
    return flux2019


def conversionProbability():
    # the conversion probability in the CAST magnet ( does not depend on g_aγ esxplicitly because we have taken it out to be able to raise it
    # later to the power of 4 without turning it >0)
    # simplified vacuum conversion prob. for small masses
    B = 9.0 * 195.353 #to convert into natural units. I need the factor 195.353 and the resulting units for B turn into eV^2. See tesla_conversion.nim
    L = 9.26 * 5.06773e+06 # resulting units in eV^-1
    return (1 / 10**9 * B * L / 2.0) ** 2 # divided by 10**9 to convert from GeV^-1 to eV^-1


# Here, dfTel and lerpTel are global variables that are defined independently of the function telescopeEff.
# So they are only done once, and then used inside the function if this is called.
# This csv file assumes parallel light. But in our case whe get light from the solar core, which is about 3 arcmin,
# so the efficiency is lower, like in the paper of 2013.
dfTel = pd.read_csv("data/llnl_xray_telescope_cast_effective_area_parallel_light_DTU_thesis.csv")
dfTel["Efficiency"] = dfTel["EffectiveArea[cm²]"] / areaBore * (8.174 / 10.15) # The numbers scale the parallel light effective area to the CAST JCAP 2015 effective area
#print(dfTel)
lerpTel = interp1d(dfTel["Energy[keV]"], dfTel["Efficiency"], bounds_error = False, fill_value = 0.0)
def telescopeEff(E):
    return lerpTel(E)

# Including the software efficiency by interpolation
#lerpSoftEffAr1 = interp1d(EnergiesSoftware, SoftwareEfficiencyAr1, bounds_error = False, fill_value = 0.0)
#def softwareEff(E):
#    return lerpSoftEffAr1(E)

def softwareEff(E):
    idx = 0
    for e in EnergiesSoftware:
        if E < e:
            break
        idx = idx + 1
    return SoftwareEfficiencyAr1[idx]

  

#Including the detector (gas+mylar) efficiecny
dfDetAr = pd.read_csv("data/ArgonAndWindowEfficiency.csv")
lerpDetAr = interp1d(dfDetAr["Photon energy [keV]"], dfDetAr["Efficiency"], bounds_error = False, fill_value = 0.0)
def detectorArEff(E):
    return lerpDetAr(E)

def totalSignal(g_aγ4):
    ## Flux integrated to total time, energy and area in counts of X-rays.
    # 1. integrate the solar flux
    ## NOTE: in practice this integration must not be done in this proc! Only perform once!
    xs = np.linspace(0.0, 10.0, 100)
    fl = np.array([solarAxionFlux(x) * telescopeEff(x) * detectorArEff(x) * softwareEff(x) for x in xs])
    #fl = np.array([solarAxionFlux(x, g_aγ) for x in xs])
    integral = sc.simpson(fl, # convert units to float for compatibility
                          xs) # convert back to units (integrated out `keV⁻¹`!)
    # 2. compute final flux by "integrating" out the time and area
    #print("Integral of axion flux * efficiency = ", integral, " time ", totalTime, " bore ", areaBore, " prob ", conversionProbability(g_aγ))
    #print("Total signal [counts] = ", integral * totalTime * areaBore * conversionProbability(g_aγ))
    return integral * totalTime * areaBore * conversionProbability() * g_aγ4

## NOTE: only important that signal and background have the same units!
# It gives the amount of signal in counts keV⁻¹ expected for each channel.
# In this case channels are energy bins, so correspond to energy E.
def signal(E, g_aγ4): # in counts keV⁻¹
    ## Returns the axion flux based on `g` and energy `E`
    return solarAxionFlux(E) * totalTime * areaBore * conversionProbability() * telescopeEff(E) * detectorArEff(E) * softwareEff(E) * g_aγ4 #be careful because this area is that of the bore, not the spot on the readout


def background(idx):
    ## Compute an interpolation of energies and background? Not necesarily in the binned approach.
    ## NOTE: For simplicity we only evaluate at the channel energies anyway. In practice
    ## one likely wants interpolation to handle all energies in the allowed range correctly!
    ## idx = Energies.lowerBound(E) # get idx of this energy
    ## Note: area of interest is the region on the chip, in which the signal is focused!
    ## This also allows us to see that the "closer" we cut to the expected axion signal on the
    ## detector, the less background we have compared to the *fixed* signal flux!
    #print("Background = ", (Background[idx] * totalTime * chipArea))
    return (Background[idx] * totalTime * chipArea) # in counts keV⁻¹  #be careful because this area is that of the spot on the readout, not the bore

def totalBackground(energies):
    #This computes the total number of background counts based on my background model.
    integralBackground = 0.0
    for index in range(energies.size):
        b = background(index)
        integralBackground += b * BinWidth
    return integralBackground
print(f"Total number of background counts =  {totalBackground(Energies)}" )

def likelihood(g_aγ4: float, energies: np.ndarray, cs: np.ndarray) -> float:
    """
    `energies` = energies corresponding to each channel
    `cs` = each element is number of counts in that energy channel
    """
    result = np.exp(-totalSignal(g_aγ4)) # `e^{-s_tot}`

    for i in range(cs.size):
        c = cs[i]       # number of candidates in this channel
        E = energies[i] # energy of this channel
        s = signal(E, g_aγ4)
        b = background(i)
        result *= pow(1 + s / b, float(c))

    return result

def logLikelihood(g_aγ4: float, energies: np.ndarray, cs: np.ndarray):
    return -np.log(likelihood(g_aγ4, energies, cs))

def likelihood2(g_aγ4: float, energies: np.ndarray, cs: np.ndarray) -> float:
    result = np.exp(-(totalSignal(g_aγ4)+totalBackground(energies))) #e^(-(s_tot + b_tot))

    for i in range(cs.size):
        c = cs[i]
        E = energies[i]
        s = signal(E, g_aγ4)
        b = background(i)
        result *= pow(s+b, c) / math.factorial(c)
    
    return result

def logLikelihood2(g_aγ4: float, energies: np.ndarray, cs: np.ndarray):
    return -np.log(likelihood2(g_aγ4, energies, cs))       

# Now we are going to compute an expected limit.
# First we need to write a function to compute a limit.
def limit(energies: np.ndarray, cs: np.ndarray, likelihoodFunction=likelihood):
    # Compute limit, CDF@95%
    # limit needs non logspace x & y data! (at least if computed in this simple way)
    g4Lin = np.linspace(0.0, 1e-40, 1000)
    likelihoodLin = [likelihoodFunction(g4, energies, cs) for g4 in g4Lin]
    LCumSum = np.cumsum(likelihoodLin)          # cumulative sum
    LMax = LCumSum.max()               # maximum of the cumulative sum
    LCdf = LCumSum / LMax              # normalize to get (empirical) CDF
    limitIdx = bisect_left(LCdf, 0.95) # limit at 95% of the CDF

    return g4Lin[limitIdx]


def limitViaIntegral(energies: np.ndarray, cs: np.ndarray):
    # Compute limit, CDF@95%
    # limit needs non logspace x & y data! (at least if computed in this simple way)
    num = 1000
    g4Lin = np.linspace(0.0, 1e-40, num)
    likelihoodLin = [likelihood(g4, energies, cs) for g4 in g4Lin]
    fullIntegral = sc.simpson(likelihoodLin, g4Lin)
    # this does not work. It was an attempt to recover the old limit by a change of variables from dg^4 to dg.
        #gLin = np.linspace(0.0, 1e-10, num)
        #likelihoodLin = np.array([likelihood(g, energies, cs) * 4*g**3 for g in gLin])
        #fullIntegral = sc.simpson(likelihoodLin, gLin)
    # Here it would be using g, so I would recover the old limit, the one that integrates over dg instead of dg^4.
    #gLin = np.linspace(0.0, 1e-10, num)
    #likelihoodLin = np.array([likelihood(g**4, energies, cs) for g in gLin])
    #fullIntegral = sc.simpson(likelihoodLin, gLin)
    #int95 = fullIntegral * 0.95
    integrals = np.array([sc.simpson(likelihoodLin[0:i], g4Lin[0:i]) for i in range(1, num)])
    integrals = integrals / fullIntegral
    limitIdx = bisect_left(integrals, 0.95)

    #plt.plot(gLin, likelihoodLin) # to make a line plot
    plt.scatter(g4Lin, likelihoodLin) # to make dots
    plt.xlabel("Coupling constant (GeV⁻4)")
    plt.ylabel("Likelihood")
    #plt.xscale('log')
    #plt.show()
    plt.axvline(x=g4Lin[limitIdx], color='r')
    plt.savefig("energy_bins_likelihood_g4_via_integral.pdf")
    plt.close()

    return g4Lin[limitIdx]

# def limit4(energies: np.ndarray, cs : np.ndarray, toPlot: bool):



def drawCandidates(energies: np.ndarray):
    cs = np.zeros(energies.size) # it produces an array filled with zeroes, i.e. `np.zeros(2) == [0, 0]` (as numpy)
    for i in range(energies.size):
        b = background(i)
        c = np.random.poisson(b) #number of candidates in bin i, based on a Poisson distribution with mean b
        cs[i] = c
    return cs
# Compute the expected limit, i.e., the median of the possible limits we can get

def MonteCarloLim(energies: np.ndarray, nmc): # nmc = number Monte Carlo
    limits = np.zeros(nmc)
    for i in range(nmc):
        cs = drawCandidates(energies)
        lim = limit(energies, cs)
        print("Toy (expected) limit i = ", i, " is = ", pow(lim, 0.25)) 
        limits[i] = lim
    return limits
    
def AsimovDataset():
    ad = np.zeros(Background.size)
    for i in range(Background.size):        
        ad[i] = background(i) * 1 # where 1 is the bin width in keV
    return ad


# Compute the real limit
lim = limit(Energies, Candidates)
print(f"\033[1;35;40m Limit at : {pow(lim, 0.25)}\033[0m")
lim2 = limit(Energies, Candidates, likelihood2)
print(f"\033[1;35;40m Limit via Nature paper approach at : {pow(lim2, 0.25)}\033[0m")
limInt = limitViaIntegral(Energies, Candidates)
#print(f"Limit via integral at : {pow(limInt, 0.25)}")
print(f"\033[1;35;40m Limit via integral at : {pow(limInt, 0.25)}\033[0m")



# print(f"Limit via integral at : {limInt}") # if I use g instead of g^4
# Compute the expected limit

limAsimov = limit(Energies, AsimovDataset())
print(f"\033[1;35;40m Expected limit via Asimov dataset at : {pow(limAsimov,0.25)}\033[0m")

MCLim = MonteCarloLim(Energies, 1) #1000 takes about 2.5 hours
np.savetxt("expected_limits_g4.csv", MCLim, delimiter = ",")
print(f"\033[1;35;40m Expected limit via median at : {pow(np.median(MCLim),0.25)}\033[0m")

# Make some plots

# Histogram of the expected limit array
plt.hist(MCLim)
plt.xlabel("Monte Carlo limits at coupling constant (GeV$^{-4}$)")
plt.ylabel("Counts")
plt.axvline(x=np.median(MCLim), color='r')
plt.savefig("MCLimits_g4.pdf")
plt.close()

# conversion probability plot
g_aγ = np.linspace(1e-13, 1e-10, 1000)
plt.plot(g_aγ, conversionProbability()*(g_aγ**2)) # I multiply by g^2 because in the function I had excluded it by dividing it out.
plt.xlabel("$g_{aγ}$ (GeV⁻¹)")
plt.ylabel("Conversion probability")
#plt.title("Conversion Probability vs g_aγ")
plt.savefig("ConversionProbability.pdf")
plt.close()

#Primakoff solar axion flux
E = np.linspace(0, 10, 200)
g_aγ = 1e-10
solarAxionFlux_values = np.array([solarAxionFlux(e) for e in E])
plt.plot(E, solarAxionFlux_values*(g_aγ**2))
plt.xlabel('Energy (keV)')
plt.ylabel('Solar axion flux for $g_{a\gamma}·10^{-10}$ (keV⁻¹ cm⁻² s⁻¹)')
plt.savefig("SolarAxionPrimakoffFlux.pdf")
plt.close()

solarAxionFlux2019_values = np.array([solarAxionFlux2019(e) for e in E])
plt.plot(E, solarAxionFlux2019_values*(g_aγ**2))
plt.xlabel('Energy (keV)')
plt.ylabel('Solar axion flux 2019 for $g_{a\gamma}·10^{-10}$ (keV⁻¹ cm⁻² s⁻¹)')
plt.savefig("SolarAxionPrimakoffFlux2019.pdf")
plt.close()

# My signal, which is the expected solar axion flux corrected by my efficiencies
g_aγ4 = 1e-40
signal_values = [signal(e, g_aγ4) for e in E]
plt.plot(E, signal_values)
plt.xlabel('Energy (keV)')
plt.ylabel('Received solar axion flux (keV⁻¹cm⁻²s⁻¹)')
plt.savefig("SolarAxionTheoreticalSignalFlux.pdf")
plt.close()

# Total signal
g_aγ4 = np.linspace(1e-52, 1e-40, 1000)
totalSignal_values = [totalSignal(g4) for g4 in g_aγ4]
plt.plot(g_aγ4, totalSignal_values)
plt.xlabel('$g_{aγ}^4$ (GeV$^{-4}$)')
plt.ylabel('Total signal (counts)')
plt.savefig("TotalSignal.pdf")
plt.close()

# Telescope efficiency
plt.plot(dfTel["Energy[keV]"], dfTel["Efficiency"])
plt.xlabel("Energy (keV)")
plt.ylabel("Efficiency")
plt.savefig("telescopeEfficiency.pdf")
plt.close()

#loglikelihood plot, i.e. -chi^2/2 vs gag^4 ??
# This does not work properly because when I raise g to the power of 4 I make it always >0
# I need to re-write the full script but making the functions `solarFlux` and `conversionProbability` 
# independent of g, and then introduce g later.
# see the new script binned_lokelihood_CAST_g4.py
g_aγs = np.linspace(-1e-40, 1e-40, 1000)
logL = [logLikelihood(g_aγ, Energies, Candidates) for g_aγ in g_aγs]
plt.plot(g_aγs, logL)
plt.xlabel(' $g_{aγ}^4$ (GeV$^{-4}$)')
plt.ylabel('-log L or Chi$^2/2$')
plt.savefig("ChiSquareg4.pdf")
plt.close()

g_aγs = np.linspace(-1e-10, 1e-40, 1000)
L = [likelihood(g_aγ, Energies, Candidates) for g_aγ in g_aγs]
plt.plot(g_aγs, L)
plt.xlabel(' $g_{aγ}^4$ (GeV$^{-4}$)')
plt.ylabel('L')
plt.savefig("likelihood_negative_g4.pdf")
plt.close()

g4Lin = np.linspace(0.0, 1e-40, 1000)
likelihoodLin = [likelihood(g4, Energies, Candidates) for g4 in g4Lin]
plt.plot(g4Lin, likelihoodLin)
plt.xlabel("Coupling constant (GeV$^-4$)")
plt.ylabel("Likelihood")
#plt.xscale('log')
#plt.show()
plt.axvline(x=lim, color='r')
plt.savefig("energy_bins_likelihood_g4_bottom.pdf")
plt.close()


likelihoodLin2 = [likelihood2(g4, Energies, Candidates) for g4 in g4Lin]
plt.plot(g4Lin, likelihoodLin2)
plt.xlabel("Coupling constant (GeV$^-4$)")
plt.ylabel("Likelihood2")
#plt.xscale('log')
#plt.show()
plt.axvline(x=lim, color='r')
plt.savefig("energy_bins_likelihood_g4_Nature_approach.pdf")
plt.close()

g_aγs = np.linspace(-1e-40, 1e-40, 1000)
logL2 = [logLikelihood2(g_aγ, Energies, Candidates) for g_aγ in g_aγs]
plt.plot(g_aγs, logL2)
plt.xlabel(' $g_{aγ}^4$ (GeV$^{-4}$)')
plt.ylabel('-log L or Chi$^2/2$')
plt.savefig("ChiSquareg4_Nature_approach.pdf")
plt.close()