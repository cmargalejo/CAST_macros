from codecs import backslashreplace_errors
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d, griddata, RectBivariateSpline, RegularGridInterpolator
from bisect import bisect_left
import scipy.integrate as sc # simpson numerical integration routine


# Step 1: Define a Dataset Class
class Dataset:
    def __init__(self, background, software_efficiency, detector_efficiency, total_time, candidates_csv, energies, energies_software, name):
        self.background = background
        self.software_efficiency = software_efficiency
        self.detector_efficiency = detector_efficiency
        self.total_time = total_time * 60 * 60
        self.candidates = pd.read_csv(candidates_csv)
        self.energies_software = energies_software
        self.energies = energies
        self.lerpBackground = interp1d(energies, background, bounds_error = False, fill_value = "extrapolate") 
        dfDetector = pd.read_csv(detector_efficiency)
        self.lerpDetector = interp1d(dfDetector["Photon energy [keV]"], dfDetector["Efficiency"], bounds_error = False, fill_value = 0.0)
        self.name = name
# Example of initializing a Dataset instance (using existing Background array and CSV file as placeholders):
# dataset_example = Dataset(Background, SoftwareEfficiencyAr1, ... , totalTime, "data/cluster_candidates_tracking.csv")

# Note: We'll need the specific arrays and values for other datasets to initialize them.

#These are real data as of 14th july 2023
#totalTime = 337.0 * 60 * 60 # 337 h of "tracking time"
areaBore = math.pi * (2.15 * 2.15) # cm² because bore diameter is 4.3 cm
chipArea = 0.16*np.pi #0.16*np.pi #36 # in cm^2. I assume all flux is focused into a circular area of radius 4 mm. Thus pi*r^2= pi* 0.4^2 cm² on the detector. Relevant area for background!

Energies = np.array([0.0, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.0]) # Energy of the center of each bin in keV
# we add energy 0 and 10 for two reasons: one, to get the correct total background via integral. And two, to get more reasonable values from the interpolation,
# because sometimes we were getting negative values for the background on one end because we were relying on extrapolation.
#Background = np.array([3.19e-8, 8.92e-7, 2.74e-6, 2.58e-6, 1.70e-6, 1.74e-6, 1.72e-6, 3.06e-6, 6.07e-6, 3.42e-6]) # for R = 10 mm converted to a rate in keV⁻¹cm⁻²s⁻¹
#BackgroundAr1 = np.array([0, 0, 1.75e-6, 1.97e-6, 1.75e-6, 1.53e-6, 1.53e-6, 1.31e-6, 2.62e-6, 6.55e-6, 2.62e-6, 2.62e-6]) # for R = 4 mm wrong?
BackgroundAr1 = np.array([0, 1.7856E-06, 1.116E-06, 1.7856E-06, 1.5624E-06, 1.3392E-06, 1.3392E-06, 2.4552E-06, 5.3568E-06, 1.116E-06, 8.92801E-07, 6.69601E-07]) # for R = 4 mm
#BackgroundAr1 = np.array([3.5712E-08, 6.42817E-07, 2.64269E-06, 1.89274E-06, 1.35706E-06, 1.28563E-06, 1.24992E-06, 2.24986E-06, 4.53543E-06, 1.64275E-06, 7.49953E-07, 6.42817E-07]) # for R = 10 mm
#These are dummy backgrounds!!!
#BackgroundAr2 = np.array([0, 0, 1.75e-6, 1.97e-6, 1.75e-6, 1.53e-6, 1.53e-6, 1.31e-6, 2.62e-6, 6.55e-6, 2.62e-6, 2.62e-6]) # for R = 4 mm wrong?
BackgroundAr2 = np.array([0, 0, 3.29448E-06, 6.58896E-06, 4.94172E-06, 1.64724E-06, 0, 0, 1.64724E-06, 3.29448E-06, 0, 3.29448E-06]) # for R = 4 mm
#BackgroundAr2 = np.array([2.63558E-07	0	3.95337E-06	3.42626E-06	1.05423E-06	2.10847E-06	1.05423E-06	3.1627E-06	3.1627E-06	1.05423E-06	5.27117E-07	1.05423E-06]) # for R = 10 mm
#BackgroundXe = np.array([0, 0, 1.75e-6, 1.00e-6, 1.00e-6, 1.53e-6, 1.53e-6, 1.31e-6, 2.62e-6, 6.55e-6, 2.62e-6, 2.62e-6]) # for R = 4 mm wrong?
BackgroundXe = np.array([0, 1.13234E-06, 1.45586E-06, 9.70575E-07, 1.45586E-06, 1.45586E-06, 1.2941E-06, 2.26468E-06, 4.20583E-06, 1.94115E-06, 1.45586E-06, 6.4705E-07]) # for R = 4 mm
#BackgroundXe = np.array([0, 8.28224E-07, 1.55292E-06, 1.24234E-06, 1.60468E-06, 1.31998E-06, 1.44939E-06, 2.76937E-06, 4.94346E-06, 2.61408E-06, 1.8635E-06, 1.19057E-06, 2.13785E-05]) # for R = 10 mm

#Candidates = np.array([0,      2,      0,     0,      0,      0,       0,      2,    8,      0]) #number of candidates in each energy bin or channel
#dfCandidates = pd.read_csv("data/cluster_candidates_tracking.csv")
EnergiesSoftware = np.array([1.75, 2.5,   3.75,    5.25,   7,   12]) #right edge of each energy bin
SoftwareEfficiencyAr1 = np.array([0.65, 0.68,   0.77,     0.71,   0.79,   0.74])
SoftwareEfficiencyAr2 = np.array([0.62, 0.79,   0.75,     0.68,   0.79,   0.71])
SoftwareEfficiencyXe = np.array([0.80, 0.94,   0.84,     0.81,   0.88,   0.87])

# Initializing the three Dataset instances

# Argon dataset 1
dataset_Ar1 = Dataset(
    background=BackgroundAr1,
    software_efficiency=SoftwareEfficiencyAr1,
    detector_efficiency="data/ArgonAndWindowEfficiency.csv",
    total_time=130.367, #in hours
    #candidates_csv="data/cluster_candidates_tracking.csv",
    candidates_csv="data/bg_df_candidates_Ar1.csv",
    energies=Energies,
    energies_software=EnergiesSoftware,
    name="Ar1"
)

# Argon dataset 2
dataset_Ar2 = Dataset(
    background=BackgroundAr2,
    software_efficiency=SoftwareEfficiencyAr2,
    detector_efficiency="data/ArgonAndWindowEfficiency.csv",
    total_time=25.5833, #in hours
    candidates_csv="data/bg_df_candidates_Ar2.csv",
    energies=Energies,
    energies_software=EnergiesSoftware,
    name="Ar2"
)

# Xe dataset
dataset_Xe = Dataset(
    background=BackgroundXe,
    software_efficiency=SoftwareEfficiencyXe,
    detector_efficiency="data/XenonAndNeonAndWindowEfficiency.csv",
    total_time=159, #in hours
    candidates_csv="data/bg_df_candidates_Xe.csv",
    energies=Energies,
    energies_software=EnergiesSoftware,
    name="Xe"
)


# we make a linear interpolation for the background and energy arrays, so that the bckg is not binned anymore
#lerpBackground = interp1d(Energies, Background, bounds_error = False, fill_value = "extrapolate") #I input an energy and it outputs the interpolated bckg. The class interp1d is deprecated!
#lerpBackground = np.interp(E, Energies, Background) # see the background fucntion below
#lerpPosition = interp2d



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

def softwareEff(dataset, E):
    idx = 0
    for e in dataset.energies_software:
        if E < e:
            break
        idx = idx + 1
    return dataset.software_efficiency[idx]

  

#Including the detector (gas+mylar) efficiecny
#dfDetAr = pd.read_csv("data/ArgonAndWindowEfficiency.csv")
#lerpDetAr = interp1d(dfDetAr["Photon energy [keV]"], dfDetAr["Efficiency"], bounds_error = False, fill_value = 0.0)
def detectorEff(dataset, E):
    return dataset.lerpDetector(E)

def totalSignal(dataset, g_aγ4):
    ## Flux integrated to total time, energy and area in counts of X-rays.
    # 1. integrate the solar flux
    ## NOTE: in practice this integration must not be done in this proc! Only perform once!
    xs = np.linspace(0.0, 10.0, 100)
    fl = np.array([solarAxionFlux(x) * telescopeEff(x) * detectorEff(dataset, x) * softwareEff(dataset, x) for x in xs])
    #fl = np.array([solarAxionFlux(x, g_aγ) for x in xs])
    integral = sc.simpson(fl, # convert units to float for compatibility
                          xs) # convert back to units (integrated out `keV⁻¹`!)
    # 2. compute final flux by "integrating" out the time and area
    #print("Integral of axion flux * efficiency = ", integral, " time ", totalTime, " bore ", areaBore, " prob ", conversionProbability(g_aγ))
    #print("Total signal [counts] = ", integral * totalTime * areaBore * conversionProbability(g_aγ))
    return integral * dataset.total_time * areaBore * conversionProbability() * g_aγ4

## NOTE: only important that signal and background have the same units!
# It gives the amount of signal in counts keV⁻¹ expected for each channel.
# In this case channels are energy bins, so correspond to energy E.
def signal(dataset, E, g_aγ4): # in counts keV⁻¹
    ## Returns the axion flux based on `g` and energy `E`
    return solarAxionFlux(E) * dataset.total_time * areaBore * conversionProbability() * telescopeEff(E) * detectorEff(dataset, E) * softwareEff(dataset, E) * g_aγ4 #be careful because this area is that of the bore, not the spot on the readout


def background(dataset, E):
    ## For the unbinned approach, I need to compute an interpolation of energies and background.    
    ## Note: area of interest is the region on the chip, in which the signal is focused!
    ## This also allows us to see that the "closer" we cut to the expected axion signal on the
    ## detector, the less background we have compared to the *fixed* signal flux!
    #print("Background = ", (lerpBackground(E) * totalTime * chipArea))
    return (dataset.lerpBackground(E) * dataset.total_time * chipArea) # in counts keV⁻¹  #be careful because this area is that of the spot on the readout, not the bore
    #return (np.interp(E, Energies, Background) * totalTime * chipArea)

def totalBackground(dataset):
    energies = np.linspace(0.0,10.0,1000)
    backgrounds = background(dataset, energies) # gives an array of 1000 backg. rates
    return sc.simpson(backgrounds, energies)




def likelihood(dataset, g_aγ4) -> float:
    result = np.exp(-(totalSignal(dataset, g_aγ4)+totalBackground(dataset))) #e^(-(s_tot + b_tot))
    cEnergies = dataset.candidates["Es"] # this will have to be modified when we have the position
    idx = 0
    for E in cEnergies: # column names in df are xs, ys, Es
        s = signal(dataset, E, g_aγ4)
        b = background(dataset, E)
        result *= s+b    
        #if result < 0.0:
        #    print("got less 0 ", result, " from ", s, " and ", b, " at ", E, " idx ", idx)
        #    quit()
        idx += 1
    return result

def likelihood2(dataset, g_aγ4) -> float: # Basti's version based on the division of Poissonian probabilities
    result = np.exp(-(totalSignal(dataset, g_aγ4))) #e^(-(s_tot ))
    cEnergies = dataset.candidates["Es"] # this will have to be modified when we have the position
    for E in cEnergies: # column names in df are xs, ys, Es
        s = signal(dataset, E, g_aγ4)
        b = background(dataset, E)
        result *= (1.0 + s/b)
        #print("result ", result, " for candidate energy ", E)
    
    return result

def logLikelihood(dataset, g_aγ4: float):
    return -np.log(likelihood(dataset, g_aγ4))  
     

# Now we are going to write a function to compute a limit.
def limit(dataset, likelihoodFunction=likelihood):
    # Compute limit, CDF@95%
    # limit needs non logspace x & y data! (at least if computed in this simple way)
    g4Lin = np.linspace(0.0, 1e-40, 1000)
    likelihoodLin = [likelihoodFunction(dataset, g4) for g4 in g4Lin]
    LCumSum = np.cumsum(likelihoodLin)          # cumulative sum
    LMax = LCumSum.max()               # maximum of the cumulative sum
    LCdf = LCumSum / LMax              # normalize to get (empirical) CDF
    limitIdx = bisect_left(LCdf, 0.95) # limit at 95% of the CDF

    return g4Lin[limitIdx]

#Now we finally compute the limit
for dataset in [dataset_Ar1, dataset_Ar2, dataset_Xe]:
    print(f"Total number of expected background counts via integral of {dataset.name} =  {totalBackground(dataset)}" )
    lim = limit(dataset, likelihood2)
    print(f"\033[1;35;40m Limit at : {pow(lim, 0.25)}\033[0m")

    #g_aγs = np.linspace(-1e-40, 1e-40, 1000)
    #logL = [logLikelihood(g_aγ, dfCandidates) for g_aγ in g_aγs]
    #plt.plot(g_aγs, logL)
    #plt.xlabel(' $g_{aγ}^4$ (GeV$^{-4}$)')
    #plt.ylabel('-log L or Chi$^2/2$')
    #plt.savefig("ChiSquareg4_Nature_approach_unbinned.pdf")
    #plt.close()

    g4Lin = np.linspace(0.0, 1e-38, 1000)
    likelihoodLin = [likelihood2(dataset, g4) for g4 in g4Lin]
    plt.plot(g4Lin, likelihoodLin)
    plt.xlabel("Coupling constant (GeV$^-4$)")
    plt.ylabel("Likelihood")
    #plt.xscale('log')
    #plt.show()
    #plt.axvline(x=lim, color='r')
    plt.savefig(f"energy_bins_likelihood_g4_unbinned_{dataset.name}.pdf")
    plt.close()