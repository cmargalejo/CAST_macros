from codecs import backslashreplace_errors
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
from bisect import bisect_left
import scipy.integrate as sc # simpson numerical integration routine

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
    B = 9.0 * 195.353 #to convert into natural units. I need the factor 195.353 and the resulting units for B turn into eV^2. See tesla_conversion.nim
    L = 9.26 * 5.06773e+06 # resulting units in eV^-1
    return (1 / 10**9 * B * L / 2.0) ** 2 # divided by 10**9 to convert from GeV^-1 to eV^-1


dfTel = pd.read_csv("data/llnl_xray_telescope_cast_effective_area_parallel_light_DTU_thesis.csv")
dfTel["Efficiency"] = dfTel["EffectiveArea[cm²]"] / areaBore * (8.174 / 10.15) # The numbers scale the parallel light effective area to the CAST JCAP 2015 effective area
#print(dfTel)
lerpTel = interp1d(dfTel["Energy[keV]"], dfTel["Efficiency"], bounds_error = False, fill_value = 0.0)
def telescopeEff(E):
    return lerpTel(E)

def softwareEff(E):
    idx = 0
    for e in EnergiesSoftware:
        if E < e:
            break
        idx = idx + 1
    return SoftwareEfficiencyAr1[idx]

#Including the detector (gas+mylar) efficiency
dfDetAr = pd.read_csv("data/ArgonAndWindowEfficiency.csv")
lerpDetAr = interp1d(dfDetAr["Photon energy [keV]"], dfDetAr["Efficiency"], bounds_error = False, fill_value = 0.0)
def detectorArEff(E):
    return lerpDetAr(E)

def totalSignal(g_aγ4):
    ## Flux integrated to total time, energy and area in counts of X-rays.
   
    xs = np.linspace(0.0, 10.0, 100)
    fl = np.array([solarAxionFlux(x) * telescopeEff(x) * detectorArEff(x) * softwareEff(x) for x in xs])
    integral = sc.simpson(fl, # convert units to float for compatibility
                          xs) # convert back to units (integrated out `keV⁻¹`!)
    return integral * totalTime * areaBore * conversionProbability() * g_aγ4

def signal(E, g_aγ4): # in counts keV⁻¹
    return solarAxionFlux(E) * totalTime * areaBore * conversionProbability() * telescopeEff(E) * detectorArEff(E) * softwareEff(E) * g_aγ4 #be careful because this area is that of the bore, not the spot on the readout



# conversion probability plot
g_aγ = np.linspace(1e-13, 1e-10, 1000)
plt.plot(g_aγ, conversionProbability()*(g_aγ**2)) # I multiply by g^2 because in the function I had excluded it by dividing it out.
plt.xlabel("$g_{aγ}$ (GeV⁻¹)")
plt.ylabel("Conversion probability")
#plt.title("Conversion Probability vs g_aγ")
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
plt.savefig("plots/SolarAxionPrimakoffFlux2019.pdf")
plt.close()

# My signal, which is the expected solar axion flux corrected by my efficiencies
g_aγ4 = 1e-40
signal_values = [signal(e, g_aγ4) for e in E]
plt.plot(E, signal_values)
plt.xlabel('Energy (keV)')
plt.ylabel('Received solar axion flux (keV⁻¹cm⁻²s⁻¹)')
plt.savefig("plots/SolarAxionTheoreticalSignalFlux.pdf")
plt.close()

# Total signal
g_aγ4 = np.linspace(1e-52, 1e-40, 1000)
totalSignal_values = [totalSignal(g4) for g4 in g_aγ4]
plt.plot(g_aγ4, totalSignal_values)
plt.xlabel('$g_{aγ}^4$ (GeV$^{-4}$)')
plt.ylabel('Total signal (counts)')
plt.savefig("plots/TotalSignal.pdf")
plt.close()

# Telescope efficiency
plt.plot(dfTel["Energy[keV]"], dfTel["Efficiency"])
plt.xlabel("Energy (keV)")
plt.ylabel("Efficiency")
plt.savefig("plots/telescopeEfficiency.pdf")
plt.close()
