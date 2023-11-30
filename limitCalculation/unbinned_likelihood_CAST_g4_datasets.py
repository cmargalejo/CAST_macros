from codecs import backslashreplace_errors
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d, griddata, RectBivariateSpline, RegularGridInterpolator
from bisect import bisect_left
import scipy.integrate as sc # simpson numerical integration routine
from scipy.ndimage import gaussian_filter


#These are real data as of 14th july 2023
#totalTime = 337.0 * 60 * 60 # 337 h of "tracking time"
areaBore = math.pi * (2.15 * 2.15) # cm² because bore diameter is 4.3 cm
x_min = -10.0 # in mm
x_max = 10.0
y_min = -10.0
y_max = 10.0
readoutArea = 1**2 * np.pi #0.16*np.pi #36 # in cm^2. 
# If I assume all flux is focused into a circular area of radius 4 mm. Thus pi*r^2= pi* 0.4^2 cm² on the detector. Relevant area for background!
# Is I assume 10 mm radius, then area = pi * 1^2 cm^2  = pi cm^2
#readoutArea = (x_max - x_min) * (y_max - y_min) / 100.0 # in cm
#print("area= ", readoutArea)

#Energies = np.array([0.0, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.0]) # Energy of the center of each bin in keV --> this should go from 0.5 t0 11.5. Think about it!
Energies = np.array([ 0.0, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.0 ]) # If I use this I need to include 0s in BackroungAr1, etc.
# we add energy 0 and 10 for two reasons: one, to get the correct total background via integral. And two, to get more reasonable values from the interpolation,
# because sometimes we were getting negative values for the background on one end because we were relying on extrapolation.
#BackgroundAr1 = np.array([0, 1.7856E-06, 1.116E-06, 1.7856E-06, 1.5624E-06, 1.3392E-06, 1.3392E-06, 2.4552E-06, 5.3568E-06, 1.116E-06, 8.92801E-07, 6.69601E-07]) # for R = 4 mm
BackgroundAr1 = np.array([0, 3.5712E-08, 6.42817E-07, 2.64269E-06, 1.89274E-06, 1.35706E-06, 1.28563E-06, 1.24992E-06, 2.24986E-06, 4.53543E-06, 1.64275E-06, 7.49953E-07, 6.42817E-07, 0]) # for R = 10 mm
#BackgroundAr2 = np.array([0, 0, 3.29448E-06, 6.58896E-06, 4.94172E-06, 1.64724E-06, 0, 0, 1.64724E-06, 3.29448E-06, 0, 3.29448E-06]) # for R = 4 mm
BackgroundAr2 = np.array([0, 2.63558E-07, 0, 3.95337E-06, 3.42626E-06, 1.05423E-06, 2.10847E-06, 1.05423E-06, 3.1627E-06, 3.1627E-06, 1.05423E-06, 5.27117E-07, 1.05423E-06, 0]) # for R = 10 mm
#BackgroundXe = np.array([0, 1.13234E-06, 1.45586E-06, 9.70575E-07, 1.45586E-06, 1.45586E-06, 1.2941E-06, 2.26468E-06, 4.20583E-06, 1.94115E-06, 1.45586E-06, 6.4705E-07]) # for R = 4 mm
BackgroundXe = np.array([0, 1e-8, 8.28224E-07, 1.55292E-06, 1.24234E-06, 1.60468E-06, 1.31998E-06, 1.44939E-06, 2.76937E-06, 4.94346E-06, 2.61408E-06, 1.8635E-06, 1.19057E-06, 0]) # for R = 10 mm

#Candidates = np.array([0,      2,      0,     0,      0,      0,       0,      2,    8,      0]) #number of candidates in each energy bin or channel
#dfCandidates = pd.read_csv("data/cluster_candidates_tracking.csv")
EnergiesSoftware = np.array([1.75, 2.5,   3.75,    5.25,   7,   12]) #right edge of each energy bin
SoftwareEfficiencyAr1 = np.array([0.65, 0.68,   0.77,     0.71,   0.79,   0.74])
SoftwareEfficiencyAr2 = np.array([0.62, 0.79,   0.75,     0.68,   0.79,   0.71])
SoftwareEfficiencyXe = np.array([0.80, 0.94,   0.84,     0.81,   0.88,   0.87])

#axion_image_filename = '/home/cristina/GitHub/CAST_macros/limitCalculation/data/llnl_raytracing_Jaime_all_energies.txt'
#axion_image_filename = '/home/cris/src/GitHub/CAST_macros/limitCalculation/data/llnl_raytracing_Jaime_all_energies.txt'
axion_image_filename = 'data/llnl_raytracing_Jaime_all_energies_raw_sum.txt'
convolved_resolution = 500

def setupAxionImageInterpolator(filename, convolved_resolution): #I should check individually that this sums up to 1, by computing the
    # Read the text file using Pandas
    df = pd.read_csv(filename, skiprows = 1, delim_whitespace = True, names = ["x", "y", "z", "zMean"])
    # Compute normalized data. Using *MEAN* of *DATA*
    #df["x"] = df["x"] - df["x"].mean() + 1.3 #because the mean is 31.3, but I have to move it only 30 mm
    #df["y"] = df["y"] - df["y"].mean() + 0.025 #because the mean is 30.25, but I have to move it only 30 mm
    df["x"] = df["x"] - 31.114 # these are the values of the centroid that match it with the position of the X-ray finger cetnroid
    df["y"] = df["y"] - 30.283 # these are the values of the centroid that match it with the position of the X-ray finger cetnroid 30.283
    # Sort data by *X* then *Y*
    df = df.sort_values(by = ["x", "y"], ascending = True)
    # Number of elements per axis
    size = int(np.sqrt(len(df)))
    # get unique elements
    unique_x = pd.unique(df["x"])
    unique_y = pd.unique(df["y"])
    #print("unique x = ", unique_x, "unique y = ", unique_y)
    # Construct array for interpolation
    zs = np.zeros([size, size])
    zSum = np.sum(df["z"])
    areaPerPixel = 0.01 * 0.01 # strip pitch in cm, final units in cm^2. But in the .txt file I have data every 0.1 mm, which is 0.01 cm
    for x in range(size):
        for y in range(size):
            ## Important: First entries are `x`, second `y`. "Matrix" convention, i.e.
            ## first index is *row*, not column!
            z = df["z"][x * size + y] 
            zs[x, y] = (z / zSum / areaPerPixel) # Sorted by `x` first, hence `x` changes NORMALIZATION TO 1 if we cover the entire readout
                                             # *after* y.
    # Convolve the data with a Gaussian kernel
    sigma = convolved_resolution / 2.35482 #by definition, convolved_resolution = FWHM = 2 * sqrt(2 * ln(2)) * sigma = 2.35482 * sigma
    sigma = sigma / 100 # to convert into um
    print("Sigma = ", sigma)
    convolved_zs = gaussian_filter(zs, sigma=sigma, mode='wrap')
    # Create a RegularGridInterpolator
    return RegularGridInterpolator((unique_x, unique_y), convolved_zs, method='linear', bounds_error = False, fill_value = 0.0) #linear, nearest

def candidate_position_transformation(positions_data, x_min, x_max, y_min, y_max, radius):
    # Create specular image by negating x-coordinates
    df = positions_data
    df["xs"] = -df["xs"]
    # Rotate the positions by 0 degrees clockwise (up to the right)
    angle_degrees = 45
    xs = df["xs"].copy()
    ys = df["ys"].copy()
    df["xs"] = xs * np.cos(np.radians(angle_degrees))  + ys * np.sin(np.radians(angle_degrees))
    df["ys"] = xs * -np.sin(np.radians(angle_degrees)) + ys * np.cos(np.radians(angle_degrees))
    # Move positions 1.5 mm to the right to align it with the X-ray finger simulations
    # If I use the centroid values above in setupAxionImageInterpolator(), the shift should be 0.0 mm
    shift_distance = 0.0
    df["xs"] = df["xs"] + shift_distance
    # Filter out CSV values that are within the wanted area. For a square:
    #df = df[(df["xs"] > x_min) & (df["xs"] < x_max) & (df["ys"] > y_min) & (df["ys"] < y_max)] # rectangle
    # Filter out CSV values that are within the wanted area. For a circle:
    df = df[(df["xs"] * df["xs"] + df["ys"] * df["ys"] < radius*radius)] #100 ir radius^2
    # Filter by energy
    df = df[(df["Es"] < 10)] # Only events below 10 keV
    return df

# Step 1: Define a Dataset Class
class Dataset:
    def __init__(self, background, software_efficiency, detector_efficiency, total_time, candidates_csv, energies, energies_software, name, axion_image):
        self.background = background
        self.software_efficiency = software_efficiency
        self.detector_efficiency = detector_efficiency
        self.total_time = total_time * 60 * 60
        df = pd.read_csv(candidates_csv)
        print("Number of candidates total readout area: ", len(df))
        self.candidates = candidate_position_transformation(df, -10.0, 10.0, -10.0, 10.0, radius=10.0) # these values define the area
        print("Number of candidates in the inner 10 mm radius circle: ", len(self.candidates))
        self.energies_software = energies_software
        self.energies = energies
        self.lerpBackground = interp1d(energies, background, bounds_error = False, fill_value = "extrapolate") 
        dfDetector = pd.read_csv(detector_efficiency)
        self.lerpDetector = interp1d(dfDetector["Photon energy [keV]"], dfDetector["Efficiency"], bounds_error = False, fill_value = 0.0)
        self.name = name
        self.axion_image_weights = setupAxionImageInterpolator(axion_image, 500)

# Initializing the three Dataset instances

# Argon dataset 1
dataset_Ar1 = Dataset(
    background=BackgroundAr1,
    software_efficiency=SoftwareEfficiencyAr1,
    detector_efficiency="data/ArgonAndWindowEfficiency.csv",
    total_time=130.367, #in hours
    #candidates_csv="data/cluster_candidates_tracking.csv",
    candidates_csv="data/bg_df_candidates_Ar1.csv",
    #candidates_csv="data/bg_df_candidates_many_candidates_Ar1.csv", #data/bg_df_candidates_zero_candidates.csv data/bg_df_candidates_many_candidates_Ar1.csv 
    energies=Energies,
    energies_software=EnergiesSoftware,
    name="Ar1",
    axion_image = axion_image_filename
)

# Argon dataset 2
dataset_Ar2 = Dataset(
    background=BackgroundAr2,
    software_efficiency=SoftwareEfficiencyAr2,
    detector_efficiency="data/ArgonAndWindowEfficiency.csv",
    total_time=25.5833, #in hours
    candidates_csv="data/bg_df_candidates_Ar2.csv",
    #candidates_csv="data/bg_df_candidates_many_candidates_Ar2.csv", #data/bg_df_candidates_zero_candidates.csv data/bg_df_candidates_many_candidates_Ar2.csv 
    energies=Energies,
    energies_software=EnergiesSoftware,
    name="Ar2",
    axion_image = axion_image_filename
)

# Xe dataset
dataset_Xe = Dataset(
    background=BackgroundXe,
    software_efficiency=SoftwareEfficiencyXe,
    detector_efficiency="data/XenonAndNeonAndWindowEfficiency.csv",
    total_time=159, #in hours
    candidates_csv="data/bg_df_candidates_Xe.csv",
    #candidates_csv="data/bg_df_candidates_many_candidates_Xe.csv", #data/bg_df_candidates_zero_candidates.csv data/bg_df_candidates_many_candidates_Xe.csv 
    energies=Energies,
    energies_software=EnergiesSoftware,
    name="Xe",
    axion_image = axion_image_filename
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
    # the conversion probability in the CAST magnet ( does not depend on g_aγ explicitly because we have taken it out to be able to raise it
    # later to the power of 4 without turning it >0)
    # simplified vacuum conversion prob. for small masses
    B = 9.0 * 195.353 #to convert into natural units. I need the factor 195.353 and the resulting units for B turn into eV^2. See tesla_conversion.nim
    L = 9.26 * 5.06773e+06 # resulting units in eV^-1
    return (1 / 10**9 * B * L / 2.0) ** 2 # divided by 10**9 to convert from GeV^-1 to eV^-1

# Here, dfTel and lerpTel are global variables that are defined independently of the function telescopeEff.
# So they are only done once, and then used inside the function if this is called.
# This csv file assumes parallel light. But in our case whe get light from the solar core, which is about 3 arcmin,
# so the efficiency is lower, like in the paper of 2013.
#dfTel = pd.read_csv("data/llnl_xray_telescope_cast_effective_area_parallel_light_DTU_thesis.csv")
#dfTel["Efficiency"] = dfTel["EffectiveArea[cm²]"] / areaBore * (8.174 / 10.15) # The numbers scale the parallel light effective area to the CAST JCAP 2015 effective area
dfTel = pd.read_csv("data/Jaime_data/2016_DEC_Final_CAST_XRT/EffectiveArea.txt", delim_whitespace = True)
dfTel["Efficiency"] = dfTel["Area(cm^2)"] / areaBore
#print(dfTel)
lerpTel = interp1d(dfTel["E(keV)"], dfTel["Efficiency"], bounds_error = False, fill_value = 0.0)
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

#maybe ranme the function to position_weight. Because it can be used to compute the weight at any position, not only candidates.
def candidate_weights(dataset, x, y): # in a separate file, load this code and make checks like in here: https://github.com/Vindaar/TimepixAnalysis/blob/master/Analysis/ingrid/mcmc_limit_calculation.nim#L4005-L4118
    # In particular using axion_image_weights from setupAxionImageInterpolator I should make sure the weights sum up to 1.
    return dataset.axion_image_weights([x, y])
    """
    #quick not-working example. Because the sum of all the results in the function candidates_weights over every position on my readout
    #should be 1, becasue the weight is the fraction of solar flux that arrives to that position.
    # chip edge and number is up to you.
   xs = np.linspace(-chip edge, chip edge, 100)
    ys = np.linspace(-chip edge, chip edge, 100
    sizeOfBlock = chipArea / (100 * 100)
    sum = 0
    for x in xs:
    for y in ys:
        sum += sizeOfBlock * candidate_weights(x, y) # candidate_weights returns counts / cm²
    assert(abs(sum - 1.0) < 1e-6)
    """

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
def signal(dataset, E, g_aγ4, x, y): # in counts keV⁻¹ cm^-2
    ## Returns the axion flux based on `g` and energy `E`
    #print("Solar axion flux = ", solarAxionFlux(E), ", area of bore = ", areaBore, ", conversion probability = ", conversionProbability(), ", telescope eff = ", telescopeEff(E), ", detector eff = ", detectorEff(dataset, E), ", software eff = ", softwareEff(dataset, E), "gag = ", g_aγ4, "candidate weights = ", candidate_weights(dataset, x, y))
    return solarAxionFlux(E) * dataset.total_time * areaBore * conversionProbability() * telescopeEff(E) * detectorEff(dataset, E) * softwareEff(dataset, E) * g_aγ4 * candidate_weights(dataset, x, y) #be careful because this area is that of the bore, not the spot on the readout


def background(dataset, E):
    ## For the unbinned approach, I need to compute an interpolation of energies and background.    
    ## Note: area of interest is the region on the chip, in which the signal is focused!
    ## This also allows us to see that the "closer" we cut to the expected axion signal on the
    ## detector, the less background we have compared to the *fixed* signal flux!
    #print("Background = ", (lerpBackground(E) * totalTime * readoutArea))
    #return (dataset.lerpBackground(E) * dataset.total_time * readoutArea ) # in counts keV⁻¹  #be careful because this area is that of the spot on the readout, not the bore
    return (dataset.lerpBackground(E) * dataset.total_time) # in counts keV⁻¹ cm^-2
    #return (np.interp(E, Energies, Background) * totalTime * readoutArea)

def totalBackground(dataset):
    energies = np.linspace(0.0,10.0,1000)
    backgrounds = background(dataset, energies) # gives an array of 1000 background rates in counts keV⁻¹ cm^-2
    #return sc.simpson(backgrounds, energies) #this was before, when I had the wrong units in background() because I was including there the readoutArea
    integrated_background = sc.simpson(backgrounds, energies)
    return integrated_background * readoutArea  # Multiply by the readout area to get total counts

def likelihood(dataset, g_aγ4) -> float:
    result = np.exp(-(totalSignal(dataset, g_aγ4)+totalBackground(dataset))) #e^(-(s_tot + b_tot))
    cEnergies = np.array(dataset.candidates["Es"]) # this will have to be modified when we have the position
    candidate_pos_x = np.array(dataset.candidates["xs"])
    candidate_pos_y = np.array(dataset.candidates["ys"])
    idx = 0
    for candidate in range(len(cEnergies)): # column names in df are xs, ys, Es
        E = cEnergies[candidate]
        x = candidate_pos_x[candidate]
        y = candidate_pos_y[candidate]
        s = signal(dataset, E, g_aγ4, x, y)
        b = background(dataset, E)
        result *= s+b    
        #if result < 0.0:
        #    print("got less 0 ", result, " from ", s, " and ", b, " at ", E, " idx ", idx)
        #    quit()
        idx += 1
    return result

def likelihood2(dataset, g_aγ4) -> float: # Basti's version based on the division of Poissonian probabilities
    result = np.exp(-(totalSignal(dataset, g_aγ4))) #e^(-(s_tot ))
    #print("Total signal result = ", result)
    cEnergies = np.array(dataset.candidates["Es"]) # this will have to be modified when we have the position
    candidate_pos_x = np.array(dataset.candidates["xs"])
    candidate_pos_y = np.array(dataset.candidates["ys"])
    for candidate in range(len(cEnergies)): # column names in df are xs, ys, Es
        E = cEnergies[candidate]
        x = candidate_pos_x[candidate]
        y = candidate_pos_y[candidate]
        s = signal(dataset, E, g_aγ4, x, y)
        b = background(dataset, E)
        result *= (1.0 + s/b)
        #print("result ", result, " for candidate energy ", E, "background = ", b, " and signal = ", s)    
    return result

def logLikelihood(dataset, g_aγ4: float):
    return -np.log(likelihood(dataset, g_aγ4))  

def logLikelihood2(dataset, g_aγ4: float):
    return -np.log(likelihood2(dataset, g_aγ4))  

def logLikelihoodBasti(dataset, g_aγ4) -> float: # Basti's version based on the division of Poissonian probabilities
    result = (-(totalSignal(dataset, g_aγ4))) #(-(s_tot ))
    #print("Total signal result = ", result)
    cEnergies = np.array(dataset.candidates["Es"]) # this will have to be modified when we have the position
    candidate_pos_x = np.array(dataset.candidates["xs"])
    candidate_pos_y = np.array(dataset.candidates["ys"])
    for candidate in range(len(cEnergies)): # column names in df are xs, ys, Es
        E = cEnergies[candidate]
        x = candidate_pos_x[candidate]
        y = candidate_pos_y[candidate]
        s = signal(dataset, E, g_aγ4, x, y)
        b = background(dataset, E)
        result += np.log(1.0 + s/b)
        #print("result ", result, " for candidate energy ", E, "background = ", b, " and signal = ", s)
    return -result

def minusLogLikelihoodCris(dataset, g_aγ4) -> float: # Basti's version based on the division of Poissonian probabilities
    result = (totalSignal(dataset, g_aγ4) + totalBackground(dataset)) #((s_tot + b_tot ))
    #print("Total signal result = ", result)
    cEnergies = np.array(dataset.candidates["Es"]) # this will have to be modified when we have the position
    candidate_pos_x = np.array(dataset.candidates["xs"])
    candidate_pos_y = np.array(dataset.candidates["ys"])
    for candidate in range(len(cEnergies)): # column names in df are xs, ys, Es
        E = cEnergies[candidate]
        x = candidate_pos_x[candidate]
        y = candidate_pos_y[candidate]
        s = signal(dataset, E, g_aγ4, x, y)
        b = background(dataset, E)
        result -= np.log(s+b)
        #print("result ", result, " for candidate energy ", E, "background = ", b, " and signal = ", s)
    return result

def minusLogLikelihoodIgor(dataset, g_aγ4) -> float:
    result = (totalSignal(dataset, g_aγ4)+totalBackground(dataset)) #(s_tot + b_tot) are theoretical values >0
    cEnergies = np.array(dataset.candidates["Es"]) # this will have to be modified when we have the position
    candidate_pos_x = np.array(dataset.candidates["xs"])
    candidate_pos_y = np.array(dataset.candidates["ys"])
    idx = 0
    #print("before for loop: ", result)
    for candidate in range(len(cEnergies)): # column names in df are xs, ys, Es
        E = cEnergies[candidate]
        x = candidate_pos_x[candidate]
        y = candidate_pos_y[candidate]
        s = signal(dataset, E, g_aγ4, x, y) #this is s*f_s. The total signal would be s: ∫ s·f_s dx = totalSignal
        b = background(dataset, E)
        result -= np.log(s+b) #in the book it is +, but it is a typo.
        #if result < 0.0:
        #    print("got less 0 ", result, " from ", s, " and ", b, " at ", E, " idx ", idx)
        #    quit()
        idx += 1
    #print("After for loop: ", result)
    return result

def likelihoodIgor(dataset, g_aγ4) -> float:
    return np.exp(-minusLogLikelihoodIgor(dataset, g_aγ4))

def likelihoodCris(dataset, g_aγ4) -> float:
    return np.exp(-minusLogLikelihoodCris(dataset, g_aγ4))

def chi2(dataset, g_aγ4) -> float:
    return 2.0*(minusLogLikelihoodCris(dataset, g_aγ4))

def chi2basti(dataset, g_aγ4) -> float:
    return 2.0*(logLikelihood2(dataset, g_aγ4))

def chi2igor(dataset, g_aγ4) -> float:
    return 2.0*(minusLogLikelihoodIgor(dataset, g_aγ4))
    
# Function to find the minimum minus log-likelihood and corresponding g^4
def findMinimumMinusLogLikelihood(likelihoodFunction, dataset, g_aγ4_values):
    min_likelihood = float('inf')
    min_g_aγ4 = None    
    for g_aγ4 in g_aγ4_values:
        likelihood = likelihoodFunction(dataset, g_aγ4)
        if likelihood < min_likelihood:
            min_likelihood = likelihood
            min_g_aγ4 = g_aγ4
    #print(f"Minimum minus log-likelihood: {min_likelihood}")
    #print(f"Corresponding g^4 value: {min_g_aγ4}")
    print(f"\033[1;35;40m Minimum minus log-likelihood: {min_likelihood}\033[0m")
    print(f"\033[1;35;40m Corresponding g^4 value: {min_g_aγ4}\033[0m")

# Now we compute the combined likelihood, which is simply the product of the likelihood of each of the 3 datasets     
def totalLikelihood(dataset1, dataset2, dataset3, g_aγ4, likelihoodFunction=likelihood) -> float:
    likelihood_1 = likelihoodFunction(dataset1, g_aγ4)
    likelihood_2 = likelihoodFunction(dataset2, g_aγ4)
    likelihood_3 = likelihoodFunction(dataset3, g_aγ4)
    
    return likelihood_1 * likelihood_2 * likelihood_3

# Now we compute the combined log likelihood, which is simply the sum of the chi^2 of each of the 3 datasets     
def totalLogLikelihood(dataset1, dataset2, dataset3, g_aγ4, logLikelihoodFunction) -> float:
    logLikelihood_1 = logLikelihoodFunction(dataset1, g_aγ4)
    logLikelihood_2 = logLikelihoodFunction(dataset2, g_aγ4)
    logLikelihood_3 = logLikelihoodFunction(dataset3, g_aγ4)
    
    return logLikelihood_1 + logLikelihood_2 + logLikelihood_3

# Now we are going to write a function to compute a limit.
def limit(dataset, likelihoodFunction=likelihood):
    # Compute limit, CDF@95%
    # limit needs non logspace x & y data! (at least if computed in this simple way)
    g4Lin = np.linspace(0.0, 5e-40, 1000)
    likelihoodLin = likelihoodFunction(dataset, g4Lin)
    LCumSum = np.cumsum(likelihoodLin)          # cumulative sum
    LMax = LCumSum.max()               # maximum of the cumulative sum
    LCdf = LCumSum / LMax              # normalize to get (empirical) CDF
    limitIdx = bisect_left(LCdf, 0.95) # limit at 95% of the CDF

    return g4Lin[limitIdx]

# And now we write a fucntion to compute the combined limit
def totalLimit(dataset1, dataset2, dataset3, likelihoodFunction=likelihood):
    # Compute limit, CDF@95%
    # limit needs non logspace x & y data! (at least if computed in this simple way)
    g4Lin = np.linspace(0.0, 1e-40, 2000) # mayvbe I need more steps? Or a smaller g range?
    likelihoodLin = totalLikelihood(dataset1, dataset2, dataset3, g4Lin, likelihoodFunction)
    LCumSum = np.cumsum(likelihoodLin)          # cumulative sum
    LMax = LCumSum.max()               # maximum of the cumulative sum
    LCdf = LCumSum / LMax              # normalize to get (empirical) CDF
    limitIdx = bisect_left(LCdf, 0.95) # limit at 95% of the CDF

    return g4Lin[limitIdx]



def main():
    print(likelihood2(dataset_Xe, -1e-38))
    print(likelihood2(dataset_Xe, 0))
    print(likelihood2(dataset_Xe, 1e-38))
    #if True: quit()

    #Now we finally compute the limit for each dataset
    for dataset in [dataset_Ar1, dataset_Ar2, dataset_Xe]:
    #for dataset in [dataset_Xe]:
        print(f"Total number of expected background counts via integral of {dataset.name} =  {totalBackground(dataset)}" )
        lim = limit(dataset, likelihood)
        print(f"\033[1;35;40m Limit likelihood at : {pow(lim, 0.25)}\033[0m")

        print(f"Total number of expected background counts via integral of {dataset.name} =  {totalBackground(dataset)}" )
        lim = limit(dataset, likelihood2)
        print(f"\033[1;35;40m Limit likelihood 2 at : {pow(lim, 0.25)}\033[0m")

        print(f"Total number of expected background counts via integral of {dataset.name} =  {totalBackground(dataset)}" )
        lim = limit(dataset, likelihoodIgor)
        print(f"\033[1;35;40m Limit Igor at : {pow(lim, 0.25)}\033[0m")

        print(f"Total number of expected background counts via integral of {dataset.name} =  {totalBackground(dataset)}" )
        lim = limit(dataset, likelihoodCris)
        print(f"\033[1;35;40m Limit Cris at : {pow(lim, 0.25)}\033[0m")

        #g_aγs = np.linspace(-1e-40, 1e-40, 1000)
        #logL = [logLikelihood(g_aγ, dfCandidates) for g_aγ in g_aγs]
        #plt.plot(g_aγs, logL)
        #plt.xlabel(' $g_{aγ}^4$ (GeV$^{-4}$)')
        #plt.ylabel('-log L or Chi$^2/2$')
        #plt.savefig("ChiSquareg4_Nature_approach_unbinned.pdf")
        #plt.close()

    """
        g4Lin = np.linspace(0.0, 5e-40, 1000)
        likelihoodLin = [likelihoodIgor(dataset, g4) for g4 in g4Lin]
        plt.plot(g4Lin, likelihoodLin)
        plt.xlabel("Coupling constant (GeV$^-4$)")
        plt.ylabel("Likelihood")
        #plt.xscale('log')
        #plt.show()
        plt.axvline(x=lim, color='r')
        plt.savefig(f"energy_bins_likelihood_Igor_g4_unbinned_{dataset.name}.pdf")
        plt.close()

        g_aγs = np.linspace(-1e-41, 1e-41, 1000)
        logL2 = [logLikelihood(dataset, g_aγ) for g_aγ in g_aγs]
        plt.plot(g_aγs, logL2)
        plt.xlabel(' $g_{aγ}^4$ (GeV$^{-4}$)')
        plt.ylabel('-log L or Chi$^2/2$')
        plt.savefig(f"ChiSquareg4_logL2_approach_e41_{dataset.name}.pdf")
        plt.close()

        g_aγs = np.linspace(-1e-41, 1e-41, 1000)
        logL2 = [minusLogLikelihoodIgor(dataset, g_aγ) for g_aγ in g_aγs]
        plt.plot(g_aγs, logL2)
        plt.xlabel(' $g_{aγ}^4$ (GeV$^{-4}$)')
        plt.ylabel('-log L or Chi$^2/2$')
        plt.savefig(f"ChiSquareg4_Igor_approach_e41_{dataset.name}.pdf")
        plt.close()
    """
    # And now we compute the total likelihood and limit
    totalLim = totalLimit(dataset_Ar1, dataset_Ar2, dataset_Xe, likelihood)
    print(f"\033[1;35;40m Combined limit likelihood at : {pow(totalLim, 0.25)}\033[0m")

    totalLim = totalLimit(dataset_Ar1, dataset_Ar2, dataset_Xe, likelihood2)
    print(f"\033[1;35;40m Combined limit likelihood2 at : {pow(totalLim, 0.25)}\033[0m")

    totalLim = totalLimit(dataset_Ar1, dataset_Ar2, dataset_Xe, likelihoodIgor)
    print(f"\033[1;35;40m Combined limit Igor  at : {pow(totalLim, 0.25)}\033[0m")

    totalLim = totalLimit(dataset_Ar1, dataset_Ar2, dataset_Xe, likelihoodCris)
    print(f"\033[1;35;40m Combined limit Cris  at : {pow(totalLim, 0.25)}\033[0m")

    # first two quick plots
    g4Lin = np.linspace(-2e-40, 1e-40, 1000)
    likelihoodLin = totalLikelihood(dataset_Ar1, dataset_Ar2, dataset_Xe, g4Lin, likelihood)
    plt.plot(g4Lin, likelihoodLin)
    plt.xlabel("Coupling constant (GeV$^{-4}$)")
    plt.ylabel("Likelihood")
    plt.show()
    plt.axvline(x=totalLim, color='r')
    plt.savefig(f"energy_bins_likelihood_g4_unbinned_all_datasets_optimized.pdf")
    plt.close()

    g_aγs = np.linspace(-3.0e-40,1.0e-40, 1000)
    logL2 = totalLogLikelihood(dataset_Ar1, dataset_Ar2, dataset_Xe, g_aγs, chi2)
    plt.plot(g_aγs, logL2)
    plt.xlabel(' $g_{aγ}^4$ (GeV$^{-4}$)')
    plt.ylabel('-2log L or Chi$^2$')
    plt.show()
    plt.savefig(f"ChiSquareg4_Cris_approach_all_datasets_zoom_optimized.pdf")
    plt.close()

    # other plots
    g4Lin = np.linspace(0.0, 1e-40, 1000)
    likelihoodLin = totalLikelihood(dataset_Ar1, dataset_Ar2, dataset_Xe, g4Lin, likelihood)
    plt.plot(g4Lin, likelihoodLin)
    plt.xlabel("Coupling constant (GeV$^-4$)")
    plt.ylabel("Likelihood")
    #plt.xscale('log')
    #plt.show()
    plt.axvline(x=totalLim, color='r')
    plt.savefig(f"energy_bins_likelihood_g4_unbinned_all_datasets.pdf")
    plt.close()

    g4Lin = np.linspace(0.0, 1e-40, 1000)
    likelihoodLin = totalLikelihood(dataset_Ar1, dataset_Ar2, dataset_Xe, g4Lin, likelihood2)
    plt.plot(g4Lin, likelihoodLin)
    plt.xlabel("Coupling constant (GeV$^-4$)")
    plt.ylabel("Likelihood")
    #plt.xscale('log')
    #plt.show()
    plt.axvline(x=totalLim, color='r')
    plt.savefig(f"energy_bins_likelihood2_g4_unbinned_all_datasets.pdf")
    plt.close()

    g4Lin = np.linspace(0.0, 1e-40, 1000)
    likelihoodLin = totalLikelihood(dataset_Ar1, dataset_Ar2, dataset_Xe, g4Lin, likelihoodIgor)
    plt.plot(g4Lin, likelihoodLin)
    plt.xlabel("Coupling constant (GeV$^-4$)")
    plt.ylabel("Likelihood")
    plt.axvline(x=totalLim, color='r')
    plt.savefig(f"energy_bins_likelihoodIgor_g4_unbinned_all_datasets.pdf")
    plt.close()
    """
    g_aγs = np.linspace(-3.0e-40, 3.0e-40, 1000)
    logL2 = [totalLogLikelihood(dataset_Ar1, dataset_Ar2, dataset_Xe, g_aγ, minusLogLikelihoodIgor) for g_aγ in g_aγs]
    plt.plot(g_aγs, logL2)
    plt.xlabel(' $g_{aγ}^4$ (GeV$^{-4}$)')
    plt.ylabel('-log L or Chi$^2/2$')
    plt.savefig(f"ChiSquareg4_Igor_approach_all_datasets_zoom.pdf")
    plt.show()
    plt.close()
    #plt.savefig(f"ChiSquareg4_Nature_approach_all_datasets.pdf")
    findMinimumMinusLogLikelihood(minusLogLikelihoodIgor, dataset, g_aγs)
    """
    g4Lin = np.linspace(0.0, 1e-40, 1000)
    likelihoodLin = totalLikelihood(dataset_Ar1, dataset_Ar2, dataset_Xe, g4Lin, likelihoodCris)
    plt.plot(g4Lin, likelihoodLin)
    plt.xlabel("Coupling constant (GeV$^-4$)")
    plt.ylabel("Likelihood")
    plt.axvline(x=totalLim, color='r')
    plt.savefig(f"energy_bins_likelihoodCris_g4_unbinned_all_datasets.pdf")
    plt.close()
    """
    g_aγs = np.linspace(-3.0e-40, 3.0e-40, 1000)
    logL2 = [totalLogLikelihood(dataset_Ar1, dataset_Ar2, dataset_Xe, g_aγ, minusLogLikelihoodCris) for g_aγ in g_aγs]
    plt.plot(g_aγs, logL2)
    plt.xlabel(' $g_{aγ}^4$ (GeV$^{-4}$)')
    plt.ylabel('-log L or Chi$^2/2$')
    #plt.show()
    plt.savefig(f"ChiSquareg4_Cris_approach_all_datasets_zoom.pdf")
    plt.close()
    #plt.savefig(f"ChiSquareg4_Nature_approach_all_datasets.pdf")
    #findMinimumMinusLogLikelihood(likelihoodCris, dataset, g_aγs)
    """
    g_aγs = np.linspace(-3.0e-40, 3.0e-40, 1000)
    logL2 = totalLogLikelihood(dataset_Ar1, dataset_Ar2, dataset_Xe, g_aγs, chi2)# 2*minusLogLikelihoodCris
    plt.plot(g_aγs, logL2)
    plt.xlabel(' $g_{aγ}^4$ (GeV$^{-4}$)')
    plt.ylabel('-2log L or Chi$^2$')
    #plt.show()
    plt.savefig(f"ChiSquareg4_Cris_approach_all_datasets_zoom.pdf")
    plt.close()
    #plt.savefig(f"ChiSquareg4_Nature_approach_all_datasets.pdf")
    findMinimumMinusLogLikelihood(chi2, dataset, g_aγs)

    """
    g_aγs = np.linspace(-1.0e-40, 1.0e-40, 1000)
    logL2 = [totalLikelihood(dataset_Ar1, dataset_Ar2, dataset_Xe, g_aγ, logLikelihood) for g_aγ in g_aγs]
    plt.plot(g_aγs, logL2)
    plt.xlabel(' $g_{aγ}^4$ (GeV$^{-4}$)')
    plt.ylabel('-log L or Chi$^2/2$')
    plt.savefig(f"ChiSquareg4_likelihood_approach_all_datasets_e40.pdf")
    plt.close()
    """
    """
    g_aγs = np.linspace(-3.0e-40, 3.0e-40, 1000)
    logL2 = [totalLogLikelihood(dataset_Ar1, dataset_Ar2, dataset_Xe, g_aγ, logLikelihood2) for g_aγ in g_aγs]
    plt.plot(g_aγs, logL2)
    plt.xlabel(' $g_{aγ}^4$ (GeV$^{-4}$)')
    plt.ylabel('-log L or Chi$^2/2$')
    #plt.show()
    plt.savefig(f"ChiSquareg4_likelihood2_approach_all_datasets_zoom.pdf")
    plt.close()
    """
    g_aγs = np.linspace(-3.0e-40, 3.0e-40, 1000)
    logL2 = totalLogLikelihood(dataset_Ar1, dataset_Ar2, dataset_Xe, g_aγs, chi2basti)
    plt.plot(g_aγs, logL2)
    plt.xlabel(' $g_{aγ}^4$ (GeV$^{-4}$)')
    plt.ylabel('-2log L or Chi$^2$')
    #plt.show()
    plt.savefig(f"ChiSquareg4_chi2basti_approach_all_datasets_zoom.pdf")
    plt.close()

    g_aγs = np.linspace(-3.0e-40, 3.0e-40, 1000)
    logL2 = totalLogLikelihood(dataset_Ar1, dataset_Ar2, dataset_Xe, g_aγs, chi2igor)
    plt.plot(g_aγs, logL2)
    plt.xlabel(' $g_{aγ}^4$ (GeV$^{-4}$)')
    plt.ylabel('-2log L or Chi$^2$')
    #plt.show()
    plt.savefig(f"ChiSquareg4_chi2igor_approach_all_datasets_zoom.pdf")
    plt.close()
    """
    g_aγs = np.linspace(-1.0e-39, 1.0e-39, 1000)
    logL2 = [totalLikelihood(dataset_Ar1, dataset_Ar2, dataset_Xe, g_aγ, logLikelihoodBasti) for g_aγ in g_aγs]
    plt.plot(g_aγs, logL2)
    plt.xlabel(' $g_{aγ}^4$ (GeV$^{-4}$)')
    plt.ylabel('-log L or Chi$^2/2$')
    plt.savefig(f"ChiSquareg4_likelihoodBasti_approach_all_datasets_e39.pdf")
    plt.close()
    """

    """
    # I plot all the likelihoods together

    # Combining the plotting commands to make a single plot with all the curves
    # Define the datasets
    datasets = [dataset_Ar1, dataset_Ar2, dataset_Xe]
    # Plotting the likelihood for each individual dataset
    for dataset in datasets:
        g4Lin_individual = np.linspace(0.0, 5e-40, 1000)
        likelihoodLin_individual = [likelihood2(dataset, g4) for g4 in g4Lin_individual]
        plt.plot(g4Lin_individual, likelihoodLin_individual, label=f"{dataset.name}")
        lim_individual = limit(dataset, likelihood2)
        plt.axvline(x=lim_individual, linestyle='--', color='grey')
    # Plotting the total likelihood for all datasets combined
    g4Lin_total = np.linspace(0.0, 5e-40, 1000)
    likelihoodLin_total = [totalLikelihood(dataset_Ar1, dataset_Ar2, dataset_Xe, g4, likelihood2) for g4 in g4Lin_total]
    plt.plot(g4Lin_total, likelihoodLin_total, label='All Datasets Combined', color='black', linewidth=3)
    plt.axvline(x=totalLim, color='r', linewidth=2)

    # Adding labels, title, legend, and saving the plot
    plt.xlabel("Coupling constant (GeV$^-4$)")
    plt.ylabel("Likelihood")
    plt.title("Likelihood vs. Coupling Constant for Different Datasets using likelihood2")
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig("combined_likelihood_plot.pdf")
    #plt.show()
    """

    """
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

    # Plotting both fluxes in one figure
    plt.figure(figsize=(10, 6))
    plt.plot(E, solarAxionFlux_values*(g_aγ**2), label='Solar axion flux for $g_{a\gamma}·10^{-10}$')
    plt.plot(E, solarAxionFlux2019_values*(g_aγ**2), label='Solar axion flux 2019 for $g_{a\gamma}·10^{-10}$', linestyle='--')
    # Labels, title, and legend
    plt.xlabel('Energy (keV)')
    plt.ylabel('Solar axion flux (keV⁻¹ cm⁻² s⁻¹)')
    plt.title('Comparison of Solar Axion Fluxes')
    plt.legend()
    # Save and show the plot
    plt.tight_layout()
    plt.savefig("CombinedSolarAxionFlux.pdf")
    plt.show()

    # Telescope efficiency
    plt.plot(dfTel["E(keV)"], dfTel["Efficiency"])
    plt.xlabel("Energy (keV)")
    plt.ylabel("Efficiency")
    plt.savefig("telescopeEfficiency.pdf")
    plt.close()
    """ 

    """ 
    # My signal, which is the expected solar axion flux corrected by my efficiencies
    g_aγ4 = 1e-40
    signal_values = [signal(dataset, e, g_aγ4, x, y) for e in E]
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


    """
if __name__ == "__main__":
    main()

