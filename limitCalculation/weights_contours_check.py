from unbinned_likelihood_CAST_g4_datasets import *


def setupAxionImageInterpolatorCircles(filename, convolved_resolution): #I should check individually that this sums up to 1, by computing the
    # Read the text file using Pandas
    df = pd.read_csv(filename, skiprows = 1, delim_whitespace = True, names = ["x", "y", "z", "zMean"])
    df["x"] = df["x"] - 31.114 # these are the values of the centroid that match it with the position of the X-ray finger centroid
    df["y"] = df["y"] - 30.283 # these are the values of the centroid that match it with the position of the X-ray finger centroid 30.283
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
    #zSum = np.sum(df["z"])
    zSums = np.zeros(5) # the number of circles+rings
    areas = np.zeros(5) # in mm^2
    weights = np.zeros(5) # rate as a percentage
    normW = np.zeros(5) # rate / area in cm^2

    # Define the circles and rings. I start by defining the radii
    # then r^2 = x^2 + y^2
    # radii in cm
    r0 = 0
    r1 = 0.1 # radius 1 = 1 mm
    r2 = 0.2 
    r3 = 0.3
    r4 = 0.4
    r5 = 1
    for x in range(size):
        for y in range(size):
            index = x * size + y # it computes the index we are looking at.
            z = df["z"][index] # and then it finds the z value at that index.
            # Now I want to find the (x,y) position for each z
            xpos = df["x"][index]
            ypos = df["y"][index]
            radius = np.sqrt((xpos*xpos) + (ypos*ypos)) / 10 # in cm
            #print("X position = ", xpos, ", Y position = ", ypos, ", radius = ", radius)
            if (radius >= r0 and radius < r1):
                zSums[0] = zSums[0] + z
                areas[0] = np.pi * (r1**2 - r0**2)
            elif (radius >= r1 and radius < r2):
                zSums[1] = zSums[1] + z
                areas[1] = np.pi * (r2**2 - r1**2)
            elif (radius >= r2 and radius < r3):
                zSums[2] = zSums[2] + z
                areas[2] = np.pi * (r3**2 - r2**2)
            elif (radius >= r3 and radius < r4):
                zSums[3] = zSums[3] + z
                areas[3] = np.pi * (r4**2 - r3**2)
            elif (radius >= r4 and radius < r5):
                zSums[4] = zSums[4] + z
                areas[4] = np.pi * (r5**2 - r4**2)
    #print("zSums = ", zSums)
    #print("areas = ", areas)
    zSum = np.sum(zSums)
    #print("zSum = ", zSum)

    #plt.bar([r1, r2, r3, r4, r5],zSums)
    #plt.show()

    weights = zSums / zSum
    #print("weights = ", weights)
    wSum = np.sum(weights)
    #print("Weights sum = ", wSum)
    normW = zSums / zSum / areas
    #print("normalised weights = ", normW)
    normWSum = np.sum(normW)
    #print("Normalised weights sum = ", normWSum)

    def findWeight(ar): # ar = [x, y]
        radius = np.sqrt((ar[0]*ar[0] + ar[1]*ar[1])) / 10.0 # in cm
        if radius <= r1:
            return normW[0]
        elif radius <= r2:
            return normW[1]
        elif radius <= r3:
            return normW[2]
        elif radius <= r4:
            return normW[3]
        elif radius <= r5:
            return normW[4]

    return findWeight
    
for dataset in [dataset_Ar1, dataset_Ar2, dataset_Xe]: # [dataset_Ar1, dataset_Ar2, dataset_Xe]
        dataset.axion_image_weights = setupAxionImageInterpolatorCircles(axion_image_filename, 500)
        lim = limit(dataset, likelihood)
        print(f"\033[1;35;40m Limit likelihood at : {pow(lim, 0.25)}\033[0m")

totalLim = totalLimit(dataset_Ar1, dataset_Ar2, dataset_Xe, likelihood)
print(f"\033[1;35;40m Combined limit likelihood at : {pow(totalLim, 0.25)}\033[0m")

# first two quick plots
g4Lin = np.linspace(0, 1e-40, 1000)
likelihoodLin = totalLikelihood(dataset_Ar1, dataset_Ar2, dataset_Xe, g4Lin, likelihood)
plt.plot(g4Lin, likelihoodLin)
plt.xlabel("Coupling constant (GeV$^{-4}$)")
plt.ylabel("Likelihood")
plt.axvline(x=totalLim, color='r')
plt.show()
#plt.savefig(f"energy_bins_likelihood_g4_unbinned_all_datasets_optimized.pdf")
plt.close()

g_aγs = np.linspace(-3.0e-40,1.0e-40, 1000)
logL2 = totalLogLikelihood(dataset_Ar1, dataset_Ar2, dataset_Xe, g_aγs, chi2)
plt.plot(g_aγs, logL2)
plt.xlabel(' $g_{aγ}^4$ (GeV$^{-4}$)')
plt.ylabel('-2log L or Chi$^2$')
plt.show()
#plt.savefig(f"ChiSquareg4_Cris_approach_all_datasets_zoom_optimized.pdf")
plt.close()
