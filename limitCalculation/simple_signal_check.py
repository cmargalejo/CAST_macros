from unbinned_likelihood_CAST_g4_datasets import *

g_aγ4 = (6e-11)**4
dataset = dataset_Ar1

# Evaluate totalSignal
total_signal_value = totalSignal(dataset, g_aγ4) # Flux integrated to total time, energy and area in counts of X-rays.
print("Total signal integrated to total time, energy and area in counts of X-rays is: ", total_signal_value, " counts in the 0 to 10 keV range.")

s3keV = signal(dataset, 3.0, g_aγ4, 0, 0) 
print("Signal at (0,0) for 3 keV: ", s3keV, " c/keV/cm2 (not per second!, i.e., during the full exposure time)") # in counts keV⁻¹ cm^-2

#In the inner 20 mm side square I have 40*40=1600 "pixels"
s = signal(dataset, 3.0, g_aγ4, 0, 0)*4.0 # 4 cm^2 is the area of my 2x2 cm rectangle
print("Signal at (0,0) for 3 keV if applied to an area of 4 cm^2: ", s, " counts / keV")

# now I do it for several squares to make it more realistic.
# i need to define the coordinates of the center of each square
n = 1000 #number of rectangles in each dimension
area = 4.0/n**2 #area per rectangle or pixel
xMin = -10
xMax = 10
xSize = 2.0 / n
coord = np.linspace(xMin + xSize / 2.0, xMax - xSize / 2.0, n)
sum = 0.0
for x in coord:
    for y in coord:
        sum += signal(dataset, 3.0, g_aγ4 , x, y) * area * 10 #keV maybe 10 keV
print("Sum = ", sum, " counts in 0 to 10 keV range") # upper limit


# Compare the two values
if np.isclose(total_signal_value, sum):
    print("The sum of individual signals is approximately equal to the total signal.")
    print(f"Total signal: {total_signal_value}, counts. Upper limit of the sum of individual signals: {sum} counts")
else:
    print("There is a discrepancy between the sum of individual signals and the total signal.")
    print(f"Total signal: {total_signal_value}, counts. Upper limit of the sum of individual signals: {sum} counts")

"""
 # Total signal
g_aγ4 = np.linspace(1e-52, 1e-40, 1000)
totalSignal_values = [totalSignal(g4) for g4 in g_aγ4]
plt.plot(g_aγ4, totalSignal_values)
plt.xlabel('$g_{aγ}^4$ (GeV$^{-4}$)')
plt.ylabel('Total signal (counts)')
plt.show()
"""