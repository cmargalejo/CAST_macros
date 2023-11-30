from unbinned_likelihood_CAST_g4_datasets import *
import numpy as np
import scipy.integrate as sc

g_aγ4 = (6e-11)**4
dataset = dataset_Ar1

# Evaluate totalSignal
total_signal_value = totalSignal(dataset, g_aγ4)
print("Total signal is: ", total_signal_value, " counts")

total_background = totalBackground(dataset) # counts
backgroundTime = 2476
trackingTime = 130.367
backgroundExpCounts = 530
ratio = trackingTime / backgroundTime
proportional_bck = backgroundExpCounts * ratio # 530 is the number of counts during 2476 h of no tracking. I hace 130.367 h of tracking, so the proporiton is 19.
print("Total background during tracking time is: ", total_background, " counts, which should match the known number of background counts, i.e.,", proportional_bck ," counts.")

# Sum up the individual signals
Es = np.linspace(0.0, 10.0, 100)  # Energy steps (was 12, the number of bins in the bck histo)
energy_bin_width = 10.0/100.0 # (was 10.0 / 12.0)

# Define the number of rectangles and the total area they cover
n = 100 #40  # Number of rectangles along one dimension in the 20 mm side square
area_per_pixel = 4.0 / n**2  # Area of each small rectangle (pixel) in cm2

# Define the range for the coordinates (assuming a square detector area)
x_min, x_max = -10, 10  # Coordinate ranges in cm
y_min, y_max = -10, 10

# Calculate the number of steps and the size of each step for x and y
n_steps = n #int((x_max - x_min) / np.sqrt(area_per_pixel))
print("steps: ", n_steps)
x_steps = np.linspace(x_min, x_max, n_steps)
y_steps = np.linspace(y_min, y_max, n_steps)

cws = 0.0
for x in x_steps:
    for y in y_steps:
        cws += candidate_weights(dataset, x, y) * area_per_pixel
print("Candidate weights sum : ", cws)

individual_signals_sum = 0 # in c/keV/cm2
for E in Es:
    for x in x_steps:
        for y in y_steps:
            s = signal(dataset, E, g_aγ4, x, y)
            #print("Signal {}, area {}, ΔE {}".format(s, area_per_pixel, energy_bin_width))
            individual_signals_sum += s * area_per_pixel * energy_bin_width # in counts
#individual_signals_sum = np.sum(signal(dataset, Es, g_aγ4, x_steps, y_steps) * area_per_pixel * energy_bin_width) # in counts
            #units are counts keV⁻¹ cm^-2, so that's why I multiply by area and ΔE, to make sure the units are comparable with those of totalSignal

bs = 0.0
for E in Es:
    bs += background(dataset, E) * energy_bin_width * np.pi
print("Background integrated to total background time: ", bs, " vs expected: ", backgroundExpCounts * ratio)
# Compare the two values
"""
if np.isclose(total_signal_value, individual_signals_sum):
    print("The sum of individual signals is approximately equal to the total signal.")
    print(f"Total signal: {total_signal_value} counts, Sum of individual signals: {individual_signals_sum} counts")
else:
    print("There is a discrepancy between the sum of individual signals and the total signal.")
    print(f"Total signal: {total_signal_value} counts, Sum of individual signals: {individual_signals_sum} counts")
    """
print(f"Total signal: {total_signal_value} counts, Sum of individual signals: {individual_signals_sum} counts")

# add check that if I integrate the background I do get the expected number of counts, i.e. 530 counts according to my excel sheet
"""
Central bin value (keV)	NEW Background per energy range R=10mm	Error  (c/keV/cm^2/s)
0.5	1	1	3.57E-08	3.57E-08
1.5	18	4	6.43E-07	1.52E-07
2.5	74	9	2.64E-06	3.07E-07
3.5	53	7	1.89E-06	2.60E-07
4.5	38	6	1.36E-06	2.20E-07
5.5	36	6	1.29E-06	2.14E-07
6.5	35	6	1.25E-06	2.11E-07
7.5	63	8	2.25E-06	2.83E-07
8.5	127	11	4.54E-06	4.02E-07
9.5	46	7	1.64E-06	2.42E-07
10.5	21	5	7.50E-07	1.64E-07
11.5	18	4	6.43E-07	1.52E-07
530	23	1.58E-06	6.85E-08
"""

# also check that background and totalBackground give 530. In the background case go bin by bin and sum them up.

