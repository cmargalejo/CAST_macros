#from unbinned_likelihood_CAST_g4_datasets import *
import pickle
import numpy as np
import matplotlib.pyplot as plt

results = None
#np.save("projected_sensitivity.dat", np.array(results))
#with open("projected_sensitivity.dat", 'rb') as file:
with open("projected_sensitivity_from_Xe.dat", 'rb') as file: # projected_sensitivity_Xe.dat
    results = pickle.load(file)

for dataset in results:
    medians = []
    mins = []
    maxs = []
    noCands = []
    times = []
    for r in dataset:
        #print("iteration number: ", r[0], ", 0 candidates limit =", r[1], ", MC limits = ", r[2])
        print("iteration number: ", r[0], ", 0 candidates limit =", r[2])
        limits = r[3]
        medians.append(np.median(limits))
        mins.append(np.min(limits))
        maxs.append(np.max(limits))
        noCands.append(r[2])
        times.append(r[1])
    print(times)

    times = np.array(times)
    #times = times / 3600.0 / 3.0 / (365.25 / 12.0) # divide by 3660 to have hours, divide by 3 because we have 3 hours per day. Then divide by number of days per month to have it in months
    days = times / 3600.0 / 3.0 # 3 hours = 10800 s = 1 day
    months = days / (365.25 / 12.0)

    # Calculate the fourth roots
    mins_root = np.array(mins) ** 0.25
    maxs_root = np.array(maxs) ** 0.25
    medians_root = np.array(medians) ** 0.25
    noCands_root = np.array(noCands) ** 0.25
  
    # Plotting
    plt.plot(months, medians_root, label='Median', color='black', linewidth=2)  # Solid, thick, black line
    plt.plot(months, mins_root, label='Min', color='black', linestyle='dashed')  # Black dashed line
    plt.plot(months, maxs_root, label='Max', color='black', linestyle='dashed')  # Black dashed line
    plt.plot(months, noCands_root, label='No Candidates', color='red', linestyle='dashed')  # Red dashed line
    #plt.axhline(y=8.53e-11, color='blue', linestyle='--', label=f'Argon observed limit')
    #plt.axhline(y=8.68e-11, color='gray', linestyle='--', label=f'Argon expected limit')
    #plt.axvline(x=1.42, color='blue', linestyle='--', label=f'Argon current exposure')
    plt.axhline(y=7.75e-11, color='blue', linestyle='--', label=f'Xenon observed limit')
    plt.axhline(y=7.78e-11, color='gray', linestyle='--', label=f'Xenon expected limit')
    plt.axvline(x=1.75, color='blue', linestyle='--', label=f'Xenon current exposure')

    plt.legend()
    plt.xlabel("Exposure (months)")
    plt.ylabel("Sensitivity (95% UL on gag)")
    plt.grid(True)
    plt.show()
    #plt.savefig(f"tests.pdf")
    #plt.close()

