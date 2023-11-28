#from unbinned_likelihood_CAST_g4_datasets import *
import pickle
import numpy as np
import matplotlib.pyplot as plt

results = None


#np.save("projected_sensitivity.dat", np.array(results))
with open("projected_sensitivity.dat", 'rb') as file:
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


    plt.plot(times, medians)
    plt.plot(times, mins)
    plt.plot(times, maxs)
    plt.plot(times, noCands)
    plt.legend()
    plt.xlabel("Time")
    plt.ylabel("Limit")
    plt.show()
    #plt.savefig(f"tests.pdf")
    #plt.close()

