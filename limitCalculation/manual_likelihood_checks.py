from unbinned_likelihood_CAST_g4_datasets import *

dataset = dataset_Ar1
c = candidate_weights(dataset_Xe,0,0)
print(c)
#cs = [()]
candidate_pos_x = np.array(dataset.candidates["xs"])
candidate_pos_y = np.array(dataset.candidates["ys"])
cEnergies = np.array(dataset.candidates["Es"]) 

df = datase t.candidates
#df = df[(df["xs"] > -2 and df["xs"] < 2 and df["ys"] > -5 and df["ys"] < 5)] 
df = df[(df["xs"] > -2)]
df = df[(df["xs"] < 2)]
df = df[(df["ys"] > -5)]
df = df[(df["ys"] < 5)] 
print(df)
plt.scatter(df["xs"],df["ys"])
plt.show()

"""
for candidate in range(len(cEnergies)): # column names in df are xs, ys, Es
    E = cEnergies[candidate]
    x = candidate_pos_x[candidate]
    y = candidate_pos_y[candidate]
    #s = signal(dataset, E, g_aγ4, x, y)
    #b = background(dataset, E)
    if (x > -2 and x < 2 and y > -5 and y < 5):
        print("Candidate E: ", E, ", x pos = ", x, ", y pos = ", y)
        #plt.scatter(x,y)
        #plt.show()
    #result *= s+b    
    #idx += 1
    #return result
"""